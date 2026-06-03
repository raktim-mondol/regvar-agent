"""Variant I/O and ranking.

Keeps coordinate conventions in one place and documented:

    * VCF / our candidate TSV are 1-based, the same convention AlphaGenome's
      `genome.Variant.position` expects, so no shift is applied to positions.
    * BED-style intervals are 0-based half-open; if we ever ingest BED we
      convert on the way in (see `bed_to_one_based`).

Getting this wrong is the single most common silent error in variant-effect
work, which is why it is isolated and unit-tested rather than scattered.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

_VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "M"}

MAX_VARIANTS: int = int(os.environ.get("REGVAR_MAX_VARIANTS", "500"))


def validate_chromosome(name: str) -> str:
    """Normalize and validate a chromosome name.

    Accepts ``chr1``–``chr22``, ``chrX``, ``chrY``, ``chrM`` (and the same
    without the ``chr`` prefix). Returns the ``chr``-prefixed form.

    Raises ``ValueError`` for anything else.
    """
    stripped = name.strip()
    core = stripped[3:] if stripped.lower().startswith("chr") else stripped
    core_upper = core.upper()
    if core_upper == "MT":
        core_upper = "M"
    if core_upper not in _VALID_CHROMS:
        raise ValueError(
            f"Invalid chromosome {name!r}. Expected chr1–chr22, chrX, chrY, or chrM."
        )
    return f"chr{core_upper}"


@dataclass(frozen=True)
class CandidateVariant:
    chromosome: str
    position: int          # 1-based
    ref: str
    alt: str
    region_id: str = ""    # e.g. enhancer / peak id from ATAC or ChIP data
    note: str = ""

    @property
    def vcf_id(self) -> str:
        return f"{self.chromosome}:{self.position}:{self.ref}>{self.alt}"

    def as_dict(self) -> dict:
        return {
            "chromosome": self.chromosome,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
        }


def bed_to_one_based(start: int) -> int:
    """BED start (0-based) -> 1-based coordinate."""
    return start + 1


def check_variant_count(
    variants: list[CandidateVariant], *, force: bool = False,
) -> list[CandidateVariant]:
    """Raise ``ValueError`` if *variants* exceeds ``MAX_VARIANTS``.

    Pass ``force=True`` to skip the check (e.g. when the user passes
    ``--force`` on the CLI).  Always returns the list unchanged so it can
    be used inline: ``variants = check_variant_count(variants)``.
    """
    if not force and len(variants) > MAX_VARIANTS:
        raise ValueError(
            f"Variant count ({len(variants)}) exceeds the safety limit of "
            f"{MAX_VARIANTS}. Pass --force to proceed anyway, or set "
            f"REGVAR_MAX_VARIANTS to raise the limit."
        )
    return variants


def read_candidates_tsv(path: str | Path) -> list[CandidateVariant]:
    """Read a tab-separated candidate file.

    Required columns: chromosome, position, ref, alt
    Optional columns: region_id, note
    """
    variants: list[CandidateVariant] = []
    with Path(path).open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"chromosome", "position", "ref", "alt"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for row_num, row in enumerate(reader, start=2):
            try:
                chrom = validate_chromosome(row["chromosome"])
            except ValueError as exc:
                raise ValueError(f"{path} line {row_num}: {exc}") from None
            variants.append(
                CandidateVariant(
                    chromosome=chrom,
                    position=int(row["position"]),
                    ref=row["ref"].strip().upper(),
                    alt=row["alt"].strip().upper(),
                    region_id=row.get("region_id", "").strip(),
                    note=row.get("note", "").strip(),
                )
            )
    return variants


def read_candidates_vcf(path: str | Path) -> list[CandidateVariant]:
    """Read a VCF or VCF.gz file and return candidate variants.

    Requires the optional ``cyvcf2`` dependency (``pip install regvar[vcf]``).

    Multi-allelic sites are expanded so each ALT allele becomes its own
    ``CandidateVariant``.  Chromosome names without a ``chr`` prefix get one
    added automatically, since AlphaGenome expects ``chr8`` not ``8``.
    The VCF ``ID`` field (e.g. an rsID) is used as ``region_id`` when it is
    not missing (``.``).
    """
    try:
        from cyvcf2 import VCF
    except ImportError as exc:
        raise ImportError(
            "cyvcf2 is required to read VCF files. "
            "Install with: pip install regvar[vcf]"
        ) from exc

    variants: list[CandidateVariant] = []
    for record in VCF(str(path)):
        chrom = validate_chromosome(record.CHROM)
        pos = record.POS  # 1-based, matching VCF convention
        ref = record.REF.upper()
        vcf_id = record.ID if record.ID and record.ID != "." else ""
        for alt_allele in record.ALT:
            variants.append(
                CandidateVariant(
                    chromosome=chrom,
                    position=pos,
                    ref=ref,
                    alt=alt_allele.upper(),
                    region_id=vcf_id,
                )
            )
    return variants


def write_scores_tsv(records: list[dict], path: str | Path) -> None:
    """Write a list of score records to a TSV file.

    Centralises the column-extraction + tab-join pattern used by the CLI,
    TUI, and agent run output.
    """
    if not records:
        return
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(records[0].keys())
    lines = ["\t".join(cols)]
    for rec in records:
        lines.append("\t".join(str(rec.get(c, "")) for c in cols))
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def rank_by_effect(tidy_scores, top_n: int = 20):
    """Given a tidy AlphaGenome score table, return the strongest predicted
    effects ranked by |quantile_score| (a calibrated, cross-assay-comparable
    measure). Falls back to |raw_score| if the quantile column is absent."""
    import pandas as pd

    df = tidy_scores.copy()
    score_col = "quantile_score" if "quantile_score" in df.columns else "raw_score"
    df["abs_effect"] = df[score_col].abs()
    keep = [c for c in (
        "query_variant", "variant_id", "assay", "output_type", "gene_name",
        "biosample_name", "raw_score", "quantile_score", "abs_effect",
    ) if c in df.columns]
    return (
        df[keep]
        .sort_values("abs_effect", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
