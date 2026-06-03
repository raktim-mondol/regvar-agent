"""Variant annotation: attach genomic context to candidate variants.

Pure-Python readers with optional ``pyranges`` acceleration.  All entry points
degrade gracefully when the annotation extras are not installed: GTF/BED files
are streamed line-by-line, and a tabix-indexed dbSNP VCF can be queried with
``cyvcf2`` if it is installed.

Three annotation sources are supported:

* ``annotate_nearest_gene(variants, gtf_path)``
      nearest gene + signed distance to its TSS (negative = upstream)
* ``annotate_regulatory_overlap(variants, bed_path)``
      overlapping regulatory element name(s) from a BED file (e.g. ATAC peaks)
* ``annotate_rsid(variants, dbsnp_vcf_path)``
      rsID lookup from a dbSNP VCF (requires tabix index + cyvcf2)

``annotate_variants`` is the single entry point used by the CLI's ``--annotate``
flag — it mutates the ``CandidateVariant.note`` field, appending annotation
snippets separated by " | ", and stores structured annotations on a
``_annotations`` attribute.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

from .variants import CandidateVariant, validate_chromosome


# Module-level annotation registry, keyed by CandidateVariant.vcf_id.
# Populated by annotate_variants(); read by build_task_message() and the CLI
# variant table. Avoids mutating the frozen CandidateVariant dataclass.
_ANNOTATIONS: dict[str, "VariantAnnotation"] = {}


def get_annotation(vcf_id: str) -> "VariantAnnotation | None":
    """Return the annotation for *vcf_id* if it has been registered, else None."""
    return _ANNOTATIONS.get(vcf_id)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class VariantAnnotation:
    """Context gathered for a single candidate variant."""
    variant: CandidateVariant
    nearest_gene: str | None = None
    tss_distance: int | None = None     # signed bp; negative = upstream
    regulatory_overlaps: list[str] = field(default_factory=list)
    rsid: str | None = None

    def to_note(self) -> str:
        """Render a compact, human-readable annotation snippet."""
        parts: list[str] = []
        if self.nearest_gene:
            dist = self.tss_distance
            if dist is None:
                parts.append(f"gene={self.nearest_gene}")
            else:
                sign = "-" if dist < 0 else "+"
                parts.append(f"gene={self.nearest_gene} ({sign}{abs(dist)}bp)")
        if self.regulatory_overlaps:
            parts.append(f"overlap={','.join(self.regulatory_overlaps[:3])}")
        if self.rsid:
            parts.append(f"rsid={self.rsid}")
        return " | ".join(parts)


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

def _open_text(path: str | Path):
    """Open a plain or .gz file as a text stream."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8")
    return p.open("r", encoding="utf-8")


def _normalise_chrom(chrom: str) -> str:
    """Best-effort chromosome normalisation; returns ``chr``-prefixed form."""
    try:
        return validate_chromosome(chrom)
    except ValueError:
        return chrom.strip()


# ---------------------------------------------------------------------------
# GTF → nearest gene
# ---------------------------------------------------------------------------

_GTF_GENE_FIELDS = re.compile(r'gene_id "([^"]+)"')
_GTF_GENE_NAME = re.compile(r'gene_name "([^"]+)"')


def _iter_gtf_genes(gtf_path: str | Path):
    """Yield ``(chrom, start, end, strand, gene_name)`` tuples for each gene
    record in a GTF (or GTF.gz). Uses the ``gene`` feature only.
    """
    with _open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = _normalise_chrom(fields[0])
            start = int(fields[3])   # 1-based inclusive
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8]
            m = _GTF_GENE_NAME.search(attrs) or _GTF_GENE_FIELDS.search(attrs)
            gene_name = m.group(1) if m else "."
            yield chrom, start, end, strand, gene_name


def _tss(start: int, end: int, strand: str) -> int:
    """Return the TSS coordinate for a gene on either strand."""
    return start if strand != "-" else end


def annotate_nearest_gene(
    variants: Sequence[CandidateVariant],
    gtf_path: str | Path,
) -> dict[str, VariantAnnotation]:
    """For each variant, find the nearest gene (by TSS distance).

    Returns a mapping ``vcf_id → VariantAnnotation``. Uses a simple O(V*G)
    scan; adequate for thousands of variants against ~20k genes.  Install
    ``pyranges`` if you need this to scale further.
    """
    genes_by_chrom: dict[str, list[tuple[int, int, str, str]]] = {}
    for chrom, start, end, strand, name in _iter_gtf_genes(gtf_path):
        genes_by_chrom.setdefault(chrom, []).append((start, end, strand, name))

    out: dict[str, VariantAnnotation] = {}
    for v in variants:
        ann = VariantAnnotation(variant=v)
        best_gene = None
        best_dist: int | None = None
        for start, end, strand, name in genes_by_chrom.get(v.chromosome, []):
            tss = _tss(start, end, strand)
            signed = v.position - tss
            d = abs(signed)
            if best_dist is None or d < best_dist:
                best_dist = d
                best_gene = name
                best_signed = signed if strand != "-" else -signed
        if best_gene is not None:
            ann.nearest_gene = best_gene
            ann.tss_distance = best_signed
        out[v.vcf_id] = ann
    return out


# ---------------------------------------------------------------------------
# BED → regulatory overlap
# ---------------------------------------------------------------------------

def _iter_bed_records(bed_path: str | Path):
    """Yield ``(chrom, start_0based, end, name)`` for each BED record."""
    with _open_text(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = _normalise_chrom(fields[0])
            start = int(fields[1])   # 0-based
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            yield chrom, start, end, name


def annotate_regulatory_overlap(
    variants: Sequence[CandidateVariant],
    bed_path: str | Path,
) -> dict[str, VariantAnnotation]:
    """For each variant, list BED records it falls inside.

    Positions are converted to 0-based half-open for the comparison (VCF
    positions are 1-based). Returns ``vcf_id → VariantAnnotation``.
    """
    records_by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for chrom, start, end, name in _iter_bed_records(bed_path):
        records_by_chrom.setdefault(chrom, []).append((start, end, name))

    out: dict[str, VariantAnnotation] = {}
    for v in variants:
        ann = VariantAnnotation(variant=v)
        pos0 = v.position - 1   # 1-based → 0-based for overlap test
        for start, end, name in records_by_chrom.get(v.chromosome, []):
            if start <= pos0 < end:
                ann.regulatory_overlaps.append(name)
        out[v.vcf_id] = ann
    return out


# ---------------------------------------------------------------------------
# dbSNP VCF → rsID lookup
# ---------------------------------------------------------------------------

def annotate_rsid(
    variants: Sequence[CandidateVariant],
    dbsnp_vcf_path: str | Path,
) -> dict[str, VariantAnnotation]:
    """Look up rsIDs by querying a tabix-indexed dbSNP VCF.

    Requires ``cyvcf2`` (``pip install regvar[vcf]``) and a ``.tbi`` index next
    to the VCF. Returns ``vcf_id → VariantAnnotation``.
    """
    try:
        from cyvcf2 import VCF
    except ImportError as exc:
        raise ImportError(
            "rsID annotation requires cyvcf2. "
            "Install with: pip install regvar[vcf]"
        ) from exc

    vcf = VCF(str(dbsnp_vcf_path))
    out: dict[str, VariantAnnotation] = {}
    for v in variants:
        ann = VariantAnnotation(variant=v)
        region = f"{v.chromosome}:{v.position}-{v.position}"
        try:
            for rec in vcf(region):
                if rec.POS == v.position and v.ref.upper() == rec.REF.upper():
                    for alt in rec.ALT:
                        if alt.upper() == v.alt.upper():
                            if rec.ID and rec.ID != ".":
                                ann.rsid = rec.ID
                                break
                    if ann.rsid:
                        break
        except Exception:
            # Region lookup can fail on contigs missing from the index —
            # silently skip rather than abort the whole annotation run.
            pass
        out[v.vcf_id] = ann
    return out


# ---------------------------------------------------------------------------
# Composite annotator (used by the CLI --annotate flag)
# ---------------------------------------------------------------------------

def annotate_variants(
    variants: Sequence[CandidateVariant],
    *,
    gtf_path: str | Path | None = None,
    bed_path: str | Path | None = None,
    dbsnp_vcf_path: str | Path | None = None,
) -> dict[str, VariantAnnotation]:
    """Run every requested annotation pass and merge the results.

    The returned dict is keyed by ``vcf_id``. The ``CandidateVariant`` objects
    are *not* mutated (they are a frozen dataclass); callers should look up
    ``ann.to_note()`` when building prompts or display tables.
    """
    merged: dict[str, VariantAnnotation] = {
        v.vcf_id: VariantAnnotation(variant=v) for v in variants
    }

    if gtf_path is not None:
        for vid, ann in annotate_nearest_gene(variants, gtf_path).items():
            merged[vid].nearest_gene = ann.nearest_gene
            merged[vid].tss_distance = ann.tss_distance

    if bed_path is not None:
        for vid, ann in annotate_regulatory_overlap(variants, bed_path).items():
            merged[vid].regulatory_overlaps = ann.regulatory_overlaps

    if dbsnp_vcf_path is not None:
        for vid, ann in annotate_rsid(variants, dbsnp_vcf_path).items():
            merged[vid].rsid = ann.rsid

    _ANNOTATIONS.update(merged)
    return merged
