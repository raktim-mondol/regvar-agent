"""Data model and parsers for regvar variant-effect reports.

Parses the TSV scores and markdown reports produced by ``regvar agent run``.

Public exports
--------------
ScoreRecord, CandidateReport  — data model
parse_scores                  — TSV score parser
extract_interpretations       — markdown interpretation extractor
assign_tiers                  — tier label scanner
synthesize_interpretation     — fallback text generator from scores
"""

from __future__ import annotations

import csv
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

from regvar.variants import CandidateVariant, read_candidates_tsv


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class ScoreRecord:
    """A single assay score for a variant."""
    assay: str
    gene_name: str
    biosample_name: str
    output_type: str
    raw_score: float
    quantile_score: float
    abs_effect: float


@dataclass
class CandidateReport:
    """All data needed to render a one-page variant report."""
    variant: CandidateVariant
    vcf_id: str
    rsid: str
    cytoband: str
    top_scores: list[ScoreRecord]
    max_effect: float
    interpretation: str
    tier: str
    validation_plan: str


# ---------------------------------------------------------------------------
# Score parsing
# ---------------------------------------------------------------------------

REQUIRED_SCORE_COLUMNS = [
    "query_variant", "gene_name", "assay", "output_type",
    "biosample_name", "raw_score", "quantile_score", "abs_effect",
]


def parse_scores(tsv_path: str | Path) -> dict[str, list[ScoreRecord]]:
    """Read scored variants TSV, group by variant, sort by |effect| desc.

    Returns
    -------
    dict mapping ``query_variant`` → list of ScoreRecord (sorted best-first).
    """
    path = Path(tsv_path)
    if not path.exists():
        raise FileNotFoundError(f"Scores file not found: {path}")

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(REQUIRED_SCORE_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing columns in {path}: {sorted(missing)}. "
                f"Expected: {', '.join(REQUIRED_SCORE_COLUMNS)}"
            )

        grouped: dict[str, list[ScoreRecord]] = {}
        for row in reader:
            vid = row["query_variant"].strip()
            sr = ScoreRecord(
                assay=row["assay"].strip(),
                gene_name=row["gene_name"].strip(),
                biosample_name=row["biosample_name"].strip(),
                output_type=row["output_type"].strip(),
                raw_score=float(row["raw_score"]),
                quantile_score=float(row["quantile_score"]),
                abs_effect=float(row["abs_effect"]),
            )
            grouped.setdefault(vid, []).append(sr)

    if not grouped:
        raise ValueError(f"No scored variants found in {path}")

    # Sort each group by abs_effect descending
    for scores in grouped.values():
        scores.sort(key=lambda s: s.abs_effect, reverse=True)

    return grouped


# ---------------------------------------------------------------------------
# Markdown extraction
# ---------------------------------------------------------------------------

_VARIANT_ID_PATTERN = re.compile(
    r"chr(\d{1,2}|[XYM]):(\d+):([ACGT]+)>([ACGT]+)", re.IGNORECASE
)
_RSID_PATTERN = re.compile(r"rs\d+")
_TIER_HEADING = re.compile(r"#+\s*(Tier\s*\d+)", re.IGNORECASE)


def extract_interpretations(md_path: str | Path) -> dict[str, str]:
    """Parse agent markdown report, extracting per-variant interpretation text.

    Strategy:
    1. Split markdown into blocks by variant ID mentions or section headings.
    2. For each variant found, grab the surrounding paragraph(s) as interpretation.
    3. Returns dict mapping variant ID string → interpretation paragraph(s).
    """
    path = Path(md_path)
    if not path.exists():
        raise FileNotFoundError(f"Report file not found: {path}")

    text = path.read_text()
    interpretations: dict[str, str] = {}

    # Find all variant IDs in the text
    variant_positions: list[tuple[int, str]] = []
    for m in _VARIANT_ID_PATTERN.finditer(text):
        vid = m.group(0)
        variant_positions.append((m.start(), vid))

    if not variant_positions:
        return interpretations

    # For each variant, extract text from its first mention to the next variant or
    # section break, up to ~500 chars.
    for i, (start, vid) in enumerate(variant_positions):
        if i + 1 < len(variant_positions):
            end = variant_positions[i + 1][0]
        else:
            end = min(start + 2000, len(text))

        chunk = text[start:end]

        # Try to extract the most relevant paragraph: look for the paragraph
        # containing the variant ID and the one after it.
        paragraphs = [p.strip() for p in chunk.split("\n\n") if p.strip()]
        relevant: list[str] = []
        for p in paragraphs:
            if len(p) > 30 and not p.startswith("|") and not p.startswith("#"):
                relevant.append(p)
            if len(" ".join(relevant)) > 600:
                break

        combined = " ".join(relevant[:3])  # up to 3 paragraphs
        if combined:
            if vid not in interpretations or len(combined) > len(interpretations[vid]):
                interpretations[vid] = combined

    return interpretations


def assign_tiers(md_path: str | Path) -> dict[str, str]:
    """Scan markdown for tier labels and map variant identifiers to tiers.

    Returns dict mapping variant ID (rsID or chr:pos format) → tier string.
    """
    path = Path(md_path)
    if not path.exists():
        return {}

    text = path.read_text()
    tiers: dict[str, str] = {}

    # Split by tier headings
    sections = _TIER_HEADING.split(text)
    current_tier = ""
    for part in sections:
        part_stripped = part.strip()
        if part_stripped.lower().startswith("tier"):
            current_tier = part_stripped
            continue
        if current_tier:
            # Find variant IDs in this tier section
            for m in _VARIANT_ID_PATTERN.finditer(part):
                vid = m.group(0)
                if vid not in tiers:
                    tiers[vid] = current_tier
            for m in _RSID_PATTERN.finditer(part):
                rsid = m.group(0)
                if rsid not in tiers:
                    tiers[rsid] = current_tier

    return tiers


def synthesize_interpretation(scores: list[ScoreRecord]) -> str:
    """Generate a fallback interpretation from score data when markdown has none."""
    if not scores:
        return "No regulatory effect scores available for this variant."

    top = scores[0]
    genes = list(dict.fromkeys(s.gene_name for s in scores[:3] if s.gene_name))
    assays = list(dict.fromkeys(s.assay for s in scores[:3]))
    biosamples = list(
        dict.fromkeys(s.biosample_name for s in scores[:3] if s.biosample_name)
    )

    gene_str = ", ".join(genes[:3]) if genes else "nearby genes"
    assay_str = ", ".join(assays[:3]) if assays else "multiple assays"
    bio_str = ", ".join(biosamples[:2]) if biosamples else "relevant biosamples"

    return (
        f"Strongest regulatory signal in {top.assay} for {top.gene_name} "
        f"(|effect|={top.abs_effect:.2f}) in {top.biosample_name}. "
        f"Additional evidence from {assay_str} across {bio_str}. "
        f"Consistent with regulatory activity at the {gene_str} locus."
    )


# ---------------------------------------------------------------------------
# Output filename helper
# ---------------------------------------------------------------------------


def _make_output_filename(
    chromosome: str,
    position: int,
    rsid: str = "",
    prefix: str = "variant_report",
) -> str:
    """Create a sanitized PDF filename for a variant report."""
    clean_parts = [re.sub(r"[^\w-]", "", prefix), re.sub(r"[^\w-]", "", chromosome), str(position)]
    if rsid:
        clean_parts.append(re.sub(r"[^\w-]", "", rsid))
    return "_".join(clean_parts) + ".pdf"


# ---------------------------------------------------------------------------
# Variant-from-ID helper
# ---------------------------------------------------------------------------


def _variant_from_id(vcf_id: str) -> CandidateVariant:
    """Derive a CandidateVariant from a VCF ID string like 'chr8:127401060:G>T'."""
    m = _VARIANT_ID_PATTERN.match(vcf_id)
    if not m:
        raise ValueError(f"Cannot parse variant ID: {vcf_id}")
    return CandidateVariant(
        chromosome=f"chr{m.group(1)}",
        position=int(m.group(2)),
        ref=m.group(3),
        alt=m.group(4),
    )


# ---------------------------------------------------------------------------
# ReportBuilder
# ---------------------------------------------------------------------------


class ReportBuilder:
    """Build and render one-page PDF variant reports."""

    def __init__(self, template_name: str = "default", tissue: str | None = None) -> None:
        self._template_name = template_name
        self._tissue = tissue
        self._env = None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_env(self):
        if self._env is None:
            import jinja2  # lazy import

            template_dir = Path(__file__).resolve().parent / "templates"
            self._env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(str(template_dir)),
                autoescape=False,
            )
        return self._env

    @staticmethod
    def _parse_rsid_table(md_text: str) -> dict[str, str]:
        """Extract {rsid → vcf_id} mappings from the ranked variant table."""
        mapping: dict[str, str] = {}
        for line in md_text.splitlines():
            m = re.search(
                r"\|\s*\**\s*(rs\d+)\s*\**\s*\|\s*\**\s*(chr\d+:\d+:[ACGT]+>[ACGT]+)\s*\**\s*\|",
                line,
            )
            if m:
                mapping[m.group(1)] = m.group(2)
        return mapping

    # ------------------------------------------------------------------
    # build
    # ------------------------------------------------------------------

    def build(
        self,
        candidates_tsv: str | Path,
        scores_tsv: str | Path,
        report_md: str | Path,
    ) -> list[CandidateReport]:
        """Parse all inputs and return a list of CandidateReports (one per variant).

        Parameters
        ----------
        candidates_tsv:
            TSV of known candidate variants (soft-fail if missing).
        scores_tsv:
            TSV of per-assay scores.
        report_md:
            Agent-generated markdown report.
        tissue:
            Optional tissue context for the footer.
        """
        # --- parse scores ---
        scores_by_vcf = parse_scores(scores_tsv)

        # --- parse markdown report ---
        md_path = Path(report_md)
        interpretations = extract_interpretations(md_path)
        tiers = assign_tiers(md_path)

        # --- extract rsID → vcf_id from the markdown rank table ---
        md_text = md_path.read_text()
        rsid_to_vcf = self._parse_rsid_table(md_text)
        vcf_to_rsid = {v: r for r, v in rsid_to_vcf.items()}

        # --- read candidates (soft-fail) ---
        candidates: list[CandidateVariant] = []
        try:
            candidates = read_candidates_tsv(candidates_tsv)
        except (FileNotFoundError, ValueError):
            pass

        candidate_by_vcf = {v.vcf_id: v for v in candidates}

        # --- build one CandidateReport per scored variant ---
        reports: list[CandidateReport] = []
        for vcf_id, scores in scores_by_vcf.items():
            variant = candidate_by_vcf.get(vcf_id)
            if variant is None:
                variant = _variant_from_id(vcf_id)

            rsid = vcf_to_rsid.get(vcf_id, "")

            interpretation = interpretations.get(vcf_id, "")
            if not interpretation:
                interpretation = synthesize_interpretation(scores)

            tier = (
                tiers.get(vcf_id)
                or tiers.get(rsid)
                or "Unranked"
            )

            max_effect = scores[0].abs_effect if scores else 0.0

            reports.append(
                CandidateReport(
                    variant=variant,
                    vcf_id=vcf_id,
                    rsid=rsid,
                    cytoband=variant.region_id or "",
                    top_scores=scores,
                    max_effect=max_effect,
                    interpretation=interpretation,
                    tier=tier,
                    validation_plan="",
                )
            )

        return reports

    # ------------------------------------------------------------------
    # render
    # ------------------------------------------------------------------

    def render(
        self,
        report: CandidateReport,
        tissue: str | None = None,
        locus_figure_path: str | None = None,
        barplot_figure_path: str | None = None,
    ) -> str:
        """Render a CandidateReport to a LaTeX string via Jinja2 template."""
        env = self._get_env()
        template = env.get_template("report.tex.j2")

        assays = {s.assay for s in report.top_scores}
        top_genes = list(
            dict.fromkeys(s.gene_name for s in report.top_scores[:3] if s.gene_name)
        )
        gene_context = ", ".join(top_genes) if top_genes else ""

        tier_colors = {
            "Tier 1": "059669",
            "Tier 2": "d97706",
            "Unranked": "6b7280",
        }
        tier_color = tier_colors.get(report.tier, "6b7280")

        tex = template.render(
            variant_id=report.vcf_id,
            rsid=report.rsid,
            cytoband=report.cytoband,
            gene_context=gene_context,
            tier=report.tier,
            tier_color=tier_color,
            top_scores=report.top_scores,
            interpretation=report.interpretation,
            validation_plan=report.validation_plan,
            assays_count=len(assays),
            tissue=tissue or "",
            date=date.today().strftime("%Y-%m-%d"),
            locus_figure=locus_figure_path or "",
            barplot_figure=barplot_figure_path or "",
        )
        return tex

    # ------------------------------------------------------------------
    # generate_pdf
    # ------------------------------------------------------------------

    def generate_pdf(
        self,
        report: CandidateReport,
        output_dir: str | Path,
        prefix: str = "variant_report",
        tissue: str | None = None,
        compile_pdf: bool = True,
        clean: bool = True,
        verbose: bool = False,
    ) -> Path:
        """Render LaTeX, generate figures, and optionally compile to PDF.

        Returns the path of the generated PDF (or .tex if ``compile=False``).
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        from regvar.plot import plot_effects, plot_sequence_context

        with tempfile.TemporaryDirectory(prefix="regvar_") as tmp_dir:
            tmp = Path(tmp_dir)

            # --- figures ---
            locus_fig = plot_sequence_context(report.variant, report.top_scores)
            locus_path = tmp / "locus_context.png"
            locus_fig.savefig(str(locus_path), dpi=150, bbox_inches="tight")

            score_dicts = [
                {
                    "assay": s.assay,
                    "gene_name": s.gene_name,
                    "biosample_name": s.biosample_name,
                    "output_type": s.output_type,
                    "raw_score": s.raw_score,
                    "quantile_score": s.quantile_score,
                    "abs_effect": s.abs_effect,
                }
                for s in report.top_scores
            ]
            bar_fig = plot_effects(score_dicts, title=report.vcf_id)
            bar_path = tmp / "barplot_effects.png"
            bar_fig.savefig(str(bar_path), dpi=150, bbox_inches="tight")

            # --- LaTeX ---
            tex_content = self.render(
                report,
                tissue=tissue,
                locus_figure_path=str(locus_path),
                barplot_figure_path=str(bar_path),
            )
            tex_path = tmp / "report.tex"
            tex_path.write_text(tex_content)

            if compile_pdf:
                for _ in range(2):
                    result = subprocess.run(
                        [
                            "pdflatex",
                            "-interaction=nonstopmode",
                            "-output-directory",
                            str(tmp),
                            str(tex_path),
                        ],
                        capture_output=True,
                        timeout=30,
                        text=True,
                    )
                    if verbose:
                        if result.stdout:
                            print(result.stdout)
                        if result.stderr:
                            print(result.stderr)

                # Issue 3: Check pdflatex return code
                if result.returncode != 0:
                    raise RuntimeError(
                        f"LaTeX compilation failed. Check log: {tex_path.with_suffix('.log')}"
                    )

                pdf_path = tmp / "report.pdf"
                if pdf_path.exists():
                    out_name = _make_output_filename(
                        report.variant.chromosome,
                        report.variant.position,
                        report.rsid,
                        prefix,
                    )
                    out_path = output_dir / out_name
                    import shutil

                    shutil.copy2(str(pdf_path), str(out_path))

                    # Issue 2: Clean up auxiliary files from output_dir
                    if clean:
                        for suffix in (".aux", ".log", ".out"):
                            aux_file = out_path.with_suffix(suffix)
                            if aux_file.exists():
                                aux_file.unlink()

                    return out_path

                raise RuntimeError("pdflatex completed but PDF was not produced.")

            # Issue 1: compile=False — copy .tex to output_dir before returning
            out_name = _make_output_filename(
                report.variant.chromosome,
                report.variant.position,
                report.rsid,
                prefix,
            )
            tex_out_path = output_dir / out_name.replace(".pdf", ".tex")
            import shutil
            shutil.copy2(str(tex_path), str(tex_out_path))
            return tex_out_path

    # ------------------------------------------------------------------
    # generate_all
    # ------------------------------------------------------------------

    def generate_all(
        self,
        reports: list[CandidateReport],
        output_dir: str | Path,
        top_n: int | None = None,
        prefix: str = "variant_report",
        tissue: str | None = None,
        compile_pdf: bool = True,
        clean: bool = True,
        verbose: bool = False,
    ) -> list[Path]:
        """Generate PDFs for *top_n* (or all) reports."""
        if top_n is not None:
            reports = reports[:top_n]

        paths: list[Path] = []
        for report in reports:
            p = self.generate_pdf(
                report,
                output_dir=output_dir,
                prefix=prefix,
                tissue=tissue,
                compile_pdf=compile_pdf,
                clean=clean,
                verbose=verbose,
            )
            paths.append(p)
        return paths
