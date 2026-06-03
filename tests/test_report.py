"""Tests for report generator — data model, parsing, rendering."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Ensure regvar is importable in test runs
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from regvar.report import (
    ScoreRecord, CandidateReport, ReportBuilder, parse_scores,
    extract_interpretations, assign_tiers, synthesize_interpretation,
)
from regvar.variants import CandidateVariant

FIXTURES = Path(__file__).resolve().parent / "fixtures"


class TestScoreRecord:
    def test_from_dict(self):
        d = {
            "assay": "ATAC-seq", "gene_name": "MYC",
            "biosample_name": "prostate", "output_type": "ATAC",
            "raw_score": "0.52", "quantile_score": "0.91", "abs_effect": "0.91",
        }
        sr = ScoreRecord(
            assay=d["assay"], gene_name=d["gene_name"],
            biosample_name=d["biosample_name"], output_type=d["output_type"],
            raw_score=float(d["raw_score"]), quantile_score=float(d["quantile_score"]),
            abs_effect=float(d["abs_effect"]),
        )
        assert sr.assay == "ATAC-seq"
        assert sr.gene_name == "MYC"
        assert sr.quantile_score == 0.91

    def test_fields_match_spec(self):
        sr = ScoreRecord(
            assay="ATAC-seq", gene_name="MYC", biosample_name="prostate",
            output_type="ATAC", raw_score=0.5, quantile_score=0.9, abs_effect=0.9,
        )
        assert hasattr(sr, "assay")
        assert hasattr(sr, "gene_name")
        assert hasattr(sr, "biosample_name")
        assert hasattr(sr, "output_type")
        assert hasattr(sr, "raw_score")
        assert hasattr(sr, "quantile_score")
        assert hasattr(sr, "abs_effect")


class TestParseScores:
    def test_groups_by_variant(self):
        result = parse_scores(FIXTURES / "example_scores.tsv")
        assert len(result) == 2
        assert "chr8:127401060:G>T" in result
        assert "chr10:46046326:A>G" in result

    def test_sorted_by_abs_effect_desc(self):
        result = parse_scores(FIXTURES / "example_scores.tsv")
        scores = result["chr8:127401060:G>T"]
        assert scores[0].abs_effect >= scores[-1].abs_effect
        assert scores[0].abs_effect == 0.91

    def test_empty_tsv_raises(self):
        import tempfile
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("query_variant\tgene_name\tassay\toutput_type\tbiosample_name\traw_score\tquantile_score\tabs_effect\n")
            tmp = Path(f.name)
        try:
            with pytest.raises(ValueError, match="No scored variants"):
                parse_scores(tmp)
        finally:
            tmp.unlink()

    def test_missing_columns_raises(self):
        import tempfile
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("variant_id\tgene_name\nchr8:1:G>T\tMYC\n")
            tmp = Path(f.name)
        try:
            with pytest.raises(ValueError, match="Missing columns"):
                parse_scores(tmp)
        finally:
            tmp.unlink()


class TestCandidateReport:
    def test_construction(self):
        variant = CandidateVariant("chr8", 127401060, "G", "T", "8q24_enh", "MYC locus")
        scores = [
            ScoreRecord("ATAC-seq", "MYC", "prostate", "ATAC", 0.52, 0.91, 0.91),
            ScoreRecord("RNA-seq", "MYC", "prostate", "RNA_SEQ", 0.41, 0.79, 0.79),
        ]
        cr = CandidateReport(
            variant=variant, vcf_id=variant.vcf_id, rsid="rs6983267",
            cytoband="8q24.21", top_scores=scores, max_effect=0.91,
            interpretation="Strong ATAC-seq signal.", tier="Tier 1",
            validation_plan="Luciferase assay in LNCaP.",
        )
        assert cr.vcf_id == "chr8:127401060:G>T"
        assert cr.max_effect == 0.91
        assert cr.tier == "Tier 1"
        assert len(cr.top_scores) == 2


class TestExtractInterpretations:
    def test_extracts_by_variant_id(self):
        result = extract_interpretations(FIXTURES / "example_report.md")
        assert "chr8:127401060:G>T" in result
        assert "chr10:46046326:A>G" in result

    def test_extracted_text_not_empty(self):
        result = extract_interpretations(FIXTURES / "example_report.md")
        text = result["chr8:127401060:G>T"]
        assert len(text) > 20
        assert "MYC" in text or "enhancer" in text or "G>T" in text

    def test_handles_minimal_report(self):
        result = extract_interpretations(FIXTURES / "example_report_minimal.md")
        assert "chr8:127401060:G>T" in result
        assert "chr10:46046326:A>G" in result


class TestAssignTiers:
    def test_extracts_tiers(self):
        result = assign_tiers(FIXTURES / "example_report.md")
        assert result.get("rs6983267") == "Tier 1" or result.get("chr8:127401060:G>T") == "Tier 1"
        assert result.get("rs10993994") == "Tier 1" or result.get("chr10:46046326:A>G") == "Tier 1"

    def test_minimal_report_returns_empty(self):
        result = assign_tiers(FIXTURES / "example_report_minimal.md")
        assert result == {}


class TestSynthesizeInterpretation:
    def test_synthesizes_from_scores(self):
        scores = [
            ScoreRecord("ATAC-seq", "MYC", "prostate", "ATAC", 0.52, 0.91, 0.91),
            ScoreRecord("RNA-seq", "MYC", "prostate", "RNA_SEQ", 0.41, 0.79, 0.79),
        ]
        text = synthesize_interpretation(scores)
        assert "ATAC-seq" in text
        assert "MYC" in text
        assert "0.91" in text
        assert len(text) > 30


class TestPlotSequenceContext:
    def test_returns_figure(self):
        from regvar.plot import plot_sequence_context
        from regvar.variants import CandidateVariant

        variant = CandidateVariant("chr8", 127401060, "G", "T", "8q24_enh", "MYC locus")
        scores = [
            ScoreRecord("ATAC-seq", "MYC", "prostate", "ATAC", 0.52, 0.91, 0.91),
            ScoreRecord("RNA-seq", "MYC", "prostate", "RNA_SEQ", 0.41, 0.79, 0.79),
        ]
        fig = plot_sequence_context(variant, scores)
        assert fig is not None
        # Should be a matplotlib Figure
        assert hasattr(fig, "savefig")

    def test_figure_dimensions(self):
        from regvar.plot import plot_sequence_context
        from regvar.variants import CandidateVariant

        variant = CandidateVariant("chr8", 127401060, "G", "T")
        scores = [ScoreRecord("ATAC-seq", "MYC", "prostate", "ATAC", 0.52, 0.91, 0.91)]
        fig = plot_sequence_context(variant, scores, figsize=(4.5, 2.0))
        size = fig.get_size_inches()
        assert size[0] == 4.5
        assert size[1] == 2.0

    def test_handles_empty_scores(self):
        from regvar.plot import plot_sequence_context
        from regvar.variants import CandidateVariant

        variant = CandidateVariant("chr8", 127401060, "G", "T")
        fig = plot_sequence_context(variant, [])
        assert fig is not None


class TestReportBuilderBuild:
    def test_build_returns_reports(self):
        builder = ReportBuilder(tissue="Prostate")
        candidates_tsv = Path(__file__).resolve().parent.parent / "examples" / "candidate_variants.tsv"
        reports = builder.build(
            candidates_tsv=candidates_tsv,
            scores_tsv=FIXTURES / "example_scores.tsv",
            report_md=FIXTURES / "example_report.md",
        )
        assert len(reports) > 0
        for r in reports:
            assert r.vcf_id
            assert r.top_scores
            assert r.max_effect > 0
            assert r.interpretation
            assert r.tier in ("Tier 1", "Tier 2", "Unranked")

    def test_build_assigns_rsids(self):
        builder = ReportBuilder()
        candidates_tsv = Path(__file__).resolve().parent.parent / "examples" / "candidate_variants.tsv"
        reports = builder.build(
            candidates_tsv=candidates_tsv,
            scores_tsv=FIXTURES / "example_scores.tsv",
            report_md=FIXTURES / "example_report.md",
        )
        rs_variants = [r for r in reports if r.rsid]
        assert len(rs_variants) > 0

    def test_build_handles_missing_candidates_tsv(self):
        builder = ReportBuilder()
        reports = builder.build(
            candidates_tsv=Path("/nonexistent/path.tsv"),
            scores_tsv=FIXTURES / "example_scores.tsv",
            report_md=FIXTURES / "example_report.md",
        )
        assert len(reports) > 0  # Should still work from scores alone


class TestReportBuilderRender:
    def test_render_produces_latex(self):
        builder = ReportBuilder()
        from regvar.variants import CandidateVariant
        variant = CandidateVariant("chr8", 127401060, "G", "T", "8q24", "MYC")
        scores = [
            ScoreRecord("ATAC-seq", "MYC", "prostate", "ATAC", 0.52, 0.91, 0.91),
            ScoreRecord("RNA-seq", "MYC", "prostate", "RNA_SEQ", 0.41, 0.79, 0.79),
        ]
        cr = CandidateReport(
            variant=variant, vcf_id=variant.vcf_id, rsid="rs6983267",
            cytoband="8q24.21", top_scores=scores, max_effect=0.91,
            interpretation="Strong ATAC-seq and RNA-seq signals at the MYC locus.",
            tier="Tier 1",
            validation_plan="Luciferase reporter assay in LNCaP cells.",
        )
        tex = builder.render(cr, tissue="Prostate")
        assert r"\documentclass" in tex
        assert "chr8:127401060:G>T" in tex
        assert "rs6983267" in tex
        assert "Tier 1" in tex
        assert "MYC" in tex
        assert "Luciferase" in tex
        # No unresolved Jinja2 variables
        assert "{{" not in tex
        assert "}}" not in tex


class TestOutputFilename:
    def test_with_rsid(self):
        from regvar.report import _make_output_filename
        name = _make_output_filename("chr8", 127401060, "rs6983267", "variant_report")
        assert "chr8" in name
        assert "127401060" in name
        assert "rs6983267" in name
        assert name.endswith(".pdf")

    def test_without_rsid(self):
        from regvar.report import _make_output_filename
        name = _make_output_filename("chr8", 127401060, "", "variant_report")
        assert "rs" not in name.lower()
        assert name.endswith(".pdf")

    def test_with_prefix(self):
        from regvar.report import _make_output_filename
        name = _make_output_filename("chr8", 127401060, "rs6983267", "patient_report")
        assert name.startswith("patient_report")


class TestIntegrationNoLatex:
    """End-to-end test that doesn't require pdflatex."""

    def test_full_pipeline_no_compile(self, tmp_path):
        builder = ReportBuilder()
        candidates_tsv = Path(__file__).resolve().parent.parent / "examples" / "candidate_variants.tsv"

        reports = builder.build(
            candidates_tsv=candidates_tsv,
            scores_tsv=FIXTURES / "example_scores.tsv",
            report_md=FIXTURES / "example_report.md",
        )
        assert len(reports) >= 2

        # Generate .tex only (no pdflatex)
        paths = builder.generate_all(
            reports=reports,
            output_dir=tmp_path,
            top_n=3,
            prefix="test_report",
            tissue="Prostate",
            compile_pdf=False,
        )
        assert len(paths) == min(3, len(reports))
        for p in paths:
            assert p.suffix == ".tex"
            assert p.exists()
            content = p.read_text()
            assert r"\documentclass" in content
            assert "{{" not in content  # No unresolved Jinja2
            assert "}}" not in content
