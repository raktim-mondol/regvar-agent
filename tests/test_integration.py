"""Integration tests — gated behind RUN_INTEGRATION_TESTS=1 and ALPHAGENOME_API_KEY.

These tests hit the real AlphaGenome and DeepSeek APIs. They are skipped by
default and only run when explicitly enabled for manual smoke-testing.

To run:
    RUN_INTEGRATION_TESTS=1 ALPHAGENOME_API_KEY=... DEEPSEEK_API_KEY=... pytest tests/test_integration.py -v
"""

import os

import pytest

from regvar.tools import tool_score_regulatory_variant

pytestmark = pytest.mark.skipif(
    not os.environ.get("RUN_INTEGRATION_TESTS"),
    reason="Set RUN_INTEGRATION_TESTS=1 and ALPHAGENOME_API_KEY to run",
)


class TestIntegrationScoreVariant:
    def test_score_known_variant(self):
        """rs6983267 at chr8:127401060 G>T — a well-studied prostate-cancer SNP."""
        result = tool_score_regulatory_variant("chr8", 127401060, "G", "T")

        assert result["variant"] == "chr8:127401060:G>T"
        assert result["n_total_scores"] > 0
        assert len(result["top_effects"]) > 0

        first = result["top_effects"][0]
        assert "gene_name" in first
        assert "raw_score" in first

    def test_score_default_assays(self):
        result = tool_score_regulatory_variant("chr8", 127401060, "G", "T")
        assert "ATAC-seq" in result["assays_scored"]
        assert "RNA-seq" in result["assays_scored"]

    def test_score_custom_assay(self):
        result = tool_score_regulatory_variant(
            "chr8", 127401060, "G", "T",
            assays=["ATAC-seq"],
        )
        assert result["assays_scored"] == ["ATAC-seq"]
        assert result["n_total_scores"] > 0


class TestIntegrationAgent:
    def test_agent_run_dry_run(self, tmp_path):
        """Verify the agent can at least load candidates and produce a task message."""
        from regvar.agent import build_task_message
        from regvar.variants import read_candidates_tsv

        tsv = tmp_path / "test.tsv"
        tsv.write_text(
            "chromosome\tposition\tref\talt\tregion_id\n"
            "chr8\t127401060\tG\tT\trs6983267\n"
        )
        candidates = read_candidates_tsv(tsv)
        msg = build_task_message(candidates)

        assert "chr8:127401060:G>T" in msg
        assert "rs6983267" in msg
