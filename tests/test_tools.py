"""Tests for the tool boundary layer (tools.py).

All tests mock the AlphaGenome client so they run without API keys or network.
"""

import json
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from regvar.tools import (
    OPENAI_TOOL_SCHEMAS,
    TOOL_DISPATCH,
    TOOL_SCHEMAS,
    run_tool,
    tool_list_supported_assays,
    tool_score_regulatory_variant,
    tool_score_variants_batch,
)
from regvar.variants import MAX_VARIANTS


def _make_tidy_df(n=3):
    """Minimal tidy DataFrame mimicking AlphaGenome score_variant output."""
    return pd.DataFrame({
        "gene_name": [f"GENE{i}" for i in range(n)],
        "raw_score": [0.1 * (i + 1) for i in range(n)],
        "quantile_score": [0.2 * (i + 1) for i in range(n)],
        "assay": ["ATAC-seq"] * n,
        "output_type": ["ATAC"] * n,
        "biosample_name": ["prostate_fibroblast"] * n,
    })


@patch("regvar.tools.get_client")
def test_tool_list_assays(mock_get_client):
    mock_client = MagicMock()
    mock_client.supported_assays.return_value = {
        "ATAC-seq": "ATAC",
        "RNA-seq": "RNA_SEQ",
    }
    mock_get_client.return_value = mock_client

    result = tool_list_supported_assays()

    assert "assays" in result
    assert result["assays"]["ATAC-seq"] == "ATAC"
    assert result["assays"]["RNA-seq"] == "RNA_SEQ"


@patch("regvar.tools.get_client")
def test_tool_score_variant_structure(mock_get_client):
    mock_client = MagicMock()
    mock_client.score_variant.return_value = _make_tidy_df()
    mock_get_client.return_value = mock_client

    result = tool_score_regulatory_variant("chr8", 127401060, "G", "T")

    assert result["variant"] == "chr8:127401060:G>T"
    assert isinstance(result["assays_scored"], list)
    assert result["n_total_scores"] == 3
    assert len(result["top_effects"]) == 3
    assert result["top_effects"][0]["gene_name"] == "GENE2"  # ranked by abs effect


@patch("regvar.tools.get_client")
def test_tool_score_variant_custom_assays(mock_get_client):
    mock_client = MagicMock()
    mock_client.score_variant.return_value = _make_tidy_df(1)
    mock_get_client.return_value = mock_client

    result = tool_score_regulatory_variant(
        "chr8", 100, "A", "C",
        assays=["RNA-seq"],
        ontology_terms=["UBERON:9999999"],
    )

    mock_client.score_variant.assert_called_once_with(
        chromosome="chr8", position=100, ref="A", alt="C",
        assays=["RNA-seq"],
        ontology_terms=["UBERON:9999999"],
    )
    assert result["assays_scored"] == ["RNA-seq"]


@patch("regvar.tools.get_client")
def test_tool_score_variants_batch_structure(mock_get_client):
    tidy = _make_tidy_df(2)
    tidy["query_variant"] = ["chr8:100:A>G"] * 2
    mock_client = MagicMock()
    mock_client.score_variants.return_value = tidy
    mock_get_client.return_value = mock_client

    variants = [{"chromosome": "chr8", "position": 100, "ref": "A", "alt": "G"}]
    result = tool_score_variants_batch(variants)

    assert result["variants_scored"] == 1
    assert result["n_total_scores"] == 2
    assert "chr8:100:A>G" in result["effects_by_variant"]


@patch("regvar.tools.get_client")
def test_tool_score_variants_batch_exceeds_limit(mock_get_client):
    too_many = [
        {"chromosome": "chr1", "position": i, "ref": "A", "alt": "G"}
        for i in range(MAX_VARIANTS + 1)
    ]
    with pytest.raises(ValueError, match="exceeding the limit"):
        tool_score_variants_batch(too_many)

    mock_get_client.assert_not_called()


@patch("regvar.tools.get_client")
def test_run_tool_known(mock_get_client):
    mock_client = MagicMock()
    mock_client.supported_assays.return_value = {"ATAC-seq": "ATAC"}
    mock_get_client.return_value = mock_client

    result_json = run_tool("list_supported_assays", {})
    parsed = json.loads(result_json)

    assert parsed["assays"]["ATAC-seq"] == "ATAC"


def test_run_tool_unknown():
    result_json = run_tool("nonexistent_tool", {})
    parsed = json.loads(result_json)

    assert "error" in parsed
    assert "nonexistent_tool" in parsed["error"]


@patch("regvar.tools.get_client")
def test_run_tool_exception(mock_get_client):
    mock_client = MagicMock()
    mock_client.score_variant.side_effect = RuntimeError("API down")
    mock_get_client.return_value = mock_client

    result_json = run_tool("score_regulatory_variant", {
        "chromosome": "chr8", "position": 100, "ref": "A", "alt": "G",
    })
    parsed = json.loads(result_json)

    assert "error" in parsed
    assert "API down" in parsed["error"]


def test_dispatch_matches_schemas():
    """TOOL_DISPATCH keys must exactly match the names in both schema formats."""
    dispatch_names = set(TOOL_DISPATCH.keys())
    anthropic_names = {s["name"] for s in TOOL_SCHEMAS}
    openai_names = {s["function"]["name"] for s in OPENAI_TOOL_SCHEMAS}

    assert dispatch_names == anthropic_names, (
        f"TOOL_DISPATCH vs TOOL_SCHEMAS mismatch: "
        f"extra={dispatch_names - anthropic_names}, missing={anthropic_names - dispatch_names}"
    )
    assert dispatch_names == openai_names, (
        f"TOOL_DISPATCH vs OPENAI_TOOL_SCHEMAS mismatch: "
        f"extra={dispatch_names - openai_names}, missing={openai_names - dispatch_names}"
    )


@patch("regvar.tools.get_client")
def test_tool_score_variant_chromosome_validation(mock_get_client):
    mock_client = MagicMock()
    mock_client.score_variant.return_value = _make_tidy_df(1)
    mock_get_client.return_value = mock_client

    result = tool_score_regulatory_variant("8", 100, "A", "G")
    assert result["variant"] == "chr8:100:A>G"

    with pytest.raises(ValueError, match="Invalid chromosome"):
        tool_score_regulatory_variant("chr99", 100, "A", "G")
