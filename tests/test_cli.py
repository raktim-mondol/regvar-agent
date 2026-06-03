"""Tests for the CLI (cli.py).

Uses Click's CliRunner for command testing. All tests run without API keys.
"""

import pytest
from click.testing import CliRunner

from regvar.cli import _parse_assays, cli
from regvar.variants import write_scores_tsv


@pytest.fixture
def runner():
    return CliRunner()


# --- _parse_assays --------------------------------------------------------

def test_parse_assays_none():
    assert _parse_assays(None) is None


def test_parse_assays_all():
    result = _parse_assays("all")
    assert isinstance(result, list)
    assert "ATAC-seq" in result
    assert "RNA-seq" in result
    assert "splicing" in result
    assert len(result) == 7


def test_parse_assays_all_case_insensitive():
    assert _parse_assays("ALL") == _parse_assays("all")
    assert _parse_assays("All") == _parse_assays("all")


def test_parse_assays_csv():
    result = _parse_assays("ATAC-seq,RNA-seq")
    assert result == ["ATAC-seq", "RNA-seq"]


def test_parse_assays_csv_with_spaces():
    result = _parse_assays(" ATAC-seq , RNA-seq ")
    assert result == ["ATAC-seq", "RNA-seq"]


def test_parse_assays_empty_string():
    result = _parse_assays("")
    assert result == []


# --- version & help -------------------------------------------------------

def test_version(runner):
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert "0.1.0" in result.output


def test_help(runner):
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "regvar" in result.output
    assert "score" in result.output
    assert "agent" in result.output


def test_score_help(runner):
    result = runner.invoke(cli, ["score", "--help"])
    assert result.exit_code == 0
    assert "--output-tsv" in result.output
    assert "--output-png" in result.output
    assert "--assays" in result.output
    assert "chromosome" in result.output.lower()


def test_agent_run_help(runner):
    result = runner.invoke(cli, ["agent", "run", "--help"])
    assert result.exit_code == 0
    assert "--dry-run" in result.output
    assert "--output-json" in result.output
    assert "--force" in result.output


def test_assays_command(runner):
    result = runner.invoke(cli, ["assays"])
    assert result.exit_code == 0
    assert "ATAC-seq" in result.output
    assert "RNA-seq" in result.output


# --- dry-run --------------------------------------------------------------

def test_dry_run_valid(runner, tmp_path):
    tsv = tmp_path / "candidates.tsv"
    tsv.write_text(
        "chromosome\tposition\tref\talt\tregion_id\tnote\n"
        "chr8\t127401060\tG\tT\tenh1\ttest\n"
    )
    result = runner.invoke(cli, ["agent", "run", str(tsv), "--dry-run"])
    assert result.exit_code == 0
    assert "Dry run" in result.output or "dry run" in result.output.lower()
    assert "Ready" in result.output or "validated" in result.output.lower()


def test_dry_run_multiple_variants(runner, tmp_path):
    tsv = tmp_path / "multi.tsv"
    tsv.write_text(
        "chromosome\tposition\tref\talt\n"
        "chr8\t127401060\tG\tT\n"
        "chr10\t46046326\tA\tG\n"
    )
    result = runner.invoke(cli, ["agent", "run", str(tsv), "--dry-run"])
    assert result.exit_code == 0
    assert "2" in result.output


def test_dry_run_invalid_path(runner):
    result = runner.invoke(cli, ["agent", "run", "/nonexistent/file.tsv", "--dry-run"])
    assert result.exit_code != 0


def test_dry_run_bad_tsv(runner, tmp_path):
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("chromosome\tposition\tref\nchr8\t1\tG\n")
    result = runner.invoke(cli, ["agent", "run", str(tsv), "--dry-run"])
    assert result.exit_code != 0


def test_dry_run_with_assays_and_tissue(runner, tmp_path):
    tsv = tmp_path / "c.tsv"
    tsv.write_text(
        "chromosome\tposition\tref\talt\n"
        "chr8\t127401060\tG\tT\n"
    )
    result = runner.invoke(cli, [
        "agent", "run", str(tsv), "--dry-run",
        "--assays", "ATAC-seq,RNA-seq",
        "--tissue", "UBERON:0002367",
        "--top-n", "5",
    ])
    assert result.exit_code == 0
    assert "ATAC-seq" in result.output
    assert "UBERON:0002367" in result.output


# --- write_scores_tsv (used by CLI) --------------------------------------

def test_write_scores_tsv(tmp_path):
    records = [
        {"gene_name": "MYC", "raw_score": 0.5, "assay": "ATAC-seq"},
        {"gene_name": "TP53", "raw_score": -0.3, "assay": "RNA-seq"},
    ]
    out = tmp_path / "scores.tsv"
    write_scores_tsv(records, out)

    lines = out.read_text().strip().split("\n")
    assert lines[0] == "gene_name\traw_score\tassay"
    assert lines[1] == "MYC\t0.5\tATAC-seq"
    assert len(lines) == 3


def test_write_scores_tsv_empty(tmp_path):
    out = tmp_path / "empty.tsv"
    write_scores_tsv([], out)
    assert not out.exists()


# --- cache commands -------------------------------------------------------

def test_cache_info(runner):
    result = runner.invoke(cli, ["cache", "info"])
    assert result.exit_code == 0
    assert "Cache directory" in result.output
    assert "Entries" in result.output


def test_cache_clear_no_entries(runner):
    result = runner.invoke(cli, ["cache", "clear", "--yes"])
    assert result.exit_code == 0


# --- backward-compat run alias -------------------------------------------

def test_run_alias_help(runner):
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0


# --- plot subcommand help -------------------------------------------------

def test_plot_effects_help(runner):
    result = runner.invoke(cli, ["plot", "effects", "--help"])
    assert result.exit_code == 0
    assert "--output" in result.output
    assert "--show" in result.output
    assert "--from-tsv" in result.output
