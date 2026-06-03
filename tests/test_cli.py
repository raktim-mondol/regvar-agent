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


# --- chat helpers --------------------------------------------------------

def test_serialize_messages_plain_dict():
    from regvar.cli import _serialize_messages
    msgs = [
        {"role": "user", "content": "hi"},
        {"role": "assistant", "content": "hello"},
    ]
    out = _serialize_messages(msgs)
    assert out == msgs


def test_serialize_messages_pydantic_object():
    """Pydantic-style message objects get dumped via model_dump."""
    from regvar.cli import _serialize_messages

    class _Msg:
        def model_dump(self, exclude_unset=True):
            return {"role": "assistant", "content": "dumped"}

    out = _serialize_messages([_Msg()])
    assert out == [{"role": "assistant", "content": "dumped"}]


def test_save_load_session_roundtrip(tmp_path):
    from regvar.cli import _load_session, _save_session
    msgs = [
        {"role": "system", "content": "sys"},
        {"role": "user", "content": "q"},
        {"role": "assistant", "content": "a"},
    ]
    path = tmp_path / "session.json"
    _save_session(msgs, path)
    loaded = _load_session(path)
    assert loaded == msgs


def test_load_session_invalid_shape(tmp_path):
    from regvar.cli import _load_session
    path = tmp_path / "bad.json"
    path.write_text('{"not_messages": []}', encoding="utf-8")
    import pytest
    with pytest.raises(ValueError, match="messages"):
        _load_session(path)


def test_render_conversation_markdown():
    from regvar.cli import _render_conversation_markdown
    msgs = [
        {"role": "system", "content": "system prompt"},
        {"role": "user", "content": "score chr8:100 A>G"},
        {"role": "assistant", "content": "here is the answer"},
    ]
    md = _render_conversation_markdown(msgs)
    assert "# Agent Chat Session" in md
    assert "## User" in md
    assert "score chr8:100 A>G" in md
    assert "here is the answer" in md


def test_render_conversation_markdown_with_tool_calls():
    from regvar.cli import _render_conversation_markdown
    msgs = [
        {
            "role": "assistant",
            "content": None,
            "tool_calls": [{
                "function": {
                    "name": "score_regulatory_variant",
                    "arguments": '{"chromosome":"chr8"}',
                },
            }],
        },
    ]
    md = _render_conversation_markdown(msgs)
    assert "score_regulatory_variant" in md


def test_agent_chat_help_lists_new_flags(runner):
    result = runner.invoke(cli, ["agent", "chat", "--help"])
    assert result.exit_code == 0
    assert "--save-session" in result.output
    assert "--load-session" in result.output
    assert "--reasoning-effort" in result.output


def test_agent_run_help_lists_annotate_and_max_concurrent(runner):
    result = runner.invoke(cli, ["agent", "run", "--help"])
    assert result.exit_code == 0
    assert "--annotate" in result.output
    assert "--max-concurrent" in result.output


# --- config default_map integration --------------------------------------

def test_config_default_map_score(runner, tmp_path, monkeypatch):
    """TOML [defaults] section propagates to subcommands at runtime.

    Click's help formatter shows the *static* default (10 for --top-n), not
    the runtime default_map override — so we exercise the value via
    ``agent run --dry-run`` which echoes the effective top_n.
    """
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("[defaults]\ntop_n = 42\n", encoding="utf-8")

    tsv = tmp_path / "candidates.tsv"
    tsv.write_text(
        "chromosome\tposition\tref\talt\n"
        "chr8\t127401060\tG\tT\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        cli,
        ["--config", str(cfg), "agent", "run", str(tsv), "--dry-run"],
    )
    assert result.exit_code == 0
    assert "Top N: 42" in result.output
