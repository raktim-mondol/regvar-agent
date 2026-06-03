"""Tests for the TOML config layer (regvar/config.py)."""

import pytest

from regvar.config import get_config, get_default, _resolve_config_path


def test_get_config_empty_when_no_file(tmp_path, monkeypatch):
    """No config anywhere → empty dict (silent fallback)."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("REGVAR_CONFIG", raising=False)
    # Avoid hitting the real XDG location
    monkeypatch.setattr(
        "regvar.config.DEFAULT_CONFIG_PATHS",
        (tmp_path / "nonexistent.toml",),
    )
    assert get_config() == {}


def test_get_config_reads_explicit_path(tmp_path):
    cfg_path = tmp_path / "custom.toml"
    cfg_path.write_text(
        '[defaults]\nassays = ["ATAC-seq", "RNA-seq"]\ntop_n = 15\n',
        encoding="utf-8",
    )
    cfg = get_config(cfg_path)
    assert cfg["defaults"]["assays"] == ["ATAC-seq", "RNA-seq"]
    assert cfg["defaults"]["top_n"] == 15


def test_get_config_missing_explicit_path_raises(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        get_config(tmp_path / "missing.toml")


def test_get_config_invalid_toml_raises(tmp_path):
    bad = tmp_path / "bad.toml"
    bad.write_text("this is [not valid toml\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Failed to parse"):
        get_config(bad)


def test_resolve_config_path_env_var(tmp_path, monkeypatch):
    cfg = tmp_path / "env.toml"
    cfg.write_text("[defaults]\n", encoding="utf-8")
    monkeypatch.setenv("REGVAR_CONFIG", str(cfg))
    assert _resolve_config_path() == cfg


def test_resolve_config_path_local_over_xdg(tmp_path, monkeypatch):
    local = tmp_path / "regvar.toml"
    local.write_text("[defaults]\n", encoding="utf-8")
    xdg = tmp_path / "xdg.toml"
    xdg.write_text("[defaults]\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("REGVAR_CONFIG", raising=False)
    import regvar.config as cfg_mod
    monkeypatch.setattr(cfg_mod, "DEFAULT_CONFIG_PATHS", (local, xdg))
    assert _resolve_config_path() == local


def test_get_default_fallback():
    assert get_default("nonexistent", "key", fallback=42) == 42


def test_get_default_from_config(tmp_path):
    p = tmp_path / "c.toml"
    p.write_text(
        '[defaults]\nmodel = "some-model"\n', encoding="utf-8"
    )
    assert get_default("defaults", "model", fallback="x", path=p) == "some-model"
    assert get_default("defaults", "missing", fallback="x", path=p) == "x"
