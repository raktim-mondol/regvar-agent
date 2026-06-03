"""User-level TOML configuration for regvar.

Two locations are checked, in order:

    1. ``$REGVAR_CONFIG`` (explicit path, if set)
    2. ``./regvar.toml``  (project-local)
    3. ``~/.config/regvar/config.toml``  (XDG-style user config)

The first file that exists wins — files are *not* merged. Missing files are
silent: ``get_config()`` always returns a dict (empty if nothing was found).

Precedence for any given setting is: CLI flag > env var > config file > default.
The config layer therefore only provides *defaults* that the CLI can fall back
to when a flag is omitted.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

# XDG-style default location.
DEFAULT_CONFIG_PATHS: tuple[Path, ...] = (
    Path("./regvar.toml"),
    Path.home() / ".config" / "regvar" / "config.toml",
)


def _load_toml(path: Path) -> dict[str, Any]:
    """Load a TOML file, preferring stdlib tomllib (3.11+) and falling back to
    the optional ``tomli`` backport for 3.10."""
    try:
        import tomllib  # type: ignore[import-not-found]
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ImportError as exc:
            raise ImportError(
                "TOML config requires Python 3.11+ or the 'tomli' package. "
                "Install with: pip install tomli"
            ) from exc
    with path.open("rb") as fh:
        return tomllib.load(fh)


def _resolve_config_path(explicit: str | Path | None = None) -> Path | None:
    """Return the first config path that exists, or None."""
    if explicit is not None:
        p = Path(explicit)
        return p if p.is_file() else None

    env = os.environ.get("REGVAR_CONFIG")
    if env:
        p = Path(env)
        if p.is_file():
            return p

    for candidate in DEFAULT_CONFIG_PATHS:
        if candidate.is_file():
            return candidate
    return None


def get_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load the regvar TOML config and return it as a plain dict.

    Returns an empty dict if no config file is found. Raises ``ValueError`` if
    an explicit path is given but does not exist, or if the file is not valid
    TOML.
    """
    if path is not None:
        p = Path(path)
        if not p.is_file():
            raise ValueError(f"Config file not found: {p}")
        try:
            return _load_toml(p)
        except Exception as exc:
            raise ValueError(f"Failed to parse config {p}: {exc}") from exc

    resolved = _resolve_config_path()
    if resolved is None:
        return {}
    try:
        return _load_toml(resolved)
    except Exception as exc:
        # A discovered config file that fails to parse is worth reporting
        # loudly — silent fallback would mask real bugs.
        raise ValueError(f"Failed to parse config {resolved}: {exc}") from exc


def get_default(section: str, key: str, *, fallback: Any = None, path: str | Path | None = None) -> Any:
    """Fetch a single value from ``[section] key`` of the config, or *fallback*.

    Useful at the CLI layer to substitute defaults, e.g.::

        model = model or get_default("defaults", "model", fallback="deepseek-v4-pro")
    """
    cfg = get_config(path)
    return cfg.get(section, {}).get(key, fallback)
