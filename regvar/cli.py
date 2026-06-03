"""regvar CLI — installable command-line tool for regulatory variant triage.

After ``pip install -e .`` (or ``pip install regvar``) the ``regvar`` command
is available in the shell.  All commands also work via ``python -m regvar``.

Command hierarchy
-----------------
::

    regvar
    ├── score  <chr> <pos> <ref> <alt>   # score a single variant
    ├── plot
    │   └── effects <chr> <pos> <ref> <alt>  # barplot of variant effects
    ├── agent
    │   ├── run  <candidates.tsv>        # full triage loop (= old `regvar run`)
    │   └── chat                         # interactive REPL with the agent
    ├── mofa
    │   ├── build <scores.tsv> <genotypes.tsv>  # build MOFA+ view matrices
    │   ├── run   <views.h5ad>           # train MOFA+ model
    │   └── compare <views.h5ad>         # compare with vs. without variant views
    ├── assays                           # list supported assays
    ├── tui                              # Textual TUI
    └── run   <candidates.tsv>          # ← backward-compat alias → agent run
"""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import click
from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.theme import Theme

# ---------------------------------------------------------------------------
# Theme & consoles
# ---------------------------------------------------------------------------

REGVAR_THEME = Theme({
    "info":    "cyan",
    "success": "bold green",
    "warning": "bold yellow",
    "error":   "bold red",
    "variant": "bold magenta",
    "assay":   "bold blue",
    "score":   "bold white",
    "dim":     "dim",
})

console        = Console(theme=REGVAR_THEME, stderr=True)
stdout_console = Console(theme=REGVAR_THEME)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _banner() -> None:
    console.print(
        Panel(
            "[bold cyan]regvar[/bold cyan] [dim]·[/dim] "
            "AlphaGenome-powered regulatory variant triage",
            border_style="cyan",
            padding=(0, 1),
        )
    )


def _save_markdown(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    console.print(f"  [success]✓[/success] Report saved to [bold]{path}[/bold]")


def _save_tsv(records: list[dict], path: Path) -> None:
    from .variants import write_scores_tsv

    if not records:
        console.print("  [warning]⚠[/warning] No records to save.", style="warning")
        return
    write_scores_tsv(records, path)
    console.print(f"  [success]✓[/success] Scores saved to [bold]{path}[/bold]")


def _save_png(records: list[dict], title: str, path: Path, top_n: int = 20) -> None:
    """Render and save a barplot via regvar.plot."""
    try:
        from .plot import plot_effects
    except ImportError as exc:
        console.print(
            "[error]Error:[/error] matplotlib is required for PNG output. "
            "Run: [bold]pip install matplotlib[/bold]"
        )
        raise SystemExit(1) from exc

    plot_effects(records, title=title, output_path=path, top_n=top_n)
    console.print(f"  [success]✓[/success] Plot saved to [bold]{path}[/bold]")


def _print_variant_table(candidates) -> None:
    table = Table(
        title="Candidate Variants",
        title_style="bold cyan",
        border_style="dim",
        show_lines=False,
        pad_edge=False,
    )
    table.add_column("#", style="dim", width=3, justify="right")
    table.add_column("Variant", style="variant")
    table.add_column("Region", style="assay")
    table.add_column("Note", style="dim", max_width=60)
    for i, c in enumerate(candidates, 1):
        note_parts = [p for p in (c.note or "", _ann_note(c.vcf_id)) if p]
        note = " | ".join(note_parts) if note_parts else "—"
        table.add_row(str(i), c.vcf_id, c.region_id or "—", note)
    console.print(table)
    console.print()


def _print_assay_table(assays: dict[str, str]) -> None:
    table = Table(
        title="Supported Assays",
        title_style="bold cyan",
        border_style="dim",
    )
    table.add_column("Assay Label", style="assay")
    table.add_column("AlphaGenome Output", style="score")
    table.add_column("Wet-Lab Readout", style="dim")
    readout_map = {
        "ATAC":         "Chromatin accessibility",
        "DNASE":        "Chromatin accessibility (alt)",
        "RNA_SEQ":      "Gene expression",
        "CHIP_HISTONE": "H3K4me3/me1, H3K27ac, H3K27me3",
        "CHIP_TF":      "Transcription-factor binding",
        "CONTACT_MAPS": "Enhancer–promoter looping",
        "SPLICE_SITES": "Splice-site usage",
    }
    for label, output in assays.items():
        table.add_row(label, output, readout_map.get(output, ""))
    console.print(table)


def _print_scores_table(records: list[dict]) -> None:
    if not records:
        console.print("  [warning]No scores returned.[/warning]")
        return
    table = Table(
        title="Top Predicted Effects",
        title_style="bold cyan",
        border_style="dim",
    )
    table.add_column("#", style="dim", width=3, justify="right")
    display_cols = [
        ("gene_name",       "Gene",        "bold white"),
        ("assay",           "Assay",       "assay"),
        ("output_type",     "Output Type", "dim"),
        ("biosample_name",  "Biosample",   "dim"),
        ("raw_score",       "Raw Score",   "score"),
        ("quantile_score",  "Quantile",    "score"),
        ("abs_effect",      "|Effect|",    "bold yellow"),
    ]
    active_cols = [
        (key, hdr, sty)
        for key, hdr, sty in display_cols
        if records and key in records[0]
    ]
    for _, hdr, sty in active_cols:
        table.add_column(hdr, style=sty)
    for i, rec in enumerate(records, 1):
        row = [str(i)]
        for key, _, _ in active_cols:
            val = rec.get(key, "")
            if isinstance(val, float):
                row.append(f"{val:+.4f}" if key != "abs_effect" else f"{val:.4f}")
            else:
                row.append(str(val))
        table.add_row(*row)
    console.print(table)


def _read_candidates(path: str):
    """Auto-detect VCF / TSV by extension and return a list of CandidateVariant."""
    from pathlib import Path as _Path

    suffix = "".join(_Path(path).suffixes).lower()
    if suffix in (".vcf", ".vcf.gz"):
        from .variants import read_candidates_vcf
        return read_candidates_vcf(path)
    from .variants import read_candidates_tsv
    return read_candidates_tsv(path)


def _get_config(ctx: click.Context | None = None) -> dict:
    """Return the resolved TOML config dict (cached on the Click context).

    Safe to call from any subcommand: if no config file is present the result
    is an empty dict, so ``cfg.get("section", {}).get("key", fallback)`` works
    without special-casing.
    """
    path: str | None = None
    if ctx is not None:
        obj = ctx.obj or {}
        if "config_cache" in obj:
            return obj["config_cache"]
        path = obj.get("config_path")
    from .config import get_config
    try:
        cfg = get_config(path)
    except ValueError as exc:
        console.print(f"[warning]⚠[/warning] {exc}", style="warning")
        cfg = {}
    if ctx is not None:
        ctx.ensure_object(dict)
        ctx.obj["config_cache"] = cfg
    return cfg


def _cfg_value(ctx: click.Context | None, section: str, key: str, fallback=None):
    """Look up ``[section] key`` in the TOML config, or return *fallback*."""
    return _get_config(ctx).get(section, {}).get(key, fallback)


def _apply_annotations(
    candidates: list,
    annotate_dir: str,
    *,
    quiet: bool = False,
) -> None:
    """Auto-detect reference files in *annotate_dir* and attach annotations.

    Recognised filenames (any subset): ``genes.gtf[.gz]``, ``peaks.bed[.gz]``,
    ``dbsnp.vcf.gz``. Missing files are silently skipped.
    """
    from pathlib import Path as _Path

    from .annotation import annotate_variants

    d = _Path(annotate_dir)
    gtf = next((p for p in (d / "genes.gtf", d / "genes.gtf.gz") if p.is_file()), None)
    bed = next((p for p in (d / "peaks.bed", d / "peaks.bed.gz") if p.is_file()), None)
    dbsnp = d / "dbsnp.vcf.gz" if (d / "dbsnp.vcf.gz").is_file() else None

    if gtf is None and bed is None and dbsnp is None:
        if not quiet:
            console.print(
                f"  [warning]⚠[/warning] --annotate dir {d} has no "
                f"genes.gtf(.gz), peaks.bed(.gz), or dbsnp.vcf.gz; skipping.",
                style="warning",
            )
        return

    if not quiet:
        found = [p.name for p in (gtf, bed, dbsnp) if p is not None]
        console.print(f"  Annotating with: {', '.join(found)}")

    try:
        annotate_variants(
            candidates,
            gtf_path=gtf,
            bed_path=bed,
            dbsnp_vcf_path=dbsnp,
        )
        annotated = sum(1 for c in candidates if _ann_note(c.vcf_id))
        if not quiet:
            console.print(
                f"  [success]✓[/success] Annotated "
                f"[bold]{annotated}/{len(candidates)}[/bold] variant(s)."
            )
    except ImportError as exc:
        if not quiet:
            console.print(f"  [warning]⚠[/warning] {exc}", style="warning")


def _ann_note(vcf_id: str) -> str:
    """Return the annotation snippet for *vcf_id*, or "" if unannotated."""
    from .annotation import get_annotation
    ann = get_annotation(vcf_id)
    return ann.to_note() if ann is not None else ""


def _parse_assays(raw: str | None) -> list[str] | None:
    """Parse the --assays CLI value.

    - ``None`` → return ``None`` (use defaults)
    - ``"all"`` → return all keys from ``ASSAY_TO_OUTPUT``
    - Comma-separated string → split and strip
    """
    if raw is None:
        return None
    stripped = raw.strip()
    if stripped.lower() == "all":
        from .alphagenome_client import ASSAY_TO_OUTPUT
        return list(ASSAY_TO_OUTPUT.keys())
    return [a.strip() for a in stripped.split(",") if a.strip()]


def _serialize_messages(messages: list[dict]) -> list[dict]:
    """Convert chat message history to JSON-safe dicts.

    OpenAI message objects (tool_calls, function) get dumped via pydantic's
    ``model_dump``; plain dicts pass through untouched.
    """
    out: list[dict] = []
    for m in messages:
        if hasattr(m, "model_dump"):
            try:
                out.append(m.model_dump(exclude_unset=True))
                continue
            except Exception:
                pass
        if hasattr(m, "dict"):
            try:
                out.append(m.dict(exclude_unset=True))
                continue
            except Exception:
                pass
        out.append(dict(m))
    return out


def _save_session(messages: list[dict], path: Path) -> None:
    """Persist message history to JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"messages": _serialize_messages(messages)}, indent=2, default=str),
        encoding="utf-8",
    )


def _load_session(path: Path) -> list[dict]:
    """Load message history from a JSON session file."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if "messages" not in data:
        raise ValueError(f"{path}: session file must contain a 'messages' list")
    msgs = data["messages"]
    if not isinstance(msgs, list):
        raise ValueError(f"{path}: session file 'messages' must be a list")
    return [dict(m) for m in msgs]


def _render_conversation_markdown(messages: list[dict]) -> str:
    """Render the conversation history as a markdown transcript."""
    lines = ["# Agent Chat Session", ""]
    for m in messages:
        if hasattr(m, "model_dump"):
            try:
                m = m.model_dump(exclude_unset=True)
            except Exception:
                m = dict(m)
        role = m.get("role", "unknown")
        content = m.get("content") or ""
        if role == "system":
            lines.append("## System\n")
            lines.append(str(content))
            lines.append("")
        elif role == "user":
            lines.append("## User\n")
            lines.append(str(content))
            lines.append("")
        elif role == "assistant":
            lines.append("## Assistant\n")
            if content:
                lines.append(str(content))
            tool_calls = m.get("tool_calls") or []
            for tc in tool_calls:
                fn = tc.get("function", {}) if isinstance(tc, dict) else {}
                name = fn.get("name", "?")
                args = fn.get("arguments", "{}")
                lines.append(f"\n*tool call:* `{name}({args})`\n")
            lines.append("")
        elif role == "tool":
            name = m.get("name", m.get("tool_call_id", "tool"))
            lines.append(f"### tool result ({name})\n")
            lines.append(f"```\n{content}\n```")
            lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------

@click.group(invoke_without_command=True)
@click.option(
    "--config",
    "config_path",
    type=click.Path(dir_okay=False),
    default=None,
    envvar="REGVAR_CONFIG",
    help="Path to a TOML config file (default: ./regvar.toml or "
         "~/.config/regvar/config.toml).",
)
@click.pass_context
@click.version_option(version="0.1.0", prog_name="regvar")
def cli(ctx: click.Context, config_path: str | None) -> None:
    """regvar — AlphaGenome-powered regulatory variant triage.

    Score candidate non-coding variants with AlphaGenome and triage them
    using a DeepSeek reasoning agent.

    \b
    Quick start:
        regvar score chr8 127401060 G T
        regvar agent run examples/candidate_variants.tsv
        regvar plot effects chr8 127401060 G T --output effects.png
        regvar assays
    """
    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config_path
    # Populate Click's default_map from TOML so omitted CLI flags fall back to
    # config values. Click auto-propagates nested sections (e.g. agent.run)
    # when the top-level map contains a group name. Top-level [defaults] values
    # must be duplicated into every leaf subcommand entry to actually reach
    # them — Click does NOT cascade defaults across nesting.
    try:
        cfg = _get_config(ctx)
    except Exception:
        cfg = {}
    if cfg and not ctx.default_map:
        defaults = dict(cfg.get("defaults", {}))
        shared = {k.replace("-", "_"): v for k, v in defaults.items()}
        top_map: dict = {}

        # Per-subcommand sections: [score], [plot_effects], [agent_run], etc.
        subcommand_entries = {
            "score":        (None,   "score"),
            "plot_effects": ("plot", "effects"),
            "agent_run":    ("agent", "run"),
            "agent_chat":   ("agent", "chat"),
        }
        for section_name, (parent_group, cmd_name) in subcommand_entries.items():
            section = cfg.get(section_name, {})
            if not section and not shared:
                continue
            merged = dict(shared)
            merged.update({k.replace("-", "_"): v for k, v in section.items()})
            if parent_group is None:
                top_map[cmd_name] = merged
            else:
                top_map.setdefault(parent_group, {})[cmd_name] = merged

        # Also accept [agent.run] / [agent.chat] TOML sections directly.
        agent_section = cfg.get("agent", {})
        if isinstance(agent_section, dict):
            for sub in ("run", "chat"):
                sub_cfg = agent_section.get(sub, {})
                if isinstance(sub_cfg, dict) and sub_cfg:
                    merged = dict(shared)
                    existing = top_map.get("agent", {}).get(sub, {})
                    merged.update(existing)
                    merged.update({k.replace("-", "_"): v for k, v in sub_cfg.items()})
                    top_map.setdefault("agent", {})[sub] = merged

        ctx.default_map = top_map
    if ctx.invoked_subcommand is None:
        _banner()
        click.echo(ctx.get_help())


# ---------------------------------------------------------------------------
# regvar score
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("chromosome")
@click.argument("position", type=int)
@click.argument("ref")
@click.argument("alt")
@click.option(
    "-a", "--assays",
    default=None,
    help="Comma-separated assays, or 'all' for every supported assay. Default: ATAC-seq, RNA-seq, ChIP-seq histone, Hi-C / pcHi-C.",
)
@click.option(
    "-t", "--tissue",
    default=None,
    help="Comma-separated UBERON/CL ontology terms. Default: UBERON:0002367 (prostate gland).",
)
@click.option(
    "-n", "--top-n",
    type=int,
    default=10,
    show_default=True,
    help="Number of top effects to show.",
)
@click.option(
    "--output-tsv",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save results as TSV.",
)
@click.option(
    "--output-png",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save a barplot of top effects as a PNG file.",
)
@click.option(
    "--output-json",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save a structured JSON with metadata and scores.",
)
@click.option("--json", "as_json", is_flag=True, help="Output raw JSON instead of a table.")
def score(
    chromosome: str,
    position: int,
    ref: str,
    alt: str,
    assays: str | None,
    tissue: str | None,
    top_n: int,
    output_tsv: str | None,
    output_png: str | None,
    output_json: str | None,
    as_json: bool,
) -> None:
    """Score a single variant directly (no agent loop).

    \b
    Examples:
        regvar score chr8 127401060 G T
        regvar score chr8 127401060 G T --assays ATAC-seq,RNA-seq --top-n 20
        regvar score chr8 127401060 G T --output-png effects.png
    """
    _banner()
    console.print(f"  Scoring [variant]{chromosome}:{position}:{ref}>{alt}[/variant]\n")

    assay_list  = _parse_assays(assays)
    tissue_list = [t.strip() for t in tissue.split(",")] if tissue else None

    from .tools import tool_score_regulatory_variant

    kwargs: dict = {
        "chromosome": chromosome,
        "position": position,
        "ref": ref,
        "alt": alt,
        "top_n": top_n,
    }
    if assay_list:
        kwargs["assays"] = assay_list
    if tissue_list:
        kwargs["ontology_terms"] = tissue_list

    with console.status("Querying AlphaGenome…", spinner="dots"):
        result = tool_score_regulatory_variant(**kwargs)

    if as_json:
        stdout_console.print_json(json.dumps(result, default=str))
    else:
        console.print(
            f"  Variant: [variant]{result['variant']}[/variant]  |  "
            f"Assays: {', '.join(result['assays_scored'])}  |  "
            f"Total scores: {result['n_total_scores']}\n"
        )
        _print_scores_table(result["top_effects"])

    if output_tsv:
        _save_tsv(result["top_effects"], Path(output_tsv))

    if output_png:
        title = f"{chromosome}:{position} {ref}>{alt}"
        _save_png(result["top_effects"], title=title, path=Path(output_png), top_n=top_n)

    if output_json:
        json_output = {
            "metadata": {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "assays": result.get("assays_scored", []),
                "tissue": tissue_list,
                "top_n": top_n,
            },
            "variants": [{
                "chromosome": chromosome,
                "position": position,
                "ref": ref,
                "alt": alt,
            }],
            "scores": result.get("top_effects", []),
        }
        json_path = Path(output_json)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(json_output, indent=2, default=str), encoding="utf-8")
        console.print(f"  [success]✓[/success] JSON saved to [bold]{json_path}[/bold]")


# ---------------------------------------------------------------------------
# regvar plot  (command group)
# ---------------------------------------------------------------------------

@cli.group()
def plot() -> None:
    """Visualise variant effect predictions.

    \b
    Subcommands:
        effects   Horizontal barplot of top predicted effects
    """


@plot.command("effects")
@click.argument("chromosome")
@click.argument("position", type=int)
@click.argument("ref")
@click.argument("alt")
@click.option(
    "-a", "--assays",
    default=None,
    help="Comma-separated assays to score, or 'all'. Default: ATAC-seq, RNA-seq, ChIP-seq histone, Hi-C / pcHi-C.",
)
@click.option(
    "-t", "--tissue",
    default=None,
    help="Comma-separated UBERON/CL ontology terms. Default: UBERON:0002367.",
)
@click.option(
    "-n", "--top-n",
    type=int,
    default=20,
    show_default=True,
    help="Number of top effects to plot.",
)
@click.option(
    "--from-tsv",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Load scores from a TSV produced by 'regvar score --output-tsv' "
         "instead of calling the API. Positional variant args are used only for the title.",
)
@click.option(
    "--assay-filter",
    default=None,
    help="Keep only rows where 'assay' matches this string (e.g. 'ATAC-seq').",
)
@click.option(
    "-o", "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save plot to this file (.png / .svg / .pdf). Required when --show is not set.",
)
@click.option(
    "--show",
    is_flag=True,
    default=False,
    help="Open an interactive matplotlib window (requires a display).",
)
def plot_effects(
    chromosome: str,
    position: int,
    ref: str,
    alt: str,
    assays: str | None,
    tissue: str | None,
    top_n: int,
    from_tsv: str | None,
    assay_filter: str | None,
    output: str | None,
    show: bool,
) -> None:
    """Barplot of top predicted effects for a single variant.

    \b
    Examples:
        regvar plot effects chr8 127401060 G T --output effects.png
        regvar plot effects chr8 127401060 G T --show
        regvar plot effects chr8 127401060 G T \\
            --from-tsv scores.tsv --assay-filter ATAC-seq --output atac.png
    """
    if not output and not show:
        raise click.UsageError(
            "Specify at least one of --output <file> or --show."
        )

    _banner()
    title = f"{chromosome}:{position} {ref}>{alt}"
    console.print(f"  Plotting effects for [variant]{title}[/variant]\n")

    # ── load data ────────────────────────────────────────────────────────────
    if from_tsv:
        import csv
        with open(from_tsv, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            records: list[dict] = list(reader)
        # Coerce numeric fields
        for rec in records:
            for field in ("raw_score", "quantile_score", "abs_effect"):
                try:
                    rec[field] = float(rec[field])
                except (KeyError, ValueError, TypeError):
                    pass
    else:
        assay_list  = _parse_assays(assays)
        tissue_list = [t.strip() for t in tissue.split(",")] if tissue else None

        from .tools import tool_score_regulatory_variant

        kwargs: dict = {
            "chromosome": chromosome,
            "position": position,
            "ref": ref,
            "alt": alt,
            "top_n": top_n,
        }
        if assay_list:
            kwargs["assays"] = assay_list
        if tissue_list:
            kwargs["ontology_terms"] = tissue_list

        with console.status("Querying AlphaGenome…", spinner="dots"):
            result = tool_score_regulatory_variant(**kwargs)
        records = result["top_effects"]

    # ── optional assay filter ────────────────────────────────────────────────
    if assay_filter:
        records = [r for r in records if str(r.get("assay", "")) == assay_filter]
        if not records:
            console.print(
                f"  [warning]⚠[/warning] No records found for assay filter "
                f"[bold]{assay_filter!r}[/bold]."
            )
            return

    # ── plot ─────────────────────────────────────────────────────────────────
    try:
        from .plot import plot_effects as _plot
    except ImportError as exc:
        console.print(
            "[error]Error:[/error] matplotlib is required. "
            "Run: [bold]pip install matplotlib[/bold]"
        )
        raise SystemExit(1) from exc

    out_path = Path(output) if output else None
    with console.status("Rendering plot…", spinner="dots"):
        _plot(records, title=title, output_path=out_path, show=show, top_n=top_n)

    if output:
        console.print(f"  [success]✓[/success] Plot saved to [bold]{output}[/bold]")
    if show:
        console.print("  [dim]Plot window opened.[/dim]")


# ---------------------------------------------------------------------------
# regvar agent  (command group)
# ---------------------------------------------------------------------------

@cli.group()
def agent() -> None:
    """Run the DeepSeek reasoning agent for variant triage.

    \b
    Subcommands:
        run    Full triage loop over a TSV of candidate variants
        chat   Interactive REPL — ask follow-up questions in plain English
    """


@agent.command("run")
@click.argument("candidates_tsv", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-a", "--assays",
    default=None,
    help="Comma-separated assays to request, or 'all' for every supported assay.",
)
@click.option(
    "-t", "--tissue",
    default=None,
    help="Comma-separated UBERON/CL ontology terms.",
)
@click.option(
    "-n", "--top-n",
    type=int,
    default=10,
    show_default=True,
    help="Number of top effects to return per variant.",
)
@click.option(
    "-m", "--model",
    default="deepseek-v4-pro",
    show_default=True,
    help="LLM model for the agent loop.",
)
@click.option(
    "-o", "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save the agent's markdown report to this file.",
)
@click.option(
    "--output-tsv",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save ranked scores as a TSV file.",
)
@click.option(
    "--output-json",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save a structured JSON with metadata, scores, and the agent report.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Validate input and show what would be scored — no API calls.",
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    default=False,
    help="Show tool calls and extra debug info.",
)
@click.option(
    "-q", "--quiet",
    is_flag=True,
    default=False,
    help="Suppress all progress output; print only the final report.",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Proceed even when the variant count exceeds the safety limit.",
)
@click.option(
    "--max-concurrent",
    type=int,
    default=None,
    help="Parallel scoring workers. Overrides [defaults] max_concurrent in config. "
         "Default (1) runs serially; higher values dispatch via a thread pool.",
)
@click.option(
    "--annotate",
    type=click.Path(exists=True, file_okay=False),
    default=None,
    help="Directory with reference files for variant annotation: genes.gtf(.gz), "
         "peaks.bed(.gz), dbsnp.vcf.gz. Any subset is accepted; annotations are "
         "appended to the agent's system prompt as context per variant.",
)
def agent_run(
    candidates_tsv: str,
    assays: str | None,
    tissue: str | None,
    top_n: int,
    model: str,
    output: str | None,
    output_tsv: str | None,
    output_json: str | None,
    dry_run: bool,
    verbose: bool,
    quiet: bool,
    force: bool,
    max_concurrent: int | None,
    annotate: str | None,
) -> None:
    """Run the agent to triage candidate regulatory variants.

    Reads variants from CANDIDATES_TSV (.tsv) or a VCF (.vcf / .vcf.gz),
    scores each with AlphaGenome, and uses a DeepSeek agent to produce a
    ranked analysis with regulatory mechanisms and a wet-lab validation plan.

    \b
    Examples:
        regvar agent run examples/candidate_variants.tsv
        regvar agent run gwas_hits.vcf.gz --dry-run
        regvar agent run examples/candidate_variants.tsv -o report.md --output-tsv scores.tsv
        regvar agent run examples/candidate_variants.tsv --output-json results.json
        regvar agent run examples/candidate_variants.tsv --max-concurrent 4
    """
    _run_agent_core(
        candidates_tsv=candidates_tsv,
        assays=assays,
        tissue=tissue,
        top_n=top_n,
        model=model,
        output=output,
        output_tsv=output_tsv,
        output_json=output_json,
        dry_run=dry_run,
        verbose=verbose,
        quiet=quiet,
        force=force,
        max_concurrent=max_concurrent,
        annotate=annotate,
    )


@agent.command("chat")
@click.option(
    "-m", "--model",
    default="deepseek-v4-pro",
    show_default=True,
    help="LLM model for the chat session.",
)
@click.option(
    "-a", "--assays",
    default=None,
    help="Preferred assays (comma-separated, or 'all'). Passed as context to the agent.",
)
@click.option(
    "-t", "--tissue",
    default=None,
    help="Preferred tissue ontology terms (comma-separated).",
)
@click.option(
    "-n", "--top-n",
    type=int,
    default=10,
    show_default=True,
    help="Default top-N effects to return per scoring call.",
)
@click.option(
    "-o", "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help="On exit, save the full conversation transcript as markdown.",
)
@click.option(
    "--reasoning-effort",
    type=click.Choice(["low", "medium", "high", "off"], case_sensitive=False),
    default="high",
    show_default=True,
    help="Reasoning effort level for the model. 'off' disables extended thinking.",
)
@click.option(
    "--save-session",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save the message history as JSON on exit (for --load-session).",
)
@click.option(
    "--load-session",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Resume a previously saved JSON session (appends to existing history).",
)
def agent_chat(
    model: str,
    assays: str | None,
    tissue: str | None,
    top_n: int,
    output: str | None,
    reasoning_effort: str,
    save_session: str | None,
    load_session: str | None,
) -> None:
    """Interactive REPL — ask the agent questions in plain English.

    Type a variant (e.g. "score chr8:127401060 G>T") or a free-form question.
    The agent has access to AlphaGenome scoring tools and will call them on
    demand. Type 'exit' or press Ctrl-C to quit.

    \b
    Examples:
        regvar agent chat
        regvar agent chat --tissue UBERON:0002367 --top-n 15
        regvar agent chat --reasoning-effort low --save-session chat.json
        regvar agent chat --load-session chat.json --output transcript.md
    """
    _banner()
    console.print(
        Panel(
            "[bold cyan]Agent Chat[/bold cyan]  [dim]·[/dim]  "
            "Interactive variant triage session\n"
            "[dim]Type a question or variant to score. "
            "Enter [bold]exit[/bold] or press Ctrl-C to quit.[/dim]",
            border_style="cyan",
            padding=(0, 1),
        )
    )

    # Pre-flight
    if not os.environ.get("DEEPSEEK_API_KEY"):
        console.print(
            "[error]Error:[/error] DEEPSEEK_API_KEY not set. "
            "Export it or add it to .env"
        )
        raise SystemExit(1)

    try:
        from openai import OpenAI
    except ImportError:
        console.print("[error]Error:[/error] pip install openai")
        raise SystemExit(1)

    from .tools import OPENAI_TOOL_SCHEMAS
    from .agent import DEEPSEEK_BASE_URL, _build_system_prompt, MAX_TURNS, run_agent_turn

    assay_list  = _parse_assays(assays)
    tissue_list = [t.strip() for t in tissue.split(",")] if tissue else None

    system_prompt = _build_system_prompt(assay_list, tissue_list, top_n)
    messages: list[dict] = [{"role": "system", "content": system_prompt}]

    if load_session:
        try:
            prior = _load_session(Path(load_session))
            # Keep current system prompt authoritative; skip any prior system msg.
            messages.extend(m for m in prior if m.get("role") != "system")
            console.print(
                f"  [success]✓[/success] Resumed session from [bold]{load_session}[/bold] "
                f"([dim]{len(prior)} message(s)[/dim])"
            )
        except (ValueError, json.JSONDecodeError, OSError) as exc:
            console.print(f"[error]Error loading session:[/error] {exc}")
            raise SystemExit(1) from exc

    effort = None if reasoning_effort.lower() == "off" else reasoning_effort.lower()

    client = OpenAI(
        api_key=os.environ["DEEPSEEK_API_KEY"],
        base_url=DEEPSEEK_BASE_URL,
    )

    console.print()

    try:
        while True:
            # ── user input ────────────────────────────────────────────────
            try:
                user_input = click.prompt(
                    click.style("you", fg="cyan", bold=True),
                    prompt_suffix=" › ",
                )
            except (EOFError, KeyboardInterrupt):
                console.print("\n  [dim]Session ended.[/dim]")
                break

            if user_input.strip().lower() in {"exit", "quit", "q"}:
                console.print("  [dim]Session ended.[/dim]")
                break
            if not user_input.strip():
                continue

            messages.append({"role": "user", "content": user_input})

            # ── agent loop ────────────────────────────────────────────────
            with console.status("Thinking…", spinner="dots"):
                answer = run_agent_turn(
                    messages=messages,
                    client=client,
                    model=model,
                    tools=OPENAI_TOOL_SCHEMAS,
                    max_tokens=4096,
                    reasoning_effort=effort,
                    max_turns=MAX_TURNS,
                )

            # ── render response ───────────────────────────────────────────
            console.print()
            console.rule("[dim]agent[/dim]", style="dim")
            stdout_console.print(Markdown(answer))
            console.rule(style="dim")
            console.print()

    except KeyboardInterrupt:
        console.print("\n  [dim]Session ended.[/dim]")
    finally:
        if output:
            md = _render_conversation_markdown(messages)
            _save_markdown(md, Path(output))
        if save_session:
            _save_session(messages, Path(save_session))
            console.print(
                f"  [success]✓[/success] Session saved to [bold]{save_session}[/bold]"
            )


# ---------------------------------------------------------------------------
# Shared agent-run implementation (used by both `agent run` and the
# backward-compat top-level `run` command)
# ---------------------------------------------------------------------------

def _run_agent_core(
    candidates_tsv: str,
    assays: str | None,
    tissue: str | None,
    top_n: int,
    model: str,
    output: str | None,
    output_tsv: str | None,
    output_json: str | None,
    dry_run: bool,
    verbose: bool,
    quiet: bool,
    force: bool = False,
    max_concurrent: int | None = None,
    annotate: str | None = None,
) -> None:
    if not quiet:
        _banner()

    # Apply parallelism override to the shared client singleton, if given.
    if max_concurrent is not None:
        from .tools import get_client
        get_client().config.max_concurrent = max(max_concurrent, 1)

    try:
        candidates = _read_candidates(candidates_tsv)
    except (ValueError, FileNotFoundError) as exc:
        console.print(f"[error]Error:[/error] {exc}")
        raise SystemExit(1) from exc

    if annotate:
        _apply_annotations(candidates, annotate, quiet=quiet)

    if not quiet:
        console.print(
            f"  Loaded [bold]{len(candidates)}[/bold] variant(s) from "
            f"[bold]{candidates_tsv}[/bold]\n"
        )
        _print_variant_table(candidates)

    from .variants import check_variant_count
    try:
        check_variant_count(candidates, force=force)
    except ValueError as exc:
        console.print(f"[error]Error:[/error] {exc}")
        raise SystemExit(1) from exc

    assay_list  = _parse_assays(assays)
    tissue_list = [t.strip() for t in tissue.split(",")] if tissue else None

    if dry_run:
        console.print(
            Panel(
                "[bold yellow]Dry run[/bold yellow] — no API calls will be made.",
                border_style="yellow",
            )
        )
        if assay_list:
            console.print(f"  Assays: {', '.join(assay_list)}")
        if tissue_list:
            console.print(f"  Tissue terms: {', '.join(tissue_list)}")
        console.print(f"  Top N: {top_n}")
        console.print(f"  Model: {model}")
        if max_concurrent is not None:
            console.print(f"  Max concurrent: {max_concurrent}")
        console.print()
        console.print("[success]Input validated. Ready to run.[/success]")
        return

    if not os.environ.get("DEEPSEEK_API_KEY"):
        console.print(
            "[error]Error:[/error] DEEPSEEK_API_KEY not set. "
            "Export it or add it to .env"
        )
        raise SystemExit(1)

    if not quiet:
        console.print("  Starting agent loop…\n", style="info")

    from .agent import run_agent
    from rich.progress import (
        BarColumn, MofNCompleteColumn, Progress, SpinnerColumn, TextColumn,
        TimeElapsedColumn,
    )

    tool_call_log: list[dict] = []
    scored_count = 0

    def _on_tool_call(name: str, args: dict, result: str) -> None:
        nonlocal scored_count
        tool_call_log.append({"tool": name, "args": args, "result": result})
        if name in ("score_regulatory_variant", "score_variants_batch"):
            if name == "score_variants_batch":
                scored_count += len(args.get("variants", []))
            else:
                scored_count += 1
        if quiet:
            return
        if verbose:
            console.print(
                f"  [dim]→[/dim] [assay]{name}[/assay]({json.dumps(args, default=str)})"
            )
        else:
            if name == "score_regulatory_variant":
                v = f"{args.get('chromosome', '?')}:{args.get('position', '?')}"
                console.print(f"  [dim]→[/dim] scoring [variant]{v}[/variant]")
            else:
                console.print(f"  [dim]→[/dim] [assay]{name}[/assay]")

    progress_columns = [
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
    ]

    start = time.time()
    if not quiet:
        with Progress(
            *progress_columns,
            console=console,
            transient=True,
        ) as progress:
            task_id = progress.add_task(
                "Scoring variants…",
                total=len(candidates),
            )
            _original_on_tool_call = _on_tool_call

            def _progress_on_tool_call(name: str, args: dict, result: str) -> None:
                _original_on_tool_call(name, args, result)
                if name in ("score_regulatory_variant", "score_variants_batch"):
                    if name == "score_variants_batch":
                        progress.advance(task_id, len(args.get("variants", [])))
                    else:
                        progress.advance(task_id, 1)
                    progress.update(task_id, description=f"Scored {scored_count}/{len(candidates)} variant(s)…")

            answer = run_agent(
                candidates,
                model=model,
                verbose=False,
                on_tool_call=_progress_on_tool_call,
                assay_hint=assay_list,
                tissue_hint=tissue_list,
                top_n_hint=top_n,
            )
    else:
        answer = run_agent(
            candidates,
            model=model,
            verbose=False,
            on_tool_call=_on_tool_call,
            assay_hint=assay_list,
            tissue_hint=tissue_list,
            top_n_hint=top_n,
        )
    elapsed = time.time() - start

    if not quiet:
        console.print()
        console.print(
            f"  [success]✓[/success] Agent finished in "
            f"[bold]{elapsed:.1f}s[/bold] "
            f"([dim]{len(tool_call_log)} tool call(s)[/dim])"
        )
        console.print()

    stdout_console.print(Markdown(answer))

    if output:
        _save_markdown(answer, Path(output))

    if output_tsv:
        score_records: list[dict] = []
        for entry in tool_call_log:
            if entry["tool"] in ("score_regulatory_variant", "score_variants_batch"):
                try:
                    parsed = (
                        json.loads(entry["result"])
                        if isinstance(entry["result"], str)
                        else entry["result"]
                    )
                    if entry["tool"] == "score_regulatory_variant":
                        for rec in parsed.get("top_effects", []):
                            rec["query_variant"] = parsed.get("variant", "")
                            score_records.append(rec)
                    else:
                        for vid, effects in parsed.get("effects_by_variant", {}).items():
                            for rec in effects:
                                rec["query_variant"] = vid
                                score_records.append(rec)
                except (json.JSONDecodeError, TypeError):
                    pass
        _save_tsv(score_records, Path(output_tsv))

    if output_json:
        variant_list = [
            {
                "chromosome": c.chromosome,
                "position": c.position,
                "ref": c.ref,
                "alt": c.alt,
                "region_id": c.region_id,
                "note": c.note,
            }
            for c in candidates
        ]
        scores: list[dict] = []
        for entry in tool_call_log:
            if entry["tool"] in ("score_regulatory_variant", "score_variants_batch"):
                try:
                    parsed = (
                        json.loads(entry["result"])
                        if isinstance(entry["result"], str)
                        else entry["result"]
                    )
                    if entry["tool"] == "score_regulatory_variant":
                        for rec in parsed.get("top_effects", []):
                            rec["query_variant"] = parsed.get("variant", "")
                            scores.append(rec)
                    else:
                        for vid, effects in parsed.get("effects_by_variant", {}).items():
                            for rec in effects:
                                rec["query_variant"] = vid
                                scores.append(rec)
                except (json.JSONDecodeError, TypeError):
                    pass
        json_output = {
            "metadata": {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "model": model,
                "assays": assay_list,
                "tissue": tissue_list,
                "top_n": top_n,
            },
            "variants": variant_list,
            "scores": scores,
            "report": answer,
        }
        json_path = Path(output_json)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(json_output, indent=2, default=str), encoding="utf-8")
        console.print(f"  [success]✓[/success] JSON saved to [bold]{json_path}[/bold]")


# ---------------------------------------------------------------------------
# regvar mofa  (command group)
# ---------------------------------------------------------------------------

@cli.group()
def mofa() -> None:
    """Build MOFA+ views from variant scores and train factor models.

    \b
    Subcommands:
        build     Convert variant scores + genotypes → MOFA+-ready matrices
        run       Train a MOFA+ model (requires mofapy2)
        compare   Train with vs. without variant views; report differences
    """


@mofa.command("build")
@click.argument("scores_tsv", type=click.Path(exists=True, dir_okay=False))
@click.argument("genotypes_tsv", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-a", "--assays",
    default=None,
    help="Comma-separated assay subset, or 'all'. Default: all assays in the scores.",
)
@click.option(
    "--by",
    type=click.Choice(["max_abs", "top_gene", "all"]),
    default="max_abs",
    show_default=True,
    help="Aggregation strategy when a variant has multiple rows per assay.",
)
@click.option(
    "-o", "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save views as AnnData .h5ad file (requires anndata).",
)
@click.option(
    "--output-csv-dir",
    type=click.Path(file_okay=False),
    default=None,
    help="Save each view as a separate CSV in this directory.",
)
def mofa_build(
    scores_tsv: str,
    genotypes_tsv: str,
    assays: str | None,
    by: str,
    output: str | None,
    output_csv_dir: str | None,
) -> None:
    """Build MOFA+-ready view matrices from variant scores and genotypes.

    SCORES_TSV is a TSV from 'regvar agent run --output-tsv'.
    GENOTYPES_TSV has sample_id + one column per variant (0/1/2/NaN).

    \b
    Examples:
        regvar mofa build scores.tsv genotypes.tsv -o views.h5ad
        regvar mofa build scores.tsv genotypes.tsv --by top_gene
        regvar mofa build scores.tsv genotypes.tsv --assays ATAC-seq,RNA-seq
        regvar mofa build scores.tsv genotypes.tsv --output-csv-dir views/
    """
    _banner()

    from .mofa_view import build_mofa_view, read_genotypes_tsv, scored_variants_from_tsv

    # ── load data ────────────────────────────────────────────────────────────
    with console.status("Loading scores…", spinner="dots"):
        try:
            scored = scored_variants_from_tsv(scores_tsv)
        except (ValueError, FileNotFoundError) as exc:
            console.print(f"[error]Error:[/error] {exc}")
            raise SystemExit(1) from exc

    console.print(
        f"  Loaded [bold]{len(scored)}[/bold] variant(s) from "
        f"[bold]{scores_tsv}[/bold]"
    )

    with console.status("Loading genotypes…", spinner="dots"):
        try:
            genotypes = read_genotypes_tsv(genotypes_tsv)
        except (ValueError, FileNotFoundError) as exc:
            console.print(f"[error]Error:[/error] {exc}")
            raise SystemExit(1) from exc

    console.print(
        f"  Loaded [bold]{len(genotypes)}[/bold] sample(s) and "
        f"[bold]{len(genotypes.columns)}[/bold] variant(s) from "
        f"[bold]{genotypes_tsv}[/bold]"
    )

    # ── build views ──────────────────────────────────────────────────────────
    assay_list = _parse_assays(assays)

    with console.status("Building MOFA+ views…", spinner="dots"):
        try:
            views = build_mofa_view(
                scored_variants=scored,
                genotypes=genotypes,
                assays=assay_list,
                by=by,
                output_path=output,
            )
        except ValueError as exc:
            console.print(f"[error]Error:[/error] {exc}")
            raise SystemExit(1) from exc

    console.print()

    # ── summary table ────────────────────────────────────────────────────────
    table = Table(
        title="MOFA+ Views Built",
        title_style="bold cyan",
        border_style="dim",
    )
    table.add_column("View", style="assay")
    table.add_column("Samples", style="score", justify="right")
    table.add_column("Features", style="score", justify="right")
    table.add_column("Missing (%)", style="dim", justify="right")

    for view_name, view_df in views.items():
        n_missing = int(view_df.isna().sum().sum())
        total = view_df.shape[0] * view_df.shape[1]
        pct = (n_missing / total * 100) if total > 0 else 0.0
        table.add_row(
            view_name,
            str(view_df.shape[0]),
            str(view_df.shape[1]),
            f"{pct:.1f}",
        )

    console.print(table)
    console.print()

    console.print(
        f"  [success]✓[/success] Built [bold]{len(views)}[/bold] view(s) "
        f"with [bold]{sum(v.shape[1] for v in views.values())}[/bold] "
        f"total features across "
        f"[bold]{next(iter(views.values())).shape[0]}[/bold] sample(s)."
    )
    console.print()

    # ── optional AnnData output ──────────────────────────────────────────────
    if output:
        console.print(f"  [success]✓[/success] Views saved to [bold]{output}[/bold]")

    # ── optional CSV output ──────────────────────────────────────────────────
    if output_csv_dir:
        out_dir = Path(output_csv_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        for view_name, view_df in views.items():
            csv_path = out_dir / f"{view_name}.tsv"
            view_df.to_csv(csv_path, sep="\t", na_rep="NaN")
            console.print(f"  [success]✓[/success] {view_name} → [bold]{csv_path}[/bold]")

    if not output and not output_csv_dir:
        console.print(
            "  [dim]Tip: Use -o views.h5ad or --output-csv-dir views/ to save.[/dim]"
        )


@mofa.command("run")
@click.argument("input_h5ad", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-k", "--factors",
    type=int,
    default=10,
    show_default=True,
    help="Number of latent factors.",
)
@click.option(
    "-o", "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save the trained model as .hdf5.",
)
@click.option(
    "--seed",
    type=int,
    default=42,
    show_default=True,
    help="Random seed.",
)
@click.option(
    "--convergence",
    type=click.Choice(["fast", "medium", "slow"]),
    default="fast",
    show_default=True,
    help="MOFA+ convergence mode.",
)
def mofa_run(
    input_h5ad: str,
    factors: int,
    output: str | None,
    seed: int,
    convergence: str,
) -> None:
    """Train a MOFA+ model on pre-built views.

    INPUT_H5AD is an .h5ad file produced by 'regvar mofa build -o'.

    Requires: pip install mofapy2 anndata

    \b
    Examples:
        regvar mofa run views.h5ad --factors 15 -o model.hdf5
        regvar mofa run views.h5ad --convergence medium
    """
    _banner()

    from .mofa_view import load_views_from_anndata

    # ── load views ───────────────────────────────────────────────────────────
    with console.status("Loading views…", spinner="dots"):
        try:
            views = load_views_from_anndata(input_h5ad)
        except ImportError as exc:
            console.print(f"[error]Error:[/error] {exc}")
            raise SystemExit(1) from exc
        except Exception as exc:
            console.print(f"[error]Error loading {input_h5ad}:[/error] {exc}")
            raise SystemExit(1) from exc

    console.print(
        f"  Loaded [bold]{len(views)}[/bold] view(s) from "
        f"[bold]{input_h5ad}[/bold]"
    )
    for vn, vdf in views.items():
        console.print(f"    [assay]{vn}[/assay]: {vdf.shape[0]} × {vdf.shape[1]}")

    # ── train ────────────────────────────────────────────────────────────────
    from .mofa_integration import train_mofa

    with console.status(
        f"Training MOFA+ with {factors} factors ({convergence} mode)…",
        spinner="dots",
    ):
        try:
            model = train_mofa(
                views,
                n_factors=factors,
                convergence_mode=convergence,
                seed=seed,
                output_path=output,
            )
        except ImportError as exc:
            console.print(f"[error]Error:[/error] {exc}")
            console.print(
                "\n  [dim]Tip: Install with [bold]pip install mofapy2 anndata[/bold] "
                "or use [bold]regvar mofa build[/bold] without training.[/dim]"
            )
            raise SystemExit(1) from exc
        except ValueError as exc:
            console.print(f"[error]Error:[/error] {exc}")
            raise SystemExit(1) from exc

    console.print()
    console.print(f"  [success]✓[/success] Training complete ({factors} factors).")

    if output:
        console.print(f"  [success]✓[/success] Model saved to [bold]{output}[/bold]")

    # ── variance explained summary ───────────────────────────────────────────
    try:
        r2 = model.get_variance_explained()
        if r2:
            table = Table(
                title="Variance Explained (R²) by View",
                title_style="bold cyan",
                border_style="dim",
            )
            table.add_column("View", style="assay")
            table.add_column("Total R²", style="score")
            for vn, per_factor in r2.items():
                table.add_row(vn, f"{sum(per_factor):.4f}")
            console.print()
            console.print(table)
    except Exception:
        pass   # variance explained not always available


@mofa.command("compare")
@click.argument("input_h5ad", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-k", "--factors",
    type=int,
    default=10,
    show_default=True,
    help="Number of latent factors for both models.",
)
@click.option(
    "-o", "--output",
    type=click.Path(dir_okay=False),
    default=None,
    help="Save the comparison report as markdown.",
)
@click.option(
    "--seed",
    type=int,
    default=42,
    show_default=True,
    help="Random seed.",
)
@click.option(
    "--variant-views",
    default=None,
    help="Comma-separated view names to treat as variant scores "
         "(default: auto-detect views with variant-like feature names).",
)
def mofa_compare(
    input_h5ad: str,
    factors: int,
    output: str | None,
    seed: int,
    variant_views: str | None,
) -> None:
    """Compare MOFA+ models with vs. without variant-effect views.

    Trains two models and quantifies how much additional structure the
    regulatory variant scores contribute to the factor decomposition.

    INPUT_H5AD is an .h5ad file produced by 'regvar mofa build -o'.

    Requires: pip install mofapy2 anndata

    \b
    Examples:
        regvar mofa compare views.h5ad -o comparison.md
        regvar mofa compare views.h5ad --factors 15
        regvar mofa compare views.h5ad --variant-views ATAC-seq,RNA-seq
    """
    _banner()

    from .mofa_view import load_views_from_anndata

    # ── load views ───────────────────────────────────────────────────────────
    with console.status("Loading views…", spinner="dots"):
        try:
            views = load_views_from_anndata(input_h5ad)
        except ImportError as exc:
            console.print(f"[error]Error:[/error] {exc}")
            raise SystemExit(1) from exc
        except Exception as exc:
            console.print(f"[error]Error loading {input_h5ad}:[/error] {exc}")
            raise SystemExit(1) from exc

    console.print(
        f"  Loaded [bold]{len(views)}[/bold] view(s) from [bold]{input_h5ad}[/bold]"
    )

    # ── determine variant vs. base views ─────────────────────────────────────
    if variant_views:
        vv_set = set(v.strip() for v in variant_views.split(","))
        unknown = vv_set - set(views.keys())
        if unknown:
            console.print(
                f"[error]Error:[/error] Unknown view(s): {', '.join(sorted(unknown))}"
            )
            raise SystemExit(1)
        base_views = {k: v for k, v in views.items() if k not in vv_set}
        variant_views_map = {k: v for k, v in views.items() if k in vv_set}
    else:
        # Auto-detect: views whose feature names contain variant-like patterns
        # (chr:pos:ref>alt_...)
        import re
        _variant_pat = re.compile(r"chr[\w]+:\d+:[ACGT]>[ACGT]")
        variant_views_map = {}
        base_views = {}
        for vn, vdf in views.items():
            n_variant_feats = sum(
                1 for c in vdf.columns if _variant_pat.match(str(c))
            )
            if n_variant_feats > 0 and n_variant_feats >= len(vdf.columns) * 0.5:
                variant_views_map[vn] = vdf
            else:
                base_views[vn] = vdf

    if not variant_views_map:
        console.print(
            "[warning]⚠[/warning] No variant-like views detected. "
            "Use --variant-views to specify them explicitly."
        )
        console.print(
            "  [dim]Tip: Variant views have feature names like "
            "'chr8:127401060:G>T_combined_ATAC'.[/dim]"
        )
        raise SystemExit(1)

    if not base_views:
        console.print(
            "[warning]⚠[/warning] All views appear to be variant views. "
            "Nothing to compare against. Add base (non-variant) views to the file."
        )
        raise SystemExit(1)

    console.print(
        f"  Base views: [assay]{', '.join(base_views)}[/assay]"
    )
    console.print(
        f"  Variant views: [assay]{', '.join(variant_views_map)}[/assay]"
    )
    console.print()

    # ── run comparison ───────────────────────────────────────────────────────
    from .mofa_integration import compare_variant_contribution

    with console.status(
        f"Training models with {factors} factors (this may take a while)…",
        spinner="dots",
    ):
        try:
            result = compare_variant_contribution(
                base_views=base_views,
                variant_views=variant_views_map,
                n_factors=factors,
                seed=seed,
                output_path=output,
            )
        except ImportError as exc:
            console.print(f"[error]Error:[/error] {exc}")
            console.print(
                "\n  [dim]Tip: Install with [bold]pip install mofapy2 anndata[/bold][/dim]"
            )
            raise SystemExit(1) from exc
        except ValueError as exc:
            console.print(f"[error]Error:[/error] {exc}")
            raise SystemExit(1) from exc

    console.print()
    console.print(f"  [success]✓[/success] Comparison complete.")
    stdout_console.print(Markdown(result["report"]))

    if output:
        console.print(f"  [success]✓[/success] Report saved to [bold]{output}[/bold]")


# ---------------------------------------------------------------------------
# regvar cache  (command group)
# ---------------------------------------------------------------------------

@cli.group()
def cache() -> None:
    """Manage the AlphaGenome score cache.

    \b
    Subcommands:
        info     Show cache directory, size, and entry count
        clear    Delete cached entries (optionally by age)
    """


@cache.command("info")
def cache_info() -> None:
    """Show cache directory, total size, and entry count."""
    from .alphagenome_client import ClientConfig

    cfg = ClientConfig()
    cache_dir = cfg.cache_dir
    entries = list(cache_dir.glob("*.pkl")) + list(cache_dir.glob("*.parquet"))
    total_bytes = sum(f.stat().st_size for f in entries)

    _banner()
    console.print(f"  Cache directory: [bold]{cache_dir}[/bold]")
    console.print(f"  Entries:         [bold]{len(entries)}[/bold]")
    if total_bytes < 1024:
        size_str = f"{total_bytes} B"
    elif total_bytes < 1024 * 1024:
        size_str = f"{total_bytes / 1024:.1f} KB"
    else:
        size_str = f"{total_bytes / (1024 * 1024):.1f} MB"
    console.print(f"  Total size:      [bold]{size_str}[/bold]")


@cache.command("clear")
@click.option(
    "--older-than",
    type=int,
    default=None,
    help="Only delete entries older than this many days.",
)
@click.option(
    "--yes", "-y",
    is_flag=True,
    default=False,
    help="Skip the confirmation prompt.",
)
def cache_clear(older_than: int | None, yes: bool) -> None:
    """Delete cached score entries.

    Without --older-than, deletes ALL cache files. With --older-than N,
    only deletes entries not modified in the last N days.
    """
    from .alphagenome_client import ClientConfig

    cfg = ClientConfig()
    cache_dir = cfg.cache_dir
    entries = list(cache_dir.glob("*.pkl")) + list(cache_dir.glob("*.parquet"))

    if older_than is not None:
        import time as _time
        cutoff = _time.time() - older_than * 86400
        entries = [f for f in entries if f.stat().st_mtime < cutoff]

    if not entries:
        _banner()
        console.print("  [dim]No matching cache entries found.[/dim]")
        return

    if not yes:
        msg = f"Delete {len(entries)} cache entr{'y' if len(entries) == 1 else 'ies'}"
        if older_than is not None:
            msg += f" older than {older_than} day(s)"
        msg += "?"
        if not click.confirm(msg, default=False):
            console.print("  [dim]Aborted.[/dim]")
            return

    for f in entries:
        f.unlink()

    _banner()
    console.print(
        f"  [success]✓[/success] Deleted [bold]{len(entries)}[/bold] "
        f"cache entr{'y' if len(entries) == 1 else 'ies'}."
    )


# ---------------------------------------------------------------------------
# regvar assays  (flat command, unchanged)
# ---------------------------------------------------------------------------

@cli.command()
def assays() -> None:
    """List the functional genomic assays AlphaGenome can score."""
    _banner()
    from .alphagenome_client import ASSAY_TO_OUTPUT
    _print_assay_table(ASSAY_TO_OUTPUT)


# ---------------------------------------------------------------------------
# regvar tui  (flat command, unchanged)
# ---------------------------------------------------------------------------

@cli.command()
def tui() -> None:
    """Launch the interactive TUI (requires: pip install textual)."""
    try:
        from .tui import RegvarApp
    except ImportError:
        console.print(
            "[error]Error:[/error] textual is not installed. "
            "Run: [bold]pip install textual>=0.55.0[/bold]"
        )
        raise SystemExit(1)
    RegvarApp().run()


# ---------------------------------------------------------------------------
# regvar run  (backward-compat alias → agent run)
# ---------------------------------------------------------------------------

@cli.command(
    hidden=True,
    context_settings={"ignore_unknown_options": True, "allow_extra_args": True},
)
@click.argument("candidates_tsv", type=click.Path(exists=True, dir_okay=False))
@click.option("-a", "--assays",      default=None)
@click.option("-t", "--tissue",      default=None)
@click.option("-n", "--top-n",       type=int, default=10,            show_default=True)
@click.option("-m", "--model",       default="deepseek-v4-pro",       show_default=True)
@click.option("-o", "--output",      type=click.Path(dir_okay=False), default=None)
@click.option("--output-tsv",        type=click.Path(dir_okay=False), default=None)
@click.option("--output-json",       type=click.Path(dir_okay=False), default=None)
@click.option("--dry-run",           is_flag=True, default=False)
@click.option("-v", "--verbose",     is_flag=True, default=False)
@click.option("-q", "--quiet",       is_flag=True, default=False)
@click.option("--force",             is_flag=True, default=False)
@click.option("--max-concurrent",    type=int, default=None)
def run(
    candidates_tsv: str,
    assays: str | None,
    tissue: str | None,
    top_n: int,
    model: str,
    output: str | None,
    output_tsv: str | None,
    output_json: str | None,
    dry_run: bool,
    verbose: bool,
    quiet: bool,
    force: bool,
    max_concurrent: int | None,
) -> None:
    """[Alias] Run the agent triage loop — same as 'regvar agent run'.

    This command is kept for backward compatibility with existing scripts.
    Prefer 'regvar agent run' for new workflows.
    """
    console.print(
        "  [dim]Note: 'regvar run' is an alias for 'regvar agent run'.[/dim]\n"
    )
    _run_agent_core(
        candidates_tsv=candidates_tsv,
        assays=assays,
        tissue=tissue,
        top_n=top_n,
        model=model,
        output=output,
        output_tsv=output_tsv,
        output_json=output_json,
        dry_run=dry_run,
        verbose=verbose,
        quiet=quiet,
        force=force,
        max_concurrent=max_concurrent,
    )


# ---------------------------------------------------------------------------
# report command group
# ---------------------------------------------------------------------------

@cli.group()
def report() -> None:
    """Generate one-page variant report PDFs from saved agent output.

    \b
    Subcommands:
        generate    Assemble per-candidate PDFs from scores TSV + markdown report
    """


@report.command("generate")
@click.argument("scores_tsv", type=click.Path(exists=True, dir_okay=False))
@click.argument("report_md", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-c", "--candidates", "candidates_tsv",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to original candidate variants TSV (for metadata).",
)
@click.option(
    "-o", "--output-dir",
    type=click.Path(file_okay=False),
    default="./reports/",
    show_default=True,
    help="Output directory for per-candidate PDFs.",
)
@click.option(
    "-n", "--top-n",
    type=int,
    default=5,
    show_default=True,
    help="Generate reports for top N variants by max effect.",
)
@click.option(
    "--all", "all_variants",
    is_flag=True,
    help="Generate reports for ALL candidates (overrides --top-n).",
)
@click.option(
    "--template",
    type=str,
    default="default",
    show_default=True,
    help="LaTeX template name (from regvar/templates/).",
)
@click.option(
    "--tissue",
    type=str,
    default=None,
    help="Tissue/cell-type context for report header.",
)
@click.option(
    "--prefix",
    type=str,
    default="variant_report",
    show_default=True,
    help="Output filename prefix.",
)
@click.option(
    "--no-compile",
    is_flag=True,
    help="Write .tex files only; skip pdflatex compilation.",
)
@click.option(
    "--clean",
    is_flag=True,
    help="Remove .aux, .log files after successful compilation.",
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Print pdflatex output to stderr.",
)
def report_generate(
    scores_tsv: str,
    report_md: str,
    candidates_tsv: str | None,
    output_dir: str,
    top_n: int,
    all_variants: bool,
    template: str,
    tissue: str | None,
    prefix: str,
    no_compile: bool,
    clean: bool,
    verbose: bool,
) -> None:
    """Generate one-page PDF reports for top regulatory variant candidates.

    Reads the scored variants TSV and agent analysis markdown produced by
    ``regvar agent run``, then assembles a one-page LaTeX PDF per candidate
    containing:

    \b
    - REF/ALT locus context figure
    - Regulatory effect magnitude barplot
    - Scoring summary table
    - Mechanistic interpretation + validation recommendation

    \b
    Examples:

        regvar report generate scores.tsv analysis.md

        regvar report generate scores.tsv analysis.md --top-n 10 --tissue "Prostate"

        regvar report generate scores.tsv analysis.md --all --no-compile --clean
    """
    console = Console(theme=REGVAR_THEME)

    # ── validate inputs ────────────────────────────────────────────────
    scores_path = Path(scores_tsv)
    report_path = Path(report_md)

    # ── resolve candidates TSV ──────────────────────────────────────────
    if candidates_tsv is None:
        default_candidates = Path("examples/candidate_variants.tsv")
        if default_candidates.exists():
            candidates_path = default_candidates
        else:
            candidates_path = Path("/nonexistent/path.tsv")
    else:
        candidates_path = Path(candidates_tsv)

    # ── check for jinja2 ────────────────────────────────────────────────
    try:
        from regvar.report import ReportBuilder
    except ImportError:
        console.print(
            "[warning]jinja2 required for report generation.[/]\n"
            "Install with: pip install regvar[report]"
        )
        raise click.Abort()

    # ── check for pdflatex ──────────────────────────────────────────────
    if not no_compile:
        import shutil
        if shutil.which("pdflatex") is None:
            console.print(
                "[warning]pdflatex not found on PATH.[/]\n"
                "Use --no-compile to generate .tex files only, "
                "or install a LaTeX distribution (texlive)."
            )

    # ── build reports ──────────────────────────────────────────────────
    console.print(f"[info]Reading scores:[/] {scores_path.name}")
    console.print(f"[info]Reading report:[/] {report_path.name}")

    builder = ReportBuilder(template_name=template, tissue=tissue)
    reports = builder.build(
        candidates_tsv=candidates_path,
        scores_tsv=scores_path,
        report_md=report_path,
    )

    if not reports:
        console.print("[warning]No reports generated. Check input files.[/]")
        return

    n_candidates = len(reports)
    effective_top_n = None if all_variants else min(top_n, n_candidates)
    n_to_generate = n_candidates if all_variants else min(top_n, n_candidates)

    console.print(
        f"[info]Candidates found:[/] {n_candidates}  "
        f"[info]Generating:[/] {n_to_generate}"
    )

    # ── generate PDFs ──────────────────────────────────────────────────
    output_path = Path(output_dir)
    compile_flag = not no_compile

    paths = builder.generate_all(
        reports=reports,
        output_dir=output_path,
        top_n=effective_top_n,
        prefix=prefix,
        tissue=tissue,
        compile_pdf=compile_flag,
        clean=clean,
        verbose=verbose,
    )

    # ── summary ─────────────────────────────────────────────────────────
    console.print(f"\n[success]✓ Generated {len(paths)} reports[/]")
    console.print(f"  Output directory: {output_path.resolve()}")
    for p in paths[:5]:
        console.print(f"  • {p.name}")
    if len(paths) > 5:
        console.print(f"  ... and {len(paths) - 5} more")

    if not compile_flag:
        console.print(
            "\n[info].tex files written. Compile with:[/]\n"
            f"  cd {output_path.resolve()} && for f in *.tex; do pdflatex $f; done"
        )
