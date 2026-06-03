"""Interactive TUI for regvar-agent, built with Textual.

Launch with:

    python -m regvar tui

The TUI mirrors the three CLI commands as interactive tabs:

    Run Agent    — equivalent to ``regvar run <file.tsv>``
    Score Variant — equivalent to ``regvar score <chr> <pos> <ref> <alt>``
    Assays       — equivalent to ``regvar assays``

All existing CLI commands are completely unaffected.
"""

from __future__ import annotations

import asyncio
import json
import time
from pathlib import Path
from typing import Any

from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal, ScrollableContainer, Vertical
from textual.css.query import NoMatches
from textual.reactive import reactive
from textual.widgets import (
    Button,
    Checkbox,
    DataTable,
    Footer,
    Header,
    Input,
    Label,
    Markdown,
    RichLog,
    Static,
    TabbedContent,
    TabPane,
)


# ---------------------------------------------------------------------------
# CSS
# ---------------------------------------------------------------------------

APP_CSS = """
/* ── Global ──────────────────────────────────────────────────────────────── */
Screen {
    background: #0d1117;
}

Header {
    background: #161b22;
    color: #58a6ff;
    text-style: bold;
}

Footer {
    background: #161b22;
    color: #8b949e;
}

/* ── Tabs ─────────────────────────────────────────────────────────────────── */
TabbedContent {
    background: #0d1117;
}

TabPane {
    padding: 1 2;
    background: #0d1117;
}

/* ── Section titles ──────────────────────────────────────────────────────── */
.section-title {
    color: #58a6ff;
    text-style: bold;
    padding: 0 0 1 0;
}

/* ── Form rows ───────────────────────────────────────────────────────────── */
.form-row {
    height: auto;
    margin-bottom: 1;
    align: left middle;
}

.form-label {
    width: 18;
    color: #8b949e;
    text-align: right;
    padding-right: 1;
}

.form-input {
    width: 1fr;
}

.form-input-short {
    width: 12;
}

.form-row-checkboxes {
    height: auto;
    margin-bottom: 1;
    padding-left: 19;
}

/* ── Buttons ─────────────────────────────────────────────────────────────── */
.btn-run {
    margin-top: 1;
    margin-left: 19;
    background: #238636;
    color: #ffffff;
    border: none;
    min-width: 20;
}

.btn-run:hover {
    background: #2ea043;
}

.btn-run.-disabled {
    background: #21262d;
    color: #484f58;
}

.btn-clear {
    margin-top: 1;
    margin-left: 2;
    background: #21262d;
    color: #8b949e;
    border: none;
    min-width: 10;
}

.btn-clear:hover {
    background: #30363d;
    color: #e6edf3;
}

/* ── Progress log ────────────────────────────────────────────────────────── */
.log-container {
    height: 10;
    border: solid #30363d;
    background: #010409;
    margin-top: 1;
}

#run-log, #score-log {
    height: 10;
    background: #010409;
}

/* ── Results / output ────────────────────────────────────────────────────── */
.results-container {
    border: solid #30363d;
    background: #0d1117;
    margin-top: 1;
    height: 1fr;
    overflow-y: auto;
}

#run-results {
    padding: 1;
}

/* ── Assays tab ──────────────────────────────────────────────────────────── */
#assays-table {
    height: 1fr;
}

/* ── Status bar ──────────────────────────────────────────────────────────── */
.status-bar {
    height: 1;
    color: #8b949e;
    padding: 0 1;
    background: #161b22;
}

.status-bar.--success {
    color: #3fb950;
}

.status-bar.--error {
    color: #f85149;
}

.status-bar.--running {
    color: #d29922;
}

/* ── Score results table ─────────────────────────────────────────────────── */
#score-table-container {
    height: 1fr;
    border: solid #30363d;
    background: #0d1117;
    margin-top: 1;
}
"""


# ---------------------------------------------------------------------------
# Run Agent tab
# ---------------------------------------------------------------------------

class RunAgentTab(TabPane):
    """Interactive form to run the full agent loop."""

    _running: reactive[bool] = reactive(False)

    def compose(self) -> ComposeResult:
        yield Static("Run Agent", classes="section-title")

        # TSV path
        with Horizontal(classes="form-row"):
            yield Label("Candidates TSV:", classes="form-label")
            yield Input(
                placeholder="examples/candidate_variants.tsv",
                id="run-tsv",
                classes="form-input",
            )

        # Assays
        with Horizontal(classes="form-row"):
            yield Label("Assays:", classes="form-label")
            yield Input(
                placeholder="ATAC-seq, RNA-seq, ChIP-seq histone, Hi-C / pcHi-C  (leave blank for defaults)",
                id="run-assays",
                classes="form-input",
            )

        # Tissue
        with Horizontal(classes="form-row"):
            yield Label("Tissue terms:", classes="form-label")
            yield Input(
                placeholder="UBERON:0002367  (leave blank for defaults)",
                id="run-tissue",
                classes="form-input",
            )

        # Top N + Model on one row
        with Horizontal(classes="form-row"):
            yield Label("Top N:", classes="form-label")
            yield Input(value="10", id="run-top-n", classes="form-input-short")
            yield Label("  Model:", classes="form-label")
            yield Input(
                value="deepseek-v4-pro",
                id="run-model",
                classes="form-input",
            )

        # Output paths
        with Horizontal(classes="form-row"):
            yield Label("Output (md):", classes="form-label")
            yield Input(
                placeholder="results/report.md  (optional)",
                id="run-output",
                classes="form-input",
            )

        with Horizontal(classes="form-row"):
            yield Label("Output (tsv):", classes="form-label")
            yield Input(
                placeholder="results/scores.tsv  (optional)",
                id="run-output-tsv",
                classes="form-input",
            )

        # Checkboxes
        with Horizontal(classes="form-row-checkboxes"):
            yield Checkbox("Dry-run", id="run-dry-run")
            yield Checkbox("Verbose", id="run-verbose")

        # Buttons
        with Horizontal(classes="form-row"):
            yield Button("▶  Run Agent", id="btn-run-agent", classes="btn-run")
            yield Button("Clear", id="btn-run-clear", classes="btn-clear")

        # Status
        yield Static("", id="run-status", classes="status-bar")

        # Progress log
        with Container(classes="log-container"):
            yield RichLog(id="run-log", highlight=True, markup=True, wrap=True)

        # Results (markdown)
        with ScrollableContainer(classes="results-container"):
            yield Markdown("", id="run-results")

    # -- event handlers -------------------------------------------------------

    @on(Button.Pressed, "#btn-run-agent")
    async def handle_run(self) -> None:
        if self._running:
            return
        tsv = self.query_one("#run-tsv", Input).value.strip()
        if not tsv:
            self._set_status("Please enter a path to a candidates TSV file.", kind="error")
            return

        tsv_path = Path(tsv)
        if not tsv_path.exists():
            self._set_status(f"File not found: {tsv}", kind="error")
            return

        self._running = True
        self.query_one("#btn-run-agent", Button).disabled = True
        log = self.query_one("#run-log", RichLog)
        log.clear()
        await self.query_one("#run-results", Markdown).update("")
        self._set_status("Running agent…", kind="running")

        # collect params
        assays_raw = self.query_one("#run-assays", Input).value.strip()
        tissue_raw = self.query_one("#run-tissue", Input).value.strip()
        top_n_raw = self.query_one("#run-top-n", Input).value.strip()
        model = self.query_one("#run-model", Input).value.strip() or "deepseek-v4-pro"
        output = self.query_one("#run-output", Input).value.strip() or None
        output_tsv = self.query_one("#run-output-tsv", Input).value.strip() or None
        dry_run = self.query_one("#run-dry-run", Checkbox).value
        verbose = self.query_one("#run-verbose", Checkbox).value

        assay_list = [a.strip() for a in assays_raw.split(",")] if assays_raw else None
        tissue_list = [t.strip() for t in tissue_raw.split(",")] if tissue_raw else None
        try:
            top_n = int(top_n_raw) if top_n_raw else 10
        except ValueError:
            top_n = 10

        self._do_run(tsv_path, assay_list, tissue_list, top_n, model,
                     output, output_tsv, dry_run, verbose)

    @work(thread=True)
    def _do_run(
        self,
        tsv_path: Path,
        assay_list: list[str] | None,
        tissue_list: list[str] | None,
        top_n: int,
        model: str,
        output: str | None,
        output_tsv: str | None,
        dry_run: bool,
        verbose: bool,
    ) -> None:
        """Run agent in a background thread so TUI stays responsive."""
        from .variants import read_candidates_tsv
        from .agent import run_agent

        suffix = "".join(tsv_path.suffixes).lower()
        if suffix in (".vcf", ".vcf.gz"):
            from .variants import read_candidates_vcf
            _reader = read_candidates_vcf
        else:
            _reader = read_candidates_tsv

        def _post(msg: str) -> None:
            self.app.call_from_thread(self._append_log, msg)

        try:
            candidates = _reader(tsv_path)
            _post(f"[cyan]Loaded {len(candidates)} variant(s) from {tsv_path}[/cyan]")

            if dry_run:
                _post("[yellow]Dry-run mode — no API calls will be made.[/yellow]")
                if assay_list:
                    _post(f"  Assays: {', '.join(assay_list)}")
                if tissue_list:
                    _post(f"  Tissue terms: {', '.join(tissue_list)}")
                _post(f"  Top N: {top_n}  |  Model: {model}")
                _post("[green]✓ Input validated. Ready to run.[/green]")
                self.app.call_from_thread(self._finish_run, "Dry-run complete.", None, None, None)
                return

            import os
            if not os.environ.get("DEEPSEEK_API_KEY"):
                self.app.call_from_thread(
                    self._finish_run,
                    "DEEPSEEK_API_KEY not set.",
                    None, None, "error",
                )
                return

            tool_call_log: list[dict] = []

            def _on_tool_call(name: str, args: dict, result: str) -> None:
                tool_call_log.append({"tool": name, "args": args, "result": result})
                if name == "score_regulatory_variant":
                    v = f"{args.get('chromosome', '?')}:{args.get('position', '?')}"
                    msg = f"  [dim]→[/dim] scoring [magenta]{v}[/magenta]"
                else:
                    msg = f"  [dim]→[/dim] [blue]{name}[/blue]"
                if verbose:
                    msg += f"  [dim]{json.dumps(args, default=str)}[/dim]"
                _post(msg)

            start = time.time()
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
            _post(f"[green]✓ Agent finished in {elapsed:.1f}s ({len(tool_call_log)} tool call(s))[/green]")

            # Save files if requested
            if output:
                p = Path(output)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(answer, encoding="utf-8")
                _post(f"[green]✓ Report saved → {output}[/green]")

            if output_tsv:
                score_records: list[dict] = []
                for entry in tool_call_log:
                    if entry["tool"] == "score_regulatory_variant":
                        try:
                            parsed = (
                                json.loads(entry["result"])
                                if isinstance(entry["result"], str)
                                else entry["result"]
                            )
                            for rec in parsed.get("top_effects", []):
                                rec["query_variant"] = parsed.get("variant", "")
                                score_records.append(rec)
                        except (json.JSONDecodeError, TypeError):
                            pass
                if score_records:
                    from .variants import write_scores_tsv
                    write_scores_tsv(score_records, output_tsv)
                    _post(f"[green]✓ Scores saved → {output_tsv}[/green]")

            self.app.call_from_thread(self._finish_run, f"Done in {elapsed:.1f}s", answer, None, None)

        except Exception as exc:  # noqa: BLE001
            _post(f"[red]Error: {exc}[/red]")
            self.app.call_from_thread(self._finish_run, f"Error: {exc}", None, None, "error")

    def _append_log(self, msg: str) -> None:
        self.query_one("#run-log", RichLog).write(msg)

    def _finish_run(
        self,
        status: str,
        answer: str | None,
        _unused: Any,
        kind: str | None,
    ) -> None:
        self._running = False
        self.query_one("#btn-run-agent", Button).disabled = False
        self._set_status(status, kind=kind or ("success" if not kind else kind))
        if answer:
            self.app.call_later(self._update_markdown, answer)

    async def _update_markdown(self, text: str) -> None:
        await self.query_one("#run-results", Markdown).update(text)

    @on(Button.Pressed, "#btn-run-clear")
    async def handle_clear(self) -> None:
        self.query_one("#run-log", RichLog).clear()
        await self.query_one("#run-results", Markdown).update("")
        self._set_status("")

    def _set_status(self, msg: str, kind: str | None = None) -> None:
        status = self.query_one("#run-status", Static)
        status.update(msg)
        status.remove_class("--success", "--error", "--running")
        if kind == "error":
            status.add_class("--error")
        elif kind == "running":
            status.add_class("--running")
        elif kind == "success":
            status.add_class("--success")


# ---------------------------------------------------------------------------
# Score Variant tab
# ---------------------------------------------------------------------------

class ScoreVariantTab(TabPane):
    """Interactive form to score a single variant."""

    _running: reactive[bool] = reactive(False)

    def compose(self) -> ComposeResult:
        yield Static("Score a Single Variant", classes="section-title")

        with Horizontal(classes="form-row"):
            yield Label("Chromosome:", classes="form-label")
            yield Input(placeholder="chr8", id="score-chrom", classes="form-input-short")
            yield Label("  Position:", classes="form-label")
            yield Input(placeholder="127401060", id="score-pos", classes="form-input-short")
            yield Label("  Ref:", classes="form-label")
            yield Input(placeholder="G", id="score-ref", classes="form-input-short")
            yield Label("  Alt:", classes="form-label")
            yield Input(placeholder="T", id="score-alt", classes="form-input-short")

        with Horizontal(classes="form-row"):
            yield Label("Assays:", classes="form-label")
            yield Input(
                placeholder="ATAC-seq, RNA-seq, ChIP-seq histone, Hi-C / pcHi-C  (leave blank for defaults)",
                id="score-assays",
                classes="form-input",
            )

        with Horizontal(classes="form-row"):
            yield Label("Tissue terms:", classes="form-label")
            yield Input(
                placeholder="UBERON:0002367  (leave blank for defaults)",
                id="score-tissue",
                classes="form-input",
            )

        with Horizontal(classes="form-row"):
            yield Label("Top N:", classes="form-label")
            yield Input(value="10", id="score-top-n", classes="form-input-short")

        with Horizontal(classes="form-row"):
            yield Label("Output TSV:", classes="form-label")
            yield Input(
                placeholder="results/scores.tsv  (optional)",
                id="score-output-tsv",
                classes="form-input",
            )

        with Horizontal(classes="form-row"):
            yield Button("▶  Score", id="btn-score", classes="btn-run")
            yield Button("Clear", id="btn-score-clear", classes="btn-clear")

        yield Static("", id="score-status", classes="status-bar")

        with Container(classes="log-container"):
            yield RichLog(id="score-log", highlight=True, markup=True, wrap=True)

        with Container(id="score-table-container"):
            yield DataTable(id="score-table", zebra_stripes=True, cursor_type="row")

    @on(Button.Pressed, "#btn-score")
    async def handle_score(self) -> None:
        if self._running:
            return

        chrom = self.query_one("#score-chrom", Input).value.strip()
        pos_raw = self.query_one("#score-pos", Input).value.strip()
        ref = self.query_one("#score-ref", Input).value.strip().upper()
        alt = self.query_one("#score-alt", Input).value.strip().upper()

        if not all([chrom, pos_raw, ref, alt]):
            self._set_status("Chromosome, position, ref, and alt are all required.", kind="error")
            return
        try:
            pos = int(pos_raw)
        except ValueError:
            self._set_status("Position must be an integer.", kind="error")
            return

        try:
            from .variants import validate_chromosome
            chrom = validate_chromosome(chrom)
        except ValueError as exc:
            self._set_status(str(exc), kind="error")
            return

        assays_raw = self.query_one("#score-assays", Input).value.strip()
        tissue_raw = self.query_one("#score-tissue", Input).value.strip()
        top_n_raw = self.query_one("#score-top-n", Input).value.strip()
        output_tsv = self.query_one("#score-output-tsv", Input).value.strip() or None

        assay_list = [a.strip() for a in assays_raw.split(",")] if assays_raw else None
        tissue_list = [t.strip() for t in tissue_raw.split(",")] if tissue_raw else None
        try:
            top_n = int(top_n_raw) if top_n_raw else 10
        except ValueError:
            top_n = 10

        self._running = True
        self.query_one("#btn-score", Button).disabled = True
        self.query_one("#score-log", RichLog).clear()
        table = self.query_one("#score-table", DataTable)
        table.clear(columns=True)
        self._set_status("Querying AlphaGenome…", kind="running")

        self._do_score(chrom, pos, ref, alt, assay_list, tissue_list, top_n, output_tsv)

    @work(thread=True)
    def _do_score(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        assay_list: list[str] | None,
        tissue_list: list[str] | None,
        top_n: int,
        output_tsv: str | None,
    ) -> None:
        from .tools import tool_score_regulatory_variant

        def _post(msg: str) -> None:
            self.app.call_from_thread(self._append_log, msg)

        try:
            kwargs: dict[str, Any] = {
                "chromosome": chrom, "position": pos,
                "ref": ref, "alt": alt, "top_n": top_n,
            }
            if assay_list:
                kwargs["assays"] = assay_list
            if tissue_list:
                kwargs["ontology_terms"] = tissue_list

            _post(f"[cyan]Scoring {chrom}:{pos}:{ref}>{alt}…[/cyan]")
            result = tool_score_regulatory_variant(**kwargs)
            _post(
                f"[green]✓ Done — {result['n_total_scores']} scores across "
                f"{len(result['assays_scored'])} assay(s)[/green]"
            )

            if output_tsv and result["top_effects"]:
                from .variants import write_scores_tsv
                write_scores_tsv(result["top_effects"], output_tsv)
                _post(f"[green]✓ Saved → {output_tsv}[/green]")

            self.app.call_from_thread(
                self._populate_table, result["top_effects"], result["variant"]
            )

        except Exception as exc:  # noqa: BLE001
            _post(f"[red]Error: {exc}[/red]")
            self.app.call_from_thread(self._finish_score, f"Error: {exc}", "error")

    def _append_log(self, msg: str) -> None:
        self.query_one("#score-log", RichLog).write(msg)

    def _populate_table(self, records: list[dict], variant_id: str) -> None:
        table = self.query_one("#score-table", DataTable)
        table.clear(columns=True)

        if not records:
            self._finish_score("No scores returned.", "error")
            return

        display_cols = [
            ("gene_name", "Gene"),
            ("assay", "Assay"),
            ("output_type", "Output Type"),
            ("biosample_name", "Biosample"),
            ("raw_score", "Raw Score"),
            ("quantile_score", "Quantile"),
            ("abs_effect", "|Effect|"),
        ]
        active = [(k, h) for k, h in display_cols if k in records[0]]

        table.add_column("#", width=4)
        for _, header in active:
            table.add_column(header)

        for i, rec in enumerate(records, 1):
            row = [str(i)]
            for key, _ in active:
                val = rec.get(key, "")
                if isinstance(val, float):
                    row.append(f"{val:+.4f}" if key != "abs_effect" else f"{val:.4f}")
                else:
                    row.append(str(val))
            table.add_row(*row)

        self._finish_score(f"Scored {variant_id} — {len(records)} top effect(s) shown.", None)

    def _finish_score(self, msg: str, kind: str | None) -> None:
        self._running = False
        self.query_one("#btn-score", Button).disabled = False
        self._set_status(msg, kind=kind or "success")

    @on(Button.Pressed, "#btn-score-clear")
    def handle_clear(self) -> None:
        self.query_one("#score-log", RichLog).clear()
        self.query_one("#score-table", DataTable).clear(columns=True)
        self._set_status("")

    def _set_status(self, msg: str, kind: str | None = None) -> None:
        status = self.query_one("#score-status", Static)
        status.update(msg)
        status.remove_class("--success", "--error", "--running")
        if kind == "error":
            status.add_class("--error")
        elif kind == "running":
            status.add_class("--running")
        elif kind == "success":
            status.add_class("--success")


# ---------------------------------------------------------------------------
# Assays tab
# ---------------------------------------------------------------------------

class AssaysTab(TabPane):
    """Display the supported assays table."""

    def compose(self) -> ComposeResult:
        yield Static("Supported AlphaGenome Assays", classes="section-title")
        yield DataTable(id="assays-table", zebra_stripes=True, cursor_type="row")

    def on_mount(self) -> None:
        from .alphagenome_client import ASSAY_TO_OUTPUT

        readout_map = {
            "ATAC": "Chromatin accessibility",
            "DNASE": "Chromatin accessibility (alt)",
            "RNA_SEQ": "Gene expression",
            "CHIP_HISTONE": "H3K4me3/me1, H3K27ac, H3K27me3",
            "CHIP_TF": "Transcription-factor binding",
            "CONTACT_MAPS": "Enhancer–promoter looping",
            "SPLICE_SITES": "Splice-site usage",
        }

        table = self.query_one("#assays-table", DataTable)
        table.add_column("#", width=4)
        table.add_column("Assay Label", width=22)
        table.add_column("AlphaGenome Output", width=20)
        table.add_column("Wet-Lab Readout")

        for i, (label, output) in enumerate(ASSAY_TO_OUTPUT.items(), 1):
            table.add_row(str(i), label, output, readout_map.get(output, ""))


# ---------------------------------------------------------------------------
# Main App
# ---------------------------------------------------------------------------

class RegvarApp(App):
    """Textual TUI for regvar — AlphaGenome regulatory variant triage."""

    TITLE = "regvar"
    SUB_TITLE = "AlphaGenome-powered regulatory variant triage  v0.1"
    CSS = APP_CSS

    BINDINGS = [
        Binding("ctrl+q", "quit", "Quit", show=True),
        Binding("f1", "focus_help", "Help", show=True),
        Binding("1", "switch_tab('run')", "Run Agent", show=False),
        Binding("2", "switch_tab('score')", "Score Variant", show=False),
        Binding("3", "switch_tab('assays')", "Assays", show=False),
    ]

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        with TabbedContent(initial="run"):
            yield RunAgentTab("▶  Run Agent", id="run")
            yield ScoreVariantTab("⬡  Score Variant", id="score")
            yield AssaysTab("☰  Assays", id="assays")
        yield Footer()

    def action_switch_tab(self, tab_id: str) -> None:
        try:
            self.query_one(TabbedContent).active = tab_id
        except NoMatches:
            pass

    def action_focus_help(self) -> None:
        self.notify(
            "1 → Run Agent tab  |  2 → Score Variant tab  |  3 → Assays tab  |  "
            "Ctrl+Q → Quit",
            title="Keyboard shortcuts",
            timeout=6,
        )


# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

def main() -> None:
    RegvarApp().run()


if __name__ == "__main__":
    main()
