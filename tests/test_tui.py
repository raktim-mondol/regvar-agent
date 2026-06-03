"""Smoke tests for the Textual TUI.

These tests use Textual's headless Pilot to verify that:
  - the app mounts without error
  - all three tabs are present and navigable
  - the Assays tab populates its table immediately on mount

No API calls are made — the Assays tab reads ASSAY_TO_OUTPUT directly,
and the Run / Score tabs only execute their workers when a button is pressed.
"""

from __future__ import annotations

import pytest

# Skip the whole module gracefully if textual is not installed.
textual = pytest.importorskip("textual", reason="textual not installed")


@pytest.mark.asyncio
async def test_app_mounts() -> None:
    """App mounts without raising any exceptions."""
    from regvar.tui import RegvarApp

    app = RegvarApp()
    async with app.run_test() as pilot:
        # Check the app title is set.
        assert app.TITLE == "regvar"


@pytest.mark.asyncio
async def test_tabs_present() -> None:
    """All three tabs are reachable by keyboard shortcut."""
    from textual.widgets import TabbedContent
    from regvar.tui import RegvarApp

    app = RegvarApp()
    async with app.run_test() as pilot:
        tc = app.query_one(TabbedContent)
        assert tc.active == "run"

        await pilot.press("2")
        assert tc.active == "score"

        await pilot.press("3")
        assert tc.active == "assays"

        await pilot.press("1")
        assert tc.active == "run"


@pytest.mark.asyncio
async def test_assays_table_populated() -> None:
    """Assays tab DataTable has rows equal to ASSAY_TO_OUTPUT length."""
    from textual.widgets import DataTable, TabbedContent
    from regvar.alphagenome_client import ASSAY_TO_OUTPUT
    from regvar.tui import RegvarApp

    app = RegvarApp()
    async with app.run_test() as pilot:
        await pilot.press("3")  # switch to Assays tab
        await pilot.pause()     # let on_mount finish
        table = app.query_one("#assays-table", DataTable)
        assert table.row_count == len(ASSAY_TO_OUTPUT)


@pytest.mark.asyncio
async def test_score_set_status_error_class() -> None:
    """_set_status(kind='error') applies the --error CSS class."""
    from textual.widgets import Static
    from regvar.tui import RegvarApp, ScoreVariantTab

    app = RegvarApp()
    async with app.run_test(size=(220, 50)) as pilot:
        await pilot.pause()
        tab = app.query_one(ScoreVariantTab)
        tab._set_status("something went wrong", kind="error")
        await pilot.pause()
        status = tab.query_one("#score-status", Static)
        assert "--error" in status.classes

