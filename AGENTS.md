# AGENTS.md — regvar-agent

Quick-ramp context for OpenCode and other agents working in this repo.

---

## Commands

```bash
# install as a real CLI tool (adds `regvar` to PATH)
pip install -e .

# ── single variant ──────────────────────────────────────────────────────────
regvar score chr8 127401060 G T                          # score & print table
regvar score chr8 127401060 G T --output-tsv scores.tsv  # save TSV
regvar score chr8 127401060 G T --output-png effects.png # save barplot
regvar score chr8 127401060 G T --json                   # raw JSON

# ── barplot ──────────────────────────────────────────────────────────────────
regvar plot effects chr8 127401060 G T --output effects.png  # save PNG
regvar plot effects chr8 127401060 G T --show                # open GUI window
regvar plot effects chr8 127401060 G T \
    --from-tsv scores.tsv --assay-filter ATAC-seq --output atac.png

# ── agent triage loop ────────────────────────────────────────────────────────
regvar agent run examples/candidate_variants.tsv          # full agent loop
regvar agent run examples/candidate_variants.tsv --dry-run # validate only
regvar agent run examples/candidate_variants.tsv \
    -o report.md --output-tsv scores.tsv

# ── interactive agent chat ───────────────────────────────────────────────────
regvar agent chat                                          # start REPL
regvar agent chat --tissue UBERON:0002367 --top-n 15

# ── utilities ────────────────────────────────────────────────────────────────
regvar assays                                              # list assays
regvar tui                                                 # Textual TUI
regvar --version

# ── backward-compat (still works, delegates to agent run) ───────────────────
regvar run examples/candidate_variants.tsv

# ── python -m equivalents (no pip install needed) ───────────────────────────
python -m regvar score chr8 127401060 G T
python -m regvar agent run examples/candidate_variants.tsv
python -m regvar.agent examples/candidate_variants.tsv     # legacy, unchanged

# test (no keys, no network)
pytest -q
```

Tests must be run from the **repo root** (`/alphaGenome/`), not from inside `regvar/`.  
`python -m regvar` also requires the repo root as cwd so the package import resolves.

---

## Environment / secrets

- Keys live in `.env` at repo root — loaded automatically via `python-dotenv` at import time in both `agent.py` and `alphagenome_client.py`.
- `.env` is gitignored. `.env.example` is committed as a template.
- Shell `export` still works as a fallback if dotenv is absent.
- Required keys: `ALPHAGENOME_API_KEY`, `DEEPSEEK_API_KEY`.

---

## Architecture in one paragraph

`regvar/agent.py` runs a DeepSeek v4 Pro tool-use loop (OpenAI-compatible API, `https://api.deepseek.com`). The loop calls tools defined in `regvar/tools.py`, which dispatches to `AlphaGenomeClient` in `regvar/alphagenome_client.py`. The client owns **all** coordinate handling (1-based VCF positions — never shift them), caching (`.pkl` files in `~/.cache/regvar/`), throttle, and retry. `mcp_server.py` at repo root wraps the same two tool functions from `tools.py` for MCP clients.

---

## Critical gotchas

- **Positions are always 1-based (VCF convention).** The client passes them directly to `genome.Variant.position`. Never add or subtract 1 anywhere else.
- **`TOOL_SCHEMAS`** in `tools.py` is the legacy Anthropic format (kept for MCP). **`OPENAI_TOOL_SCHEMAS`** is the format used by `agent.py` and DeepSeek. Don't mix them.
- `_client` in `tools.py` is a module-level singleton — cache and rate-limit state is shared across all tool calls in a session. Don't instantiate a second `AlphaGenomeClient` inside tools.
- AlphaGenome API is rate-limited and shared. The client throttles to 1 call/sec by default. Don't bypass `_call_with_retry`.
- `mcp_server.py` must be run from the repo root so `from regvar.tools import ...` resolves.

---

## Package boundaries

```
regvar/               Python package (import as `from regvar import ...`)
  __init__.py         exports: AlphaGenomeClient, ClientConfig, ASSAY_TO_OUTPUT,
                               CandidateVariant, read_candidates_tsv, rank_by_effect
  __main__.py         `python -m regvar` entry — delegates to cli.py
  cli.py              click + rich CLI: run / score / assays commands
  alphagenome_client.py   sole owner of AlphaGenome API calls + caching
  variants.py             TSV I/O, coordinate utils, effect ranking
  tools.py                tool implementations + both schema formats
  agent.py                DeepSeek agent loop (also usable via python -m regvar.agent)
mcp_server.py         MCP server — thin wrapper over regvar.tools, not part of the package
examples/             input data (swap with real VCF candidates)
tests/                unit tests — pure Python, no API calls
```

---

## Adding a new assay

1. Add entry to `ASSAY_TO_OUTPUT` dict in `alphagenome_client.py`.
2. Add the same label string to the `description` field in both `TOOL_SCHEMAS` and `OPENAI_TOOL_SCHEMAS` in `tools.py`.
3. No other files need changing.

---

## Testing notes

- All 6 tests are pure-Python unit tests (no mocks, no API calls, no fixtures beyond `tmp_path`).
- Tests cover: coordinate conversion, TSV parsing, missing-column error, effect ranking, fallback scoring.
- There are no integration tests that hit the real API — manual smoke test: run the agent with both keys set and confirm tool calls appear in stderr.
