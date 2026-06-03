# regvar-agent: Improvement & Feature Plan

> Prioritised roadmap for code quality, new features, robustness, testing, and advanced capabilities.
> Phase 1 is complete; Phases 2–5 are pending.

---

## Phase 1: Code Quality Foundation ✅ DONE

Eliminated duplication and dead code that blocked clean feature work.

| Task | Status | What changed |
|------|--------|-------------|
| **1A** Shared TSV serialization | ✅ | Added `write_scores_tsv()` to `variants.py`; replaced 3 duplicate sites in `cli.py` and `tui.py` |
| **1B** Reusable agent turn function | ✅ | Extracted `run_agent_turn()` in `agent.py`; both `run_agent()` and `agent_chat` now call it |
| **1C** Dead code removal | ✅ | Deleted unused `iter_chunks()` and duplicate `Path` import |
| **1D** Centralized `.env` loading | ✅ | Moved dotenv to `__init__.py`; removed from `agent.py` and `alphagenome_client.py` |
| **1E** Rename `compile` parameter | ✅ | `compile` → `compile_pdf` in `report.py` to avoid shadowing Python builtin |
| **1F** Type-annotate TOOL_DISPATCH | ✅ | Replaced lambdas with named functions; added return type to `run_tool` |

---

## Phase 2: High-Impact Features (P1 Roadmap Items)

### 2A. VCF Input Support

**Problem:** Currently only accepts a custom TSV format. Most real-world GWAS / fine-mapping outputs are in VCF/VCF.gz, requiring manual conversion.

**Change:**
1. Add `read_candidates_vcf(path)` to `regvar/variants.py` using lazy-imported `cyvcf2` (fast, C-backed)
2. Handle multi-allelic sites by expanding ALT alleles into separate `CandidateVariant` rows
3. Auto-prepend `chr` prefix if missing (AlphaGenome expects `chr8` not `8`)
4. Use VCF `ID` field as `region_id` when present (often contains rsIDs)
5. Auto-detect `.vcf` / `.vcf.gz` by extension in `agent_run`, `run`, and TUI
6. Add `vcf` optional dependency group in `pyproject.toml`

**Files to modify:**
- `regvar/variants.py` — add `read_candidates_vcf()`
- `regvar/cli.py` — detect file extension in `agent_run`, `run`, dry-run
- `regvar/tui.py` — auto-detect in `RunAgentTab._do_run`
- `pyproject.toml` — add `vcf = ["cyvcf2"]` optional dependency

**New file:**
- `tests/test_vcf.py` — unit tests with inline VCF fixture

---

### 2B. Expose `score_variants()` as a Batch Tool

**Problem:** `AlphaGenomeClient.score_variants()` exists but isn't wired to the tool layer. The agent must score variants one at a time — N tool calls for N variants.

**Change:**
1. Add `tool_score_variants_batch()` in `regvar/tools.py` accepting a list of variant dicts
2. Add OpenAI and Anthropic tool schemas for the new function
3. Register in `TOOL_DISPATCH`
4. Expose in `mcp_server.py` via `@mcp.tool()`
5. Update `SYSTEM_PROMPT` in `agent.py` to mention the batch tool

**Files to modify:**
- `regvar/tools.py` — new tool function, schemas, dispatch entry
- `regvar/agent.py` — update system prompt
- `mcp_server.py` — add MCP wrapper

---

### 2C. Progress Bar for Agent Runs

**Problem:** For large variant sets, the agent loop runs silently for minutes with no visual progress.

**Change:**
Wrap `_run_agent_core` in `rich.progress.Progress` (already a dependency). Use variant count for a determinate bar; update task description per tool call. The TUI already has its own `RichLog`-based progress — this only affects the CLI path.

**Files to modify:**
- `regvar/cli.py` — wrap agent loop with `rich.progress.Progress`

---

### 2D. `--assays all` Shortcut

**Problem:** `DNase-seq`, `ChIP-seq TF`, and `splicing` are defined in `ASSAY_TO_OUTPUT` but excluded from defaults. No way to request all assays without listing them individually.

**Change:**
Add `_parse_assays(raw: str | None) -> list[str] | None` helper in `cli.py` that:
- `None` → returns `None` (use defaults)
- `"all"` → returns `list(ASSAY_TO_OUTPUT.keys())`
- Comma-separated string → splits and strips

Use in all commands with `--assays` flag: `score`, `plot effects`, `agent run`, `agent chat`.

**Files to modify:**
- `regvar/cli.py` — add helper, apply everywhere `--assays` is parsed

---

### 2E. Structured JSON Output

**Problem:** Only markdown reports and ad-hoc TSV output. No machine-readable format for downstream programmatic use or database ingestion.

**Change:**
Add `--output-json` flag to `agent run` and `score` commands. Build a structured JSON object:

```json
{
  "metadata": {
    "timestamp": "2026-06-03T14:30:00",
    "model": "deepseek-v4-pro",
    "assays": ["ATAC-seq", "RNA-seq"],
    "tissue": ["UBERON:0002367"]
  },
  "variants": [...],
  "scores": [...],
  "report": "# Agent Analysis\n..."
}
```

Reuses the existing `tool_call_log` collected in `_run_agent_core`.

**Files to modify:**
- `regvar/cli.py` — add `--output-json` to `agent_run` and `score`

---

## Phase 3: Robustness & Input Validation

### 3A. Chromosome Name Validation

**Problem:** No validation on chromosome names — any string is accepted, which could produce silent AlphaGenome API errors.

**Change:**
Add `validate_chromosome(name: str) -> str` to `regvar/variants.py`:
- Accept `chr1`–`chr22`, `chrX`, `chrY`, `chrM` (and without `chr` prefix, prepending it)
- Reject anything else with a clear error
- Return the normalized form with `chr` prefix

Call from TSV/VCF readers, `tool_score_regulatory_variant`, and TUI input validation.

**Files to modify:**
- `regvar/variants.py` — add `validate_chromosome()`
- `regvar/tools.py` — validate in `tool_score_regulatory_variant`
- `regvar/tui.py` — validate in `ScoreVariantTab.handle_score`

---

### 3B. Max Variant Count Guard

**Problem:** No upper limit on variant count. Accidentally passing a 1M-variant file would try to score all of them sequentially.

**Change:**
- Add `MAX_VARIANTS = 500` constant (env-overridable via `REGVAR_MAX_VARIANTS`) in `variants.py`
- Readers warn when count exceeds limit
- Batch tool refuses beyond limit
- CLI commands print a warning and require `--force` to proceed

**Files to modify:**
- `regvar/variants.py` — add constant and warning in readers
- `regvar/cli.py` — add `--force` flag to `agent_run`
- `regvar/tools.py` — add guard in batch tool

---

### 3C. Cache Management CLI

**Problem:** The pickle cache at `~/.cache/regvar/` has no management interface. Users can't inspect, clear, or check its size.

**Change:**
Add `regvar cache` command group:

```
regvar cache info              # show dir, total size, entry count
regvar cache clear             # delete all .pkl files
regvar cache clear --older-than 30  # delete entries older than 30 days
```

Pure filesystem operations reading `ClientConfig.cache_dir`.

**Files to modify:**
- `regvar/cli.py` — add `cache` command group

---

### 3D. Replace Pickle Cache with Parquet

**Problem:** Pickle-based caching is a deserialization risk if the cache directory is compromised.

**Change:**
Switch `_save_cache` / `_load_cache` in `alphagenome_client.py`:
- Write: always use `DataFrame.to_parquet()`
- Read: try parquet first, fall back to pickle for backward compatibility
- After a transition period, drop pickle support

**Files to modify:**
- `regvar/alphagenome_client.py` — update cache methods

---

## Phase 4: Test Coverage Expansion

Currently 58 tests exist for `variants.py`, `mofa_view.py`, `report.py`, and `tui.py`. The following modules have **zero test coverage**:

### 4A. Tests for `tools.py`

| Test | What it covers |
|------|---------------|
| `test_tool_list_assays` | Returns expected dict from `ASSAY_TO_OUTPUT` |
| `test_tool_score_variant_structure` | Mock client → verify result has `variant`, `assays_scored`, `n_total_scores`, `top_effects` |
| `test_run_tool_known` | Dispatches correctly for known tool names |
| `test_run_tool_unknown` | Returns JSON error string for unknown tool |
| `test_run_tool_exception` | Catches exceptions, returns JSON error |
| `test_dispatch_matches_schemas` | `TOOL_DISPATCH` keys match schema names |

**Mock strategy:** `unittest.mock.patch` on `tools.get_client` → returns mock `AlphaGenomeClient` with known DataFrame.

### 4B. Tests for `agent.py`

| Test | What it covers |
|------|---------------|
| `test_build_task_message` | Produces expected format from `CandidateVariant` list |
| `test_build_system_prompt_no_hints` | Base prompt only |
| `test_build_system_prompt_with_hints` | Assay/tissue/top_n hints appended |
| `test_run_agent_single_turn` | Mock OpenAI → returns final answer immediately |
| `test_run_agent_multi_turn` | Mock OpenAI → tool call → dispatch → final answer |
| `test_run_agent_max_turns` | Mock OpenAI → always tool_calls → returns MAX_TURNS message |
| `test_on_tool_call_callback` | Callback invoked with correct args |

**Mock strategy:** Mock `OpenAI` class, return scripted responses.

### 4C. Tests for `cli.py`

| Test | What it covers |
|------|---------------|
| `test_version` | `regvar --version` prints version |
| `test_help` | `regvar --help` shows help |
| `test_score_help` | `regvar score --help` shows options |
| `test_dry_run_valid` | `agent run --dry-run` with valid TSV succeeds |
| `test_dry_run_invalid_path` | `agent run --dry-run` with bad path fails |
| `test_assays` | `regvar assays` prints assay table |
| `test_write_scores_tsv` | Produces correct TSV output |
| `test_parse_assays_all` | `_parse_assays("all")` returns full list |
| `test_parse_assays_csv` | `_parse_assays("ATAC-seq,RNA-seq")` splits correctly |
| `test_parse_assays_none` | `_parse_assays(None)` returns None |

**Mock strategy:** Click's `CliRunner` for command testing.

### 4D. Integration Test Scaffold

Create `tests/test_integration.py` gated behind `RUN_INTEGRATION_TESTS=1` + `ALPHAGENOME_API_KEY`:

```python
@pytest.mark.skipif(
    not os.environ.get("RUN_INTEGRATION_TESTS"),
    reason="Set RUN_INTEGRATION_TESTS=1 and ALPHAGENOME_API_KEY to run"
)
class TestIntegration:
    def test_score_known_variant(self):
        # rs6983267 at chr8:127401060 G>T
        result = tool_score_regulatory_variant("chr8", 127401060, "G", "T")
        assert result["n_total_scores"] > 0
        assert "chr8:127401060:G>T" in result["variant"]
```

**New files:**
- `tests/test_tools.py`
- `tests/test_agent.py`
- `tests/test_cli.py`
- `tests/test_integration.py`

---

## Phase 5: Advanced Features (P2–P3)

### 5A. Config File Support

TOML config at `~/.config/regvar/config.toml` or `./regvar.toml`:

```toml
[defaults]
assays = ["ATAC-seq", "RNA-seq", "ChIP-seq histone"]
tissue = ["UBERON:0002367"]
model = "deepseek-v4-pro"
top_n = 10

[cache]
dir = "~/.cache/regvar"
max_variants = 500
```

Precedence: CLI flags > config file > env vars > defaults.

**Files to modify:**
- `regvar/__init__.py` — add `get_config()` using `tomllib` (3.11+) or `tomli` (3.10)
- `regvar/cli.py` — read config defaults
- `pyproject.toml` — add `tomli` for Python <3.11

### 5B. Parallel Scoring with Semaphore

Add `max_concurrent` to `ClientConfig` (default 1). Rewrite `score_variants()` with `concurrent.futures.ThreadPoolExecutor` + `Semaphore`. Per-call throttle remains as a floor. Backward-compatible (default concurrency=1 = current behaviour).

**Files to modify:**
- `regvar/alphagenome_client.py` — add config field, rewrite `score_variants()`

### 5C. Chat REPL Improvements

Add to `agent chat`:
- `--output` flag → save conversation to markdown
- `--reasoning-effort` parameter (low/medium/high)
- `--save-session` / `--load-session` → JSON persistence of message history

Depends on Phase 1B's extracted `run_agent_turn()`.

**Files to modify:**
- `regvar/cli.py` — add flags and persistence logic to `agent_chat`

### 5D. Variant Annotation Module

New `regvar/annotation.py` with lazy imports (`pyranges` / `pybedtools`):
- Nearest gene + TSS distance from a GTF file
- Regulatory element overlap from a BED file
- rsID lookup from a dbSNP VCF

Annotations feed into the agent's system prompt as additional context per variant.

**New file:**
- `regvar/annotation.py`

**Files to modify:**
- `regvar/variants.py` — optionally attach annotations to `CandidateVariant`
- `regvar/cli.py` — add `--annotate` flag with path to reference files
- `pyproject.toml` — add `annotation = ["pyranges", "pybedtools"]` optional dependency

---

## Dependency Graph

```
Phase 1 ✅ (done)
  │
  ├─→ Phase 3A + 3B (validation guards — do before VCF/batch)
  │     │
  │     └─→ Phase 2A → 2B → 2C → 2D → 2E (features, highest impact first)
  │
  └─→ Phase 4 (tests — can run in parallel with any phase)
        │
        └─→ Phase 3C + 3D (cache improvements)
              │
              └─→ Phase 5 (advanced: 5C > 5B > 5A > 5D)
```

## Recommended Execution Order

1. ~~Phase 1 (all)~~ ✅ **Done**
2. Phase 3A + 3B (chromosome validation + max variant guard)
3. Phase 2A → 2B → 2C → 2D → 2E (features)
4. Phase 4 (test coverage for all modules)
5. Phase 3C + 3D (cache management + parquet migration)
6. Phase 5 (advanced features: 5C > 5B > 5A > 5D)

## Verification

- `pytest -q` from repo root after each phase — all tests must pass
- Phase 2A: test with a small `.vcf` file alongside `.tsv` — both should produce identical candidate lists
- Phase 2B: verify agent makes 1 batch tool call instead of N individual calls
- Phase 4: `pytest -q` should show 80+ tests, all passing without API keys
