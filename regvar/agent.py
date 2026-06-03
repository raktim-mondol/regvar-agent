"""Agentic variant-triage workflow.

Demonstrates an LLM agent driving the AlphaGenome tools to do the kind of work
the role describes: "streamline complex data processing... and accelerate
hypothesis generation", with the wet lab in mind.

Flow:
    1. We hand the agent a list of candidate regulatory variants (e.g. SNPs
       falling in stromal-fibroblast ATAC/ChIP peaks near a prostate-cancer gene).
    2. The agent calls `score_regulatory_variant` for each (real tool execution,
       not hallucinated numbers).
    3. The agent synthesises: which variants are predicted to be functional, on
       which assay/cell type, the likely regulatory mechanism, and a concrete,
       prioritised set of wet-lab validation experiments.

Uses DeepSeek v4 Pro via the OpenAI-compatible API.

    export DEEPSEEK_API_KEY=...
    export ALPHAGENOME_API_KEY=...
    python -m regvar.agent examples/candidate_variants.tsv
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any, Callable

from .tools import OPENAI_TOOL_SCHEMAS, run_tool
from .variants import read_candidates_tsv

DEFAULT_MODEL = "deepseek-v4-pro"
DEEPSEEK_BASE_URL = "https://api.deepseek.com"
MAX_TURNS = 24

SYSTEM_PROMPT = """You are a computational genomics assistant embedded in an \
epigenetics lab studying the prostate tumour microenvironment and stromal \
fibroblasts. You have tools that run AlphaGenome variant-effect predictions \
across functional assays the lab generates (ATAC-seq, RNA-seq, ChIP-seq histone \
marks, Hi-C/pcHi-C contact). You can score a single variant with \
score_regulatory_variant, or score multiple variants at once with \
score_variants_batch (preferred when scoring a full candidate list).

For the candidate regulatory variants given:
  1. Score each with the tool. Do not invent scores; use the tool's output.
  2. Rank variants by predicted functional impact and note which assay and cell \
type drives the signal.
  3. For the strongest candidates, state the likely regulatory mechanism \
(e.g. disrupted accessibility -> loss of enhancer activity -> reduced target \
expression; or altered enhancer-promoter contact).
  4. Propose a prioritised, concrete validation plan a wet-lab scientist could \
run (e.g. luciferase reporter, CRISPRi of the element, allele-specific ATAC, \
4C/capture-C around the locus), tying each experiment to the specific prediction \
it would confirm.

Be explicit about uncertainty and that AlphaGenome predictions are hypotheses \
for experimental testing, not clinical conclusions."""


def build_task_message(candidates) -> str:
    lines = ["Candidate regulatory variants:"]
    try:
        from .annotation import get_annotation
    except ImportError:
        get_annotation = None  # type: ignore[assignment]
    for c in candidates:
        note_parts = []
        if c.region_id or c.note:
            note_parts.append(c.region_id)
            if c.note:
                note_parts.append(c.note)
        if get_annotation is not None:
            ann = get_annotation(c.vcf_id)
            if ann is not None:
                snippet = ann.to_note()
                if snippet:
                    note_parts.append(snippet)
        ctx = f"  ({'; '.join(p for p in note_parts if p)})" if note_parts else ""
        lines.append(f"  - {c.vcf_id}{ctx}")
    lines.append("\nScore them and produce the prioritised analysis and validation plan.")
    return "\n".join(lines)


def _build_system_prompt(
    assay_hint: list[str] | None = None,
    tissue_hint: list[str] | None = None,
    top_n_hint: int | None = None,
) -> str:
    """Build the system prompt, appending user-supplied preference hints."""
    parts = [SYSTEM_PROMPT]
    addenda: list[str] = []
    if assay_hint:
        addenda.append(f"Preferred assays: {', '.join(assay_hint)}.")
    if tissue_hint:
        addenda.append(f"Preferred tissue ontology terms: {', '.join(tissue_hint)}.")
    if top_n_hint is not None:
        addenda.append(f"Return the top {top_n_hint} effects per variant.")
    if addenda:
        parts.append("\n\nUser preferences:\n" + "\n".join(f"  - {a}" for a in addenda))
    return "".join(parts)


# Type alias for the optional progress callback.
OnToolCall = Callable[[str, dict, str], None] | None


def run_agent_turn(
    messages: list,
    client: "OpenAI",
    model: str = DEFAULT_MODEL,
    tools: list | None = None,
    max_tokens: int = 8192,
    reasoning_effort: str | None = "high",
    on_tool_call: OnToolCall = None,
    verbose: bool = False,
    max_turns: int = MAX_TURNS,
) -> str:
    """Run one complete tool-use cycle: API call → dispatch tools → repeat until
    the model returns a final (non-tool-call) response.

    Mutates *messages* in place (appends assistant + tool messages).
    Returns the assistant's final text content.

    Parameters
    ----------
    messages : list
        Conversation history (mutated in place).
    client : OpenAI
        An initialised OpenAI client pointed at the target API.
    model : str
        LLM model name.
    tools : list, optional
        Tool schemas to pass to the API. Defaults to OPENAI_TOOL_SCHEMAS.
    max_tokens : int
        Max tokens for the response.
    reasoning_effort : str or None
        Reasoning effort level. Pass None to disable extended thinking.
    on_tool_call : callable, optional
        ``fn(name, args, result_json)`` invoked after each tool execution.
    verbose : bool
        If True, print tool calls to stderr.
    max_turns : int
        Maximum API calls before giving up.
    """
    if tools is None:
        tools = OPENAI_TOOL_SCHEMAS

    for _ in range(max_turns):
        kwargs: dict[str, Any] = dict(
            model=model,
            messages=messages,
            tools=tools,
            tool_choice="auto",
            max_tokens=max_tokens,
        )
        if reasoning_effort is not None:
            kwargs["reasoning_effort"] = reasoning_effort
            kwargs["extra_body"] = {"thinking": {"type": "enabled"}}

        resp = client.chat.completions.create(**kwargs)

        msg = resp.choices[0].message
        messages.append(msg.model_dump(exclude_unset=True))

        if resp.choices[0].finish_reason != "tool_calls":
            return msg.content or ""

        for tc in msg.tool_calls or []:
            fn_name = tc.function.name
            try:
                fn_args = json.loads(tc.function.arguments)
            except json.JSONDecodeError:
                fn_args = {}
            if verbose:
                print(f"  [tool] {fn_name}({fn_args})", file=sys.stderr)
            result = run_tool(fn_name, fn_args)
            if on_tool_call is not None:
                on_tool_call(fn_name, fn_args, result)
            messages.append({
                "role": "tool",
                "tool_call_id": tc.id,
                "content": result,
            })

    return "Stopped: reached MAX_TURNS without a final answer."


def run_agent(
    candidates,
    model: str = DEFAULT_MODEL,
    verbose: bool = True,
    on_tool_call: OnToolCall = None,
    assay_hint: list[str] | None = None,
    tissue_hint: list[str] | None = None,
    top_n_hint: int | None = None,
) -> str:
    """Run the OpenAI-format tool-use loop and return the agent's final synthesis text.

    Parameters
    ----------
    candidates : list[CandidateVariant]
        Variants to score.
    model : str
        LLM model name.
    verbose : bool
        If True, print tool calls to stderr (legacy behaviour).
    on_tool_call : callable, optional
        ``fn(name, args, result_json)`` invoked after each tool execution.
        Used by the CLI for progress display and result collection.
    assay_hint, tissue_hint, top_n_hint : optional
        User preferences from CLI flags — appended to the system prompt so the
        agent can take them into account.
    """
    try:
        from openai import OpenAI
    except ImportError as exc:
        raise SystemExit("pip install openai") from exc

    api_key = os.environ.get("DEEPSEEK_API_KEY")
    if not api_key:
        raise SystemExit("Set DEEPSEEK_API_KEY environment variable.")

    client = OpenAI(api_key=api_key, base_url=DEEPSEEK_BASE_URL)

    system_prompt = _build_system_prompt(assay_hint, tissue_hint, top_n_hint)

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": build_task_message(candidates)},
    ]

    return run_agent_turn(
        messages=messages,
        client=client,
        model=model,
        on_tool_call=on_tool_call,
        verbose=verbose,
    )


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if not argv:
        print("usage: python -m regvar.agent <candidates.tsv>", file=sys.stderr)
        return 2
    candidates = read_candidates_tsv(argv[0])
    print(f"Loaded {len(candidates)} candidate variants. Running agent...\n", file=sys.stderr)
    answer = run_agent(candidates)
    print(answer)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
