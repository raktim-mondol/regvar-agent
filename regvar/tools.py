"""The tool boundary an LLM agent calls.

This is the "integrating LLM-based agents into coding workflows" half of the
project. Rather than letting a model free-form generate AlphaGenome API calls in
a loop (expensive, and easy to get coordinates wrong), we expose a small set of
*validated* tools. The agent reasons about which variants matter and how to
interpret them; this layer guarantees every call is correct, cached, and
rate-limited.

The same callables are reused by the MCP server (mcp_server.py), so a single
implementation backs both an in-process agent loop and any external MCP client
such as Claude Code.

TOOL_SCHEMAS          — Anthropic Messages API format (legacy, kept for MCP)
OPENAI_TOOL_SCHEMAS   — OpenAI / DeepSeek function-calling format (used by agent.py)
"""

from __future__ import annotations

import json
from typing import Any

from .alphagenome_client import AlphaGenomeClient, PROSTATE_GLAND
from .variants import MAX_VARIANTS, rank_by_effect, validate_chromosome

# Module-level client so cache/throttle state is shared across tool calls.
_client: AlphaGenomeClient | None = None


def get_client() -> AlphaGenomeClient:
    global _client
    if _client is None:
        _client = AlphaGenomeClient()
    return _client


# --- tool implementations -------------------------------------------------

def tool_list_supported_assays() -> dict[str, Any]:
    """Return the assays AlphaGenome can score and the lab assay each maps to."""
    return {"assays": get_client().supported_assays()}


def tool_score_regulatory_variant(
    chromosome: str,
    position: int,
    ref: str,
    alt: str,
    assays: list[str] | None = None,
    ontology_terms: list[str] | None = None,
    top_n: int = 10,
) -> dict[str, Any]:
    """Score a single regulatory variant across assays and return the strongest
    predicted effects as records the model can read directly."""
    chromosome = validate_chromosome(chromosome)
    assays = assays or ["ATAC-seq", "RNA-seq", "ChIP-seq histone", "Hi-C / pcHi-C"]
    ontology_terms = ontology_terms or [PROSTATE_GLAND]
    tidy = get_client().score_variant(
        chromosome=chromosome, position=position, ref=ref, alt=alt,
        assays=assays, ontology_terms=ontology_terms,
    )
    ranked = rank_by_effect(tidy, top_n=top_n)
    return {
        "variant": f"{chromosome}:{position}:{ref}>{alt}",
        "assays_scored": assays,
        "n_total_scores": int(len(tidy)),
        "top_effects": ranked.to_dict(orient="records"),
    }


def tool_score_variants_batch(
    variants: list[dict[str, Any]],
    assays: list[str] | None = None,
    ontology_terms: list[str] | None = None,
    top_n: int = 10,
) -> dict[str, Any]:
    """Score multiple regulatory variants in a single batch call.

    Each element of *variants* must have keys: chromosome, position, ref, alt.
    Returns a dict with per-variant ranked effects, which is much cheaper than
    N individual ``score_regulatory_variant`` calls when the agent needs to
    score an entire candidate list.
    """
    for v in variants:
        v["chromosome"] = validate_chromosome(v["chromosome"])
    if len(variants) > MAX_VARIANTS:
        raise ValueError(
            f"Batch contains {len(variants)} variants, exceeding the limit of "
            f"{MAX_VARIANTS}. Reduce the batch size."
        )
    assays = assays or ["ATAC-seq", "RNA-seq", "ChIP-seq histone", "Hi-C / pcHi-C"]
    ontology_terms = ontology_terms or [PROSTATE_GLAND]
    tidy = get_client().score_variants(
        variants, assays=assays, ontology_terms=ontology_terms,
    )
    results = {}
    if not tidy.empty:
        for variant_id, group in tidy.groupby("query_variant"):
            ranked = rank_by_effect(group, top_n=top_n)
            results[variant_id] = ranked.to_dict(orient="records")
    return {
        "variants_scored": len(results),
        "assays_scored": assays,
        "n_total_scores": int(len(tidy)),
        "effects_by_variant": results,
    }


# --- Anthropic tool-use schemas (kept for MCP server compatibility) --------

TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "list_supported_assays",
        "description": (
            "List the functional genomic assays AlphaGenome can predict and the "
            "wet-lab assay each corresponds to. Call this first if unsure which "
            "assays to request."
        ),
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "score_regulatory_variant",
        "description": (
            "Predict the functional consequence of a single non-coding/regulatory "
            "variant using AlphaGenome, across the requested assays, in the given "
            "tissue ontology. Returns the strongest predicted effects. Position is "
            "1-based (VCF convention). Use this to test whether a candidate variant "
            "in an enhancer/promoter is predicted to alter accessibility, histone "
            "marks, expression, or enhancer-promoter contact."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "chromosome": {"type": "string", "description": "e.g. 'chr8'"},
                "position": {"type": "integer", "description": "1-based position"},
                "ref": {"type": "string", "description": "reference allele"},
                "alt": {"type": "string", "description": "alternate allele"},
                "assays": {
                    "type": "array", "items": {"type": "string"},
                    "description": (
                        "Subset of: ATAC-seq, DNase-seq, RNA-seq, "
                        "'ChIP-seq histone', 'ChIP-seq TF', 'Hi-C / pcHi-C', splicing"
                    ),
                },
                "ontology_terms": {
                    "type": "array", "items": {"type": "string"},
                    "description": "UBERON/CL terms; defaults to prostate gland.",
                },
                "top_n": {"type": "integer", "description": "How many top effects to return."},
            },
            "required": ["chromosome", "position", "ref", "alt"],
        },
    },
    {
        "name": "score_variants_batch",
        "description": (
            "Score multiple regulatory variants in one call. Prefer this over "
            "calling score_regulatory_variant repeatedly when you have a list of "
            "candidates. Each variant needs chromosome (e.g. 'chr8'), position "
            "(1-based), ref, and alt. Returns per-variant ranked effects."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "variants": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "chromosome": {"type": "string"},
                            "position": {"type": "integer"},
                            "ref": {"type": "string"},
                            "alt": {"type": "string"},
                        },
                        "required": ["chromosome", "position", "ref", "alt"],
                    },
                    "description": "List of variant dicts to score.",
                },
                "assays": {
                    "type": "array", "items": {"type": "string"},
                    "description": (
                        "Subset of: ATAC-seq, DNase-seq, RNA-seq, "
                        "'ChIP-seq histone', 'ChIP-seq TF', 'Hi-C / pcHi-C', splicing"
                    ),
                },
                "ontology_terms": {
                    "type": "array", "items": {"type": "string"},
                    "description": "UBERON/CL terms; defaults to prostate gland.",
                },
                "top_n": {"type": "integer", "description": "How many top effects per variant."},
            },
            "required": ["variants"],
        },
    },
]

# --- OpenAI / DeepSeek function-calling schemas ----------------------------
# Passed in the `tools` field of the OpenAI chat.completions API.
# Structure: {"type": "function", "function": {"name": ..., "description": ..., "parameters": ...}}

OPENAI_TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "type": "function",
        "function": {
            "name": "list_supported_assays",
            "description": (
                "List the functional genomic assays AlphaGenome can predict and the "
                "wet-lab assay each corresponds to. Call this first if unsure which "
                "assays to request."
            ),
            "parameters": {"type": "object", "properties": {}, "required": []},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "score_regulatory_variant",
            "description": (
                "Predict the functional consequence of a single non-coding/regulatory "
                "variant using AlphaGenome, across the requested assays, in the given "
                "tissue ontology. Returns the strongest predicted effects. Position is "
                "1-based (VCF convention). Use this to test whether a candidate variant "
                "in an enhancer/promoter is predicted to alter accessibility, histone "
                "marks, expression, or enhancer-promoter contact."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "chromosome": {"type": "string", "description": "e.g. 'chr8'"},
                    "position": {"type": "integer", "description": "1-based position"},
                    "ref": {"type": "string", "description": "reference allele"},
                    "alt": {"type": "string", "description": "alternate allele"},
                    "assays": {
                        "type": "array", "items": {"type": "string"},
                        "description": (
                            "Subset of: ATAC-seq, DNase-seq, RNA-seq, "
                            "'ChIP-seq histone', 'ChIP-seq TF', 'Hi-C / pcHi-C', splicing"
                        ),
                    },
                    "ontology_terms": {
                        "type": "array", "items": {"type": "string"},
                        "description": "UBERON/CL terms; defaults to prostate gland.",
                    },
                    "top_n": {"type": "integer", "description": "How many top effects to return."},
                },
                "required": ["chromosome", "position", "ref", "alt"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "score_variants_batch",
            "description": (
                "Score multiple regulatory variants in one call. Prefer this over "
                "calling score_regulatory_variant repeatedly when you have a list of "
                "candidates. Each variant needs chromosome (e.g. 'chr8'), position "
                "(1-based), ref, and alt. Returns per-variant ranked effects."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "variants": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "chromosome": {"type": "string"},
                                "position": {"type": "integer"},
                                "ref": {"type": "string"},
                                "alt": {"type": "string"},
                            },
                            "required": ["chromosome", "position", "ref", "alt"],
                        },
                        "description": "List of variant dicts to score.",
                    },
                    "assays": {
                        "type": "array", "items": {"type": "string"},
                        "description": (
                            "Subset of: ATAC-seq, DNase-seq, RNA-seq, "
                            "'ChIP-seq histone', 'ChIP-seq TF', 'Hi-C / pcHi-C', splicing"
                        ),
                    },
                    "ontology_terms": {
                        "type": "array", "items": {"type": "string"},
                        "description": "UBERON/CL terms; defaults to prostate gland.",
                    },
                    "top_n": {"type": "integer", "description": "How many top effects per variant."},
                },
                "required": ["variants"],
            },
        },
    },
]

# Dispatch table the agent loop uses to execute a tool call by name.
def _dispatch_list_assays(**_kw: Any) -> dict[str, Any]:
    return tool_list_supported_assays()


def _dispatch_score_variant(**kw: Any) -> dict[str, Any]:
    return tool_score_regulatory_variant(**kw)


def _dispatch_score_batch(**kw: Any) -> dict[str, Any]:
    return tool_score_variants_batch(**kw)


TOOL_DISPATCH: dict[str, Any] = {
    "list_supported_assays": _dispatch_list_assays,
    "score_regulatory_variant": _dispatch_score_variant,
    "score_variants_batch": _dispatch_score_batch,
}


def run_tool(name: str, arguments: dict[str, Any]) -> str:
    """Execute a tool by name and return a JSON string (the format tool_result
    blocks expect)."""
    if name not in TOOL_DISPATCH:
        return json.dumps({"error": f"unknown tool {name}"})
    try:
        return json.dumps(TOOL_DISPATCH[name](**arguments), default=str)
    except Exception as exc:  # noqa: BLE001 - report cleanly back to the model
        return json.dumps({"error": str(exc)})
