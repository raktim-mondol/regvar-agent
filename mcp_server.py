"""MCP server exposing the AlphaGenome scoring tools.

This is the second, more reusable way to "integrate LLM agents into coding
workflows": instead of an in-process loop, expose the same validated tools over
the Model Context Protocol so *any* MCP-capable client -- Claude Code in the
terminal, an IDE agent, a notebook assistant -- can call AlphaGenome with the
caching/throttling/coordinate guarantees baked in. The lab gets one server; every
agent surface reuses it.

    pip install "mcp[cli]"
    export ALPHAGENOME_API_KEY=...
    python mcp_server.py            # stdio transport, ready for an MCP client

Register in an MCP client config, e.g.:
    {
      "mcpServers": {
        "alphagenome": {"command": "python", "args": ["/path/to/mcp_server.py"]}
      }
    }
"""

from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from regvar.tools import tool_list_supported_assays, tool_score_regulatory_variant, tool_score_variants_batch

mcp = FastMCP("alphagenome-regvar")


@mcp.tool()
def list_supported_assays() -> dict:
    """List functional genomic assays AlphaGenome can predict and the wet-lab
    assay each maps to (ATAC, RNA-seq, ChIP histone, Hi-C/pcHi-C, ...)."""
    return tool_list_supported_assays()


@mcp.tool()
def score_regulatory_variant(
    chromosome: str,
    position: int,
    ref: str,
    alt: str,
    assays: list[str] | None = None,
    ontology_terms: list[str] | None = None,
    top_n: int = 10,
) -> dict:
    """Predict the functional effect of a regulatory variant with AlphaGenome.

    position is 1-based (VCF convention). Returns the strongest predicted
    effects across the requested assays in the given tissue ontology
    (defaults to prostate gland). Results are cached and rate-limited.
    """
    return tool_score_regulatory_variant(
        chromosome=chromosome, position=position, ref=ref, alt=alt,
        assays=assays, ontology_terms=ontology_terms, top_n=top_n,
    )


@mcp.tool()
def score_variants_batch(
    variants: list[dict],
    assays: list[str] | None = None,
    ontology_terms: list[str] | None = None,
    top_n: int = 10,
) -> dict:
    """Score multiple regulatory variants in a single batch call.

    Each variant dict must have keys: chromosome (e.g. 'chr8'), position
    (1-based), ref, alt. Returns per-variant ranked effects across the
    requested assays. Prefer this over calling score_regulatory_variant
    repeatedly when you have a list of candidates.
    """
    return tool_score_variants_batch(
        variants=variants, assays=assays,
        ontology_terms=ontology_terms, top_n=top_n,
    )


if __name__ == "__main__":
    mcp.run()
