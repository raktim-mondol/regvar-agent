"""Plotting helpers for regvar CLI.

Provides a single public function ``plot_effects`` that takes the
``top_effects`` list returned by ``tool_score_regulatory_variant`` and
renders a horizontal barplot of absolute variant effects, colour-coded by
assay type.

    from regvar.plot import plot_effects
    plot_effects(result["top_effects"], title="chr8:127401060 G>T",
                 output_path=Path("effects.png"))
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

# Assay → colour (matplotlib named colours — dark-theme friendly)
_ASSAY_COLOURS: dict[str, str] = {
    "ATAC-seq":        "#4fc3f7",   # light blue
    "DNase-seq":       "#81d4fa",   # pale blue
    "RNA-seq":         "#a5d6a7",   # green
    "ChIP-seq histone": "#ffb74d",  # amber
    "ChIP-seq TF":     "#ff8a65",   # orange
    "Hi-C / pcHi-C":  "#ce93d8",   # lavender
    "splicing":        "#f48fb1",   # pink
}
_FALLBACK_COLOUR = "#90caf9"


def _make_label(rec: dict[str, Any]) -> str:
    """Build a readable y-axis tick label from a score record."""
    parts: list[str] = []
    if rec.get("biosample_name"):
        parts.append(str(rec["biosample_name"])[:40])
    if rec.get("gene_name"):
        parts.append(str(rec["gene_name"]))
    if rec.get("assay"):
        parts.append(f"[{rec['assay']}]")
    return "  ".join(parts) if parts else "unknown"


def plot_effects(
    records: list[dict[str, Any]],
    title: str = "Variant Effects",
    output_path: Path | None = None,
    show: bool = False,
    top_n: int = 20,
    figsize: tuple[float, float] | None = None,
) -> "matplotlib.figure.Figure":  # noqa: F821 (lazy import)
    """Render a horizontal barplot of predicted variant effects.

    Parameters
    ----------
    records:
        List of dicts from ``tool_score_regulatory_variant``'s ``top_effects``.
    title:
        Plot title (typically the variant identifier).
    output_path:
        If given, save the figure to this path (PNG/SVG/PDF auto-detected).
    show:
        If True, call ``plt.show()`` (opens a GUI window).
    top_n:
        Maximum number of tracks to display.
    figsize:
        Override the automatic figure size.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object (useful for embedding in notebooks).
    """
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Run: pip install matplotlib"
        ) from exc

    # ── data prep ────────────────────────────────────────────────────────────
    recs = records[:top_n]
    if not recs:
        raise ValueError("No records to plot.")

    labels   = [_make_label(r) for r in recs]
    values   = [float(r.get("abs_effect", 0.0)) for r in recs]
    raw_vals = [float(r.get("raw_score", 0.0)) for r in recs]
    assays   = [str(r.get("assay", "")) for r in recs]
    colours  = [_ASSAY_COLOURS.get(a, _FALLBACK_COLOUR) for a in assays]

    # Sort by descending absolute effect (records may already be sorted, but be safe).
    order = sorted(range(len(values)), key=lambda i: values[i])
    labels   = [labels[i]  for i in order]
    values   = [values[i]  for i in order]
    raw_vals = [raw_vals[i] for i in order]
    assays   = [assays[i]  for i in order]
    colours  = [colours[i] for i in order]

    # ── figure ───────────────────────────────────────────────────────────────
    matplotlib.rcParams.update({
        "figure.facecolor":  "#0d1117",
        "axes.facecolor":    "#161b22",
        "axes.edgecolor":    "#30363d",
        "axes.labelcolor":   "#c9d1d9",
        "text.color":        "#c9d1d9",
        "xtick.color":       "#8b949e",
        "ytick.color":       "#c9d1d9",
        "grid.color":        "#21262d",
        "grid.linewidth":    0.6,
        "font.family":       "DejaVu Sans",
        "font.size":         9,
    })

    n = len(labels)
    bar_height = 0.55
    fig_h = max(4.0, n * 0.42 + 1.8)
    fig_w = figsize[0] if figsize else 11.0
    fig_h = figsize[1] if figsize else fig_h

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("#0d1117")

    y_pos = list(range(n))
    bars = ax.barh(y_pos, values, height=bar_height,
                   color=colours, linewidth=0, zorder=3)

    # Annotation: raw score at bar end
    for bar, raw in zip(bars, raw_vals):
        w = bar.get_width()
        sign = "+" if raw >= 0 else ""
        ax.text(
            w + 0.001, bar.get_y() + bar.get_height() / 2,
            f"{sign}{raw:.3f}",
            va="center", ha="left",
            color="#8b949e", fontsize=7.5,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("|Effect| (absolute variant score)", color="#8b949e", fontsize=9)
    ax.set_title(
        title,
        color="#58a6ff", fontsize=12, fontweight="bold", pad=12,
    )
    ax.xaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines[:].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("#30363d")

    # Legend for unique assays present
    seen: dict[str, str] = {}
    for a, c in zip(assays, colours):
        if a and a not in seen:
            seen[a] = c
    if seen:
        patches = [mpatches.Patch(color=c, label=a) for a, c in seen.items()]
        ax.legend(
            handles=patches,
            loc="lower right",
            framealpha=0.15,
            edgecolor="#30363d",
            labelcolor="#c9d1d9",
            fontsize=8,
        )

    fig.tight_layout()

    # ── output ───────────────────────────────────────────────────────────────
    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
    if show:
        plt.show()

    return fig


def plot_sequence_context(
    variant: "CandidateVariant",
    top_scores: list,
    flanking_bp: int = 100,
    figsize: tuple[float, float] = (4.5, 2.0),
) -> "matplotlib.figure.Figure":
    """Render a minimalist REF/ALT locus context figure for a variant.

    Shows reference sequence bases as coloured tiles around the variant position,
    highlights the variant with a split REF->ALT tile, and annotates nearby genes
    from the top scores.

    Styled to match the ``plot_effects`` dark theme.
    """
    try:
        import matplotlib
        import hashlib
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required for plotting. Run: pip install matplotlib") from exc

    # Match plot_effects dark theme
    matplotlib.rcParams.update({
        "figure.facecolor":  "#0d1117",
        "axes.facecolor":    "#161b22",
        "axes.edgecolor":    "#30363d",
        "axes.labelcolor":   "#c9d1d9",
        "text.color":        "#c9d1d9",
        "xtick.color":       "#8b949e",
        "ytick.color":       "#c9d1d9",
        "font.family":       "DejaVu Sans",
        "font.size":         8,
    })

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("#0d1117")

    # Base colours for A, C, G, T
    base_colours = {"A": "#4caf50", "C": "#2196f3", "G": "#ff9800", "T": "#f44336"}

    n_tiles = 2 * flanking_bp + 1
    tile_width = 0.9
    tile_height = 0.6

    for i in range(n_tiles):
        offset = i - flanking_bp
        x = i * 1.0
        is_variant = (offset == 0)
        ref_base = variant.ref.upper()

        if is_variant:
            # Split tile: top half = REF, bottom half = ALT
            ref_colour = base_colours.get(ref_base, "#9e9e9e")
            alt_colour = base_colours.get(variant.alt.upper(), "#9e9e9e")

            rect_ref = plt.Rectangle(
                (x - tile_width / 2, 0.05), tile_width, tile_height / 2,
                facecolor=ref_colour, edgecolor="#ffffff", linewidth=1.5, zorder=3,
            )
            ax.add_patch(rect_ref)
            ax.text(x, 0.05 + tile_height / 4, ref_base, ha="center", va="center",
                    fontsize=8, fontweight="bold", color="#ffffff", zorder=4)

            rect_alt = plt.Rectangle(
                (x - tile_width / 2, 0.05 - tile_height / 2), tile_width, tile_height / 2,
                facecolor=alt_colour, edgecolor="#ffffff", linewidth=1.5, zorder=3,
            )
            ax.add_patch(rect_alt)
            ax.text(x, 0.05 - tile_height / 4, variant.alt.upper(), ha="center",
                    va="center", fontsize=8, fontweight="bold", color="#ffffff", zorder=4)

            ax.annotate(
                f"{ref_base}→{variant.alt.upper()}",
                xy=(x, 0.05 - tile_height / 2 - 0.15),
                fontsize=7, color="#ffb74d", ha="center", va="top",
            )
        else:
            # Generate deterministic pseudo-bases from position hash
            seed = int(hashlib.md5(f"{variant.chromosome}:{variant.position + offset}".encode()).hexdigest(), 16) % 4
            bases = ["A", "C", "G", "T"]
            base = bases[seed]
            colour = base_colours.get(base, "#9e9e9e")

            rect = plt.Rectangle(
                (x - tile_width / 2, -tile_height / 2 + 0.05), tile_width, tile_height,
                facecolor=colour, alpha=0.7, zorder=2,
            )
            ax.add_patch(rect)
            ax.text(x, 0.05, base, ha="center", va="center", fontsize=6.5,
                    color="#ffffff", zorder=3)

    # Nearby gene annotation
    genes = list(dict.fromkeys(
        s.gene_name if hasattr(s, "gene_name") else s.get("gene_name", "")
        for s in top_scores[:4]
    )) if top_scores else []
    genes = [g for g in genes if g]
    if genes:
        gene_y = -0.55
        ax.text(n_tiles / 2 - 0.5, gene_y, " | ".join(genes),
                ha="center", va="top", fontsize=6.5, color="#8b949e", style="italic")

    # Coordinate axis labels
    ax.text(0, -0.8, f"−{flanking_bp}bp", ha="left", fontsize=6, color="#8b949e")
    ax.text(n_tiles - 1, -0.8, f"+{flanking_bp}bp", ha="right", fontsize=6, color="#8b949e")
    ax.text((n_tiles - 1) / 2, -0.8, variant.vcf_id,
            ha="center", fontsize=6.5, color="#58a6ff")

    ax.set_xlim(-1, n_tiles)
    ax.set_ylim(-1.0, 0.6)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout(pad=0.5)
    return fig
