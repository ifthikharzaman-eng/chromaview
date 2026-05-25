"""Publication-quality chromatogram image export.

Generates PNG or SVG chromatogram images using matplotlib for
headless/batch export. For interactive use, the PyQtGraph widget's
own export is preferred (higher fidelity to what the user sees).
"""
from __future__ import annotations

from pathlib import Path
from typing import Union, Optional

import numpy as np

from ..core.models import Chromatogram, BASE_COLORS


def export_image(
    chrom: Chromatogram,
    filepath: Union[str, Path],
    start_base: int = 0,
    end_base: Optional[int] = None,
    width: int = 16,
    height: int = 4,
    dpi: int = 300,
    show_bases: bool = True,
    show_quality: bool = False,
    title: Optional[str] = None,
) -> Path:
    """Export a chromatogram region as a publication-quality image.

    Uses matplotlib for rendering (works in headless environments).

    Args:
        chrom: Source chromatogram.
        filepath: Output path (.png or .svg).
        start_base: First base index to show (0-based).
        end_base: Last base index (exclusive). None = all bases.
        width, height: Figure size in inches.
        dpi: Resolution (for PNG).
        show_bases: Annotate peaks with base letters.
        show_quality: Add quality bar chart below trace.
        title: Optional figure title.

    Returns:
        Path to the written file.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
    except ImportError:
        raise ImportError(
            "matplotlib is required for image export. "
            "Install with: pip install matplotlib"
        )

    filepath = Path(filepath)
    if end_base is None:
        end_base = chrom.num_bases

    # Determine sample range from peak locations
    if len(chrom.peak_locations) == 0:
        return filepath

    sample_start = int(chrom.peak_locations[max(0, start_base)] - 50) if start_base < len(chrom.peak_locations) else 0
    sample_end = int(chrom.peak_locations[min(end_base - 1, len(chrom.peak_locations) - 1)] + 50)
    sample_start = max(0, sample_start)

    # ── Create figure ────────────────────────────────────────────
    n_axes = 2 if show_quality else 1
    fig, axes = plt.subplots(
        n_axes, 1,
        figsize=(width, height * n_axes * 0.7),
        gridspec_kw={"height_ratios": [3, 1] if show_quality else [1]},
        squeeze=False,
    )
    ax_trace = axes[0, 0]

    # ── Plot traces ──────────────────────────────────────────────
    for base, color in BASE_COLORS.items():
        trace = chrom.traces.get(base)
        if trace is None:
            continue
        x = np.arange(len(trace))
        mask = (x >= sample_start) & (x <= sample_end)
        ax_trace.plot(x[mask], trace[mask], color=color, linewidth=0.6, label=base, alpha=0.85)

    # ── Annotate bases ───────────────────────────────────────────
    if show_bases:
        y_max = ax_trace.get_ylim()[1]
        for i in range(start_base, min(end_base, len(chrom.basecalls))):
            loc = int(chrom.peak_locations[i])
            base = chrom.basecalls[i]
            color = BASE_COLORS.get(base.upper(), "#666666")
            ax_trace.text(
                loc, -y_max * 0.06, base,
                ha="center", va="top", fontsize=6, fontweight="bold",
                color=color, family="monospace",
            )

    ax_trace.set_xlim(sample_start, sample_end)
    ax_trace.set_ylabel("Fluorescence", fontsize=9)
    ax_trace.legend(loc="upper right", fontsize=7, framealpha=0.7)
    ax_trace.spines["top"].set_visible(False)
    ax_trace.spines["right"].set_visible(False)

    if title:
        ax_trace.set_title(title, fontsize=11, fontweight="bold")

    # ── Quality bars ─────────────────────────────────────────────
    if show_quality and n_axes > 1:
        ax_qual = axes[1, 0]
        positions = chrom.peak_locations[start_base:end_base]
        qualities = chrom.quality_scores[start_base:end_base]

        colors = []
        for q in qualities:
            if q >= 30:
                colors.append("#2ecc71")
            elif q >= 20:
                colors.append("#f39c12")
            else:
                colors.append("#e74c3c")

        ax_qual.bar(positions, qualities, width=4, color=colors, alpha=0.8)
        ax_qual.set_xlim(sample_start, sample_end)
        ax_qual.set_ylabel("Phred", fontsize=9)
        ax_qual.set_xlabel("Position", fontsize=9)
        ax_qual.axhline(y=20, color="#aaa", linestyle="--", linewidth=0.5)
        ax_qual.spines["top"].set_visible(False)
        ax_qual.spines["right"].set_visible(False)
        ax_qual.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(str(filepath), dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return filepath
