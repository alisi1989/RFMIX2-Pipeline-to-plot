#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Dict, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd


LOGGER = logging.getLogger("GAP_Plot")

DEFAULT_COLORS = [
    "#a32e2e", "#0a0ae0", "#bfa004", "#d18311", "#22ba9d",
    "#839dfc", "#9a5dc1", "#26962b", "#707070", "#00cfff", "#790ee0",
    "#ff4d6d", "#2d6a4f", "#f77f00", "#4ea8de",
]

VALUE_RE = re.compile(r"^\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)")
COLOR_RE = re.compile(r"(#[0-9a-fA-F]{6}|#[0-9a-fA-F]{3})")


def parse_gap_table(df: pd.DataFrame) -> Tuple[pd.Series, pd.DataFrame, Dict[str, str]]:
    if df.shape[1] < 2:
        raise ValueError("No ancestry columns found in GAP file.")

    individuals = df.iloc[:, 0].astype(str)
    ancestries = list(df.columns[1:])
    values = pd.DataFrame(index=df.index)
    colors: Dict[str, str] = {}

    for i, ancestry in enumerate(ancestries):
        col = df[ancestry].astype(str)
        parsed = col.str.extract(VALUE_RE)[0]
        if parsed.isna().any():
            missing = ", ".join(individuals[parsed.isna()].head(5).tolist())
            raise ValueError(f"Could not parse numeric ancestry values for {ancestry}: {missing}")
        values[ancestry] = parsed.astype(float)

        color_match = col.str.extract(COLOR_RE)[0].dropna()
        colors[ancestry] = color_match.iloc[0] if not color_match.empty else DEFAULT_COLORS[i % len(DEFAULT_COLORS)]

    return individuals, values, colors


def plot_admixture(
    df: pd.DataFrame,
    output_file: str,
    title: str,
    normalize: bool,
    width: float | None,
    height: float,
) -> None:
    individuals, proportions_df, ancestry_colors = parse_gap_table(df)

    row_sums = proportions_df.sum(axis=1)
    off_sum = (row_sums - 1.0).abs() > 0.01
    if off_sum.any():
        LOGGER.warning(
            "%d individual(s) have ancestry sums outside 1 +/- 0.01. Range: %.4f-%.4f",
            int(off_sum.sum()),
            float(row_sums.min()),
            float(row_sums.max()),
        )

    if normalize:
        safe_sums = row_sums.replace(0, pd.NA)
        proportions_df = proportions_df.div(safe_sums, axis=0).fillna(0.0)

    n_individuals = len(individuals)
    fig_width = width if width is not None else max(8.0, min(42.0, 3.5 + n_individuals * 0.34))
    fig, ax = plt.subplots(figsize=(fig_width, height))

    bottom = [0.0] * n_individuals
    for ancestry in proportions_df.columns:
        proportion = proportions_df[ancestry].to_numpy()
        ax.bar(
            individuals,
            proportion,
            bottom=bottom,
            color=ancestry_colors[ancestry],
            label=ancestry,
            width=0.96,
            linewidth=0,
        )
        bottom = [x + y for x, y in zip(bottom, proportion)]

    ax.set_ylabel("Global Ancestry Proportion")
    ax.set_ylim(0, 1.0)
    ax.set_title(title)
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1), frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fontsize = 8 if n_individuals <= 40 else 6
    plt.xticks(rotation=90, fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(output_file, format=Path(output_file).suffix.lstrip(".") or "pdf", bbox_inches="tight")
    plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create an admixture chart from a GAP file.")
    parser.add_argument("--input", help="Input .gap file", required=True)
    parser.add_argument("--output", help="Output graph file, usually .pdf", required=True)
    parser.add_argument("--title", default="Global Ancestry Admixture", help="Plot title")
    parser.add_argument("--normalize", action="store_true", help="Normalize each individual's ancestry values to sum to 1 before plotting.")
    parser.add_argument("--width", type=float, help="Figure width in inches. Default scales with sample count.")
    parser.add_argument("--height", type=float, default=6.0, help="Figure height in inches.")
    parser.add_argument("--quiet", action="store_true", help="Only warnings and errors.")
    return parser.parse_args(argv)


def setup_logging(args: argparse.Namespace) -> None:
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO, format="[%(levelname)s] %(message)s")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    setup_logging(args)

    df = pd.read_csv(args.input, sep="\t")
    plot_admixture(df, args.output, args.title, args.normalize, args.width, args.height)
    LOGGER.info("Wrote %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
