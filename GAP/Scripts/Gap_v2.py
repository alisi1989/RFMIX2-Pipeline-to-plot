#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Example of usage: ###
# Basic run (combine chr1-22 into global ancestry)
# python GAP.py --prefix results/sample_ --output-dir out --final-output out/global_ancestry.gap

# Sort individuals by African ancestry
# python GAP.py --prefix results/sample_ --output-dir out --final-output out/global_ancestry.gap --sort-ancestry AFR

# Debug mode
# python GAP.py --prefix results/sample_ --output-dir out --final-output out/global_ancestry.gap --debug

"""
Global Ancestry Painting (GAP)

This script processes RFMix v2 global ancestry output files (`.rfmix.Q`),
combines per-chromosome data for each individual, calculates global ancestry
proportions, and generates a final `.gap` file with optional color annotations.

Features:
---------
1. Combines multiple chromosomes into a single global ancestry profile per individual.
2. Computes mean ancestry proportions across chromosomes.
3. Allows optional sorting of individuals by a chosen ancestry column.
4. Annotates ancestry values with colors (defaults provided).
5. Produces a clean `.gap` file ready for downstream visualization.

Authors: Alessandro Lisi and Michael C. Campbell
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

# Default color palette (supports ancestry0 -> ancestry14)
DEFAULT_COLORS = [
    "#a32e2e", "#0a0ae0", "#bfa004", "#d18311", "#22ba9d",
    "#839dfc", "#9a5dc1", "#26962b", "#707070", "#00cfff", "#790ee0",
    "#ff4d6d", "#2d6a4f", "#f77f00", "#4ea8de",
]
import re

LOGGER = logging.getLogger("GAP")

CHROMOSOME_LENGTHS_HG38 = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
    "M": 16569,
    "MT": 16569,
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def ensure_dir(path: Path) -> None:
    """Ensure the output directory exists."""
    path.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def _normalize_chr_token(ch: str) -> str:
    ch = str(ch)
    ch = ch.replace("chr", "").replace("CHR", "")
    ch = ch.upper()
    return "MT" if ch == "M" else ch


def _guess_chrom_from_filename(name: str) -> str:
    """Extract chrom token from filenames like '*chr22.rfmix.Q'."""
    m = re.findall(r"chr([A-Za-z0-9]+)", name, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Could not parse chromosome from filename: {name}")
    return m[-1]


def discover_files(prefix: str, chroms: Sequence[str]) -> List[Path]:
    """Discover `.rfmix.Q` files (one per chromosome) for the given prefix.

    New behavior: RFMix2 Q files typically contain ALL samples per chromosome.
    We therefore return a list of chromosome Q files to combine, rather than grouping by individual.
    """
    prefix_path = Path(prefix)
    if prefix_path.is_file():
        if prefix_path.name.endswith(".rfmix.Q"):
            chrom = _guess_chrom_from_filename(prefix_path.name)
            req = {_normalize_chr_token(c).upper() for c in chroms}
            if _normalize_chr_token(chrom).upper() in req:
                return [prefix_path]
        LOGGER.warning("Prefix points to a file that is not a requested .rfmix.Q: %s", prefix_path)
        return []

    search_dir = prefix_path.parent if str(prefix_path.parent) not in ("", ".") else Path(".")
    prefix_base = prefix_path.name

    req = {_normalize_chr_token(c).upper() for c in chroms}

    hits = sorted(search_dir.glob(f"{prefix_base}*chr*.rfmix.Q"))
    out: List[Path] = []
    for fp in hits:
        if not fp.is_file():
            continue
        chrom = _guess_chrom_from_filename(fp.name)
        if _normalize_chr_token(chrom).upper() not in req:
            continue
        out.append(fp)

    # sort by chromosome numeric order if possible
    def chrom_key(p: Path) -> Tuple[int, str]:
        c = _normalize_chr_token(_guess_chrom_from_filename(p.name)).upper()
        if c.isdigit():
            return (int(c), c)
        # X/Y/M at the end
        order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        return (order.get(c, 99), c)

    out = sorted(out, key=chrom_key)
    return out



def read_q_file(fp: Path) -> Tuple[List[str], List[str], pd.DataFrame]:
    """Read a single RFMix `.rfmix.Q` file.

    Expected format:
    - first two lines are headers
    - subsequent lines: sample_id + ancestry proportions
    Returns: (header_lines, ancestry_cols, df)
    """
    with fp.open("r", encoding="utf-8") as fh:
        h1 = next(fh).strip()
        h2 = next(fh).strip()

    # Second header line contains ancestry column names (excluding sample)
    parts = h2.split("\t")
    ancestry_cols = parts[1:]

    df = pd.read_csv(
        fp,
        sep="\t",
        skiprows=2,
        names=["sample"] + ancestry_cols,
        dtype={"sample": str, **{c: float for c in ancestry_cols}},
    )
    return [h1, h2], ancestry_cols, df


def weight_for_chromosome(chrom: str, weighting: str) -> float:
    chrom_norm = _normalize_chr_token(chrom)
    if weighting == "equal":
        return 1.0
    if weighting == "chrom-length":
        try:
            return float(CHROMOSOME_LENGTHS_HG38[chrom_norm])
        except KeyError as exc:
            raise ValueError(f"No hg38 chromosome length configured for chromosome '{chrom}'") from exc
    raise ValueError(f"Unsupported weighting mode: {weighting}")


def calculate_mean(
    q_files: Sequence[Path],
    requested_chroms: Sequence[str],
    weighting: str = "chrom-length",
    require_all_chroms: bool = False,
) -> Tuple[List[str], str, pd.DataFrame]:
    """Calculate per-sample global ancestry proportions across Q files.

    Works for both common layouts:
    - one multi-sample `.rfmix.Q` per chromosome
    - one single-sample `.rfmix.Q` per sample/chromosome

    Each sample is divided by the chromosomes actually observed for that sample.
    When `require_all_chroms` is set, samples missing any requested chromosome are
    dropped with a warning.

    Returns:
        (header_lines, header_line_string, df_mean)
    where df_mean has columns: sample + ancestry columns.
    """
    if not q_files:
        raise ValueError("No .rfmix.Q files provided")

    header_lines: List[str] = []
    ancestry_cols_ref: List[str] = []
    sums: Dict[str, np.ndarray] = {}
    weights: Dict[str, float] = {}
    observed: Dict[str, set[str]] = {}
    duplicate_sample_chroms: set[Tuple[str, str]] = set()

    for i, fp in enumerate(q_files):
        headers, ancestry_cols, df = read_q_file(fp)
        chrom = _normalize_chr_token(_guess_chrom_from_filename(fp.name))
        weight = weight_for_chromosome(chrom, weighting)

        if i == 0:
            header_lines = headers
            ancestry_cols_ref = ancestry_cols
        else:
            if ancestry_cols != ancestry_cols_ref:
                raise ValueError(
                    f"Ancestry columns differ across files. First: {ancestry_cols_ref}, {fp.name}: {ancestry_cols}"
                )

        for _, row in df.iterrows():
            sample = str(row["sample"])
            values = row[ancestry_cols_ref].astype(float).to_numpy()
            key = (sample, chrom)
            if chrom in observed.get(sample, set()):
                duplicate_sample_chroms.add(key)
            observed.setdefault(sample, set()).add(chrom)
            sums[sample] = sums.get(sample, np.zeros(len(ancestry_cols_ref), dtype=float)) + values * weight
            weights[sample] = weights.get(sample, 0.0) + weight

    if duplicate_sample_chroms:
        preview = ", ".join(f"{sample}/chr{chrom}" for sample, chrom in sorted(duplicate_sample_chroms)[:10])
        LOGGER.warning("Duplicate sample/chromosome rows were combined additively: %s", preview)

    req = {_normalize_chr_token(c) for c in requested_chroms}
    records: List[Dict[str, object]] = []
    skipped: List[str] = []
    missing_summaries: Dict[Tuple[str, ...], int] = {}

    for sample in sorted(sums):
        missing = tuple(sorted(req - observed.get(sample, set()), key=lambda c: (0, int(c)) if c.isdigit() else (1, c)))
        if missing:
            missing_summaries[missing] = missing_summaries.get(missing, 0) + 1
            if require_all_chroms:
                skipped.append(sample)
                continue
        row_values = sums[sample] / weights[sample]
        records.append({"sample": sample, **dict(zip(ancestry_cols_ref, row_values))})

    if missing_summaries:
        for missing, count in sorted(missing_summaries.items(), key=lambda kv: (-kv[1], kv[0]))[:8]:
            LOGGER.warning("%d sample(s) missing requested chromosome(s): %s", count, ", ".join(missing))
    if skipped:
        LOGGER.warning("Skipped %d sample(s) because --require-all-chroms was set.", len(skipped))

    mean_df = pd.DataFrame.from_records(records, columns=["sample"] + ancestry_cols_ref)
    header_line = header_lines[1]  # the tab header with ancestry names
    return header_lines, header_line, mean_df



def combine_to_bed(
    prefix: str,
    chroms: Sequence[str],
    out_dir: Path,
    sort_ancestry: str = None,
    weighting: str = "chrom-length",
    require_all_chroms: bool = False,
) -> Path:
    """Combine RFMix `.rfmix.Q` files (multi-sample) into a single bed-like file with mean values."""
    q_files = discover_files(prefix, chroms)
    if not q_files:
        raise FileNotFoundError(f"No .rfmix.Q files found for prefix '{prefix}' and chromosomes {list(chroms)}")
    LOGGER.info("Discovered %d .rfmix.Q file(s). Weighting mode: %s", len(q_files), weighting)

    headers, header_line, df_mean = calculate_mean(q_files, chroms, weighting=weighting, require_all_chroms=require_all_chroms)

    # Optional sorting by a chosen ancestry column
    if sort_ancestry:
        ancestry_cols = header_line.split("\t")[1:]
        if sort_ancestry in ancestry_cols:
            df_mean = df_mean.sort_values(by=sort_ancestry, ascending=False)
        else:
            LOGGER.warning("Ancestry '%s' not found; sorting disabled.", sort_ancestry)

    out_path = out_dir / (Path(prefix).name.strip("_") + ".bed")
    with out_path.open("w", encoding="utf-8") as out:
        out.write(headers[0] + "\n")
        out.write(header_line + "\n")
        df_mean.to_csv(out, sep="\t", index=False, header=False, float_format="%.5f")

    LOGGER.info("Intermediate BED file created: %s", out_path)
    return out_path


def add_colors(df: pd.DataFrame, ancestry_colors: Dict[str, str]) -> pd.DataFrame:
    """
    Annotate ancestry values with colors.

    Example: 0.12345 -> "0.12345 (#a32e2e)"
    """
    for anc, col in ancestry_colors.items():
        if anc in df.columns:
            df[anc] = df[anc].fillna(0.0).apply(lambda x: f"{x:.5f} ({col})")
    return df


def load_color_config(path: str | None) -> Dict[str, str]:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    return {str(k): str(v) for k, v in cfg.items()}


def process_with_colors(in_bed: Path, out_gap: Path, colors: List[str], color_config: str | None = None) -> None:
    """
    Apply colors to the intermediate `.bed` file and produce the final `.gap`.
    Output: file con UNA sola riga di header (sample + ancestry columns).
    """
    # Salta la prima riga di commento e leggi solo l'header utile
    with in_bed.open("r") as f:
        _ = f.readline()  # <-- ignoriamo h1 (#sample)
        header = f.readline().strip()

    ancestry_cols = header.split("\t")[1:]

    # Se non hai abbastanza colori, li ripeti ciclicamente
    if len(colors) < len(ancestry_cols):
        times = len(ancestry_cols) // len(colors) + 1
        colors = (colors * times)[:len(ancestry_cols)]

    ancestry_colors = dict(zip(ancestry_cols, colors))
    ancestry_colors.update(load_color_config(color_config))

    # Leggi il dataframe saltando le prime due righe
    df = pd.read_csv(in_bed, sep="\t", skiprows=2, names=["sample"] + ancestry_cols)
    df = add_colors(df, ancestry_colors)

    # Scrivi solo l'header corretto + i dati
    with out_gap.open("w", encoding="utf-8") as out:
        out.write(header + "\n")
        df.to_csv(out, sep="\t", index=False, header=False)
    LOGGER.info("Final GAP file created: %s", out_gap)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Global Ancestry Painting (GAP)")
    parser.add_argument("--prefix", required=True, help="Prefix of RFMix .rfmix.Q files")
    parser.add_argument("--chr", nargs="+", default=[str(i) for i in range(1, 23)],
                        help="Chromosomes to include (default: 1-22)")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--sort-ancestry", help="Ancestry column to sort individuals by (descending)")
    parser.add_argument("--final-output", required=True, help="Final GAP file path (.gap)")
    parser.add_argument("--color-config", help="JSON file mapping ancestry column names to colors.")
    parser.add_argument(
        "--weighting",
        choices=["chrom-length", "equal"],
        default="chrom-length",
        help="How to combine chromosome-level Q values into genome-wide ancestry.",
    )
    parser.add_argument(
        "--require-all-chroms",
        action="store_true",
        help="Drop samples that are missing any requested chromosome.",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--debug", action="store_true", help="Enable debug logging")
    group.add_argument("--quiet", action="store_true", help="Suppress info logging")
    return parser.parse_args()


def setup_logging(args: argparse.Namespace) -> None:
    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")


def main() -> int:
    args = parse_args()
    setup_logging(args)

    out_dir = Path(args.output_dir)
    ensure_dir(out_dir)

    bed_file = combine_to_bed(
        args.prefix,
        args.chr,
        out_dir,
        args.sort_ancestry,
        weighting=args.weighting,
        require_all_chroms=args.require_all_chroms,
    )

    final_out = Path(args.final_output)
    if final_out.suffix != ".gap":
        final_out = final_out.with_suffix(".gap")

    process_with_colors(bed_file, final_out, DEFAULT_COLORS, color_config=args.color_config)

    try:
        bed_file.unlink()
        LOGGER.debug("Removed intermediate file: %s", bed_file)
    except Exception as e:
        LOGGER.warning("Could not remove intermediate file %s: %s", bed_file, e)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
