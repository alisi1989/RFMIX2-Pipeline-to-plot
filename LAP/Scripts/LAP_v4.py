#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###    Example of usage: ###
# Caso base (tutti i chr1..22)
# python rfmix2_plot_pipeline_pro.py \
# --prefix /path/to/out/prefix_ \
# --output-dir out_bed

#  Con threads e feature su chr2
#  python rfmix2_plot_pipeline_pro.py \
#  --prefix results/sample_ \
#  --chr chr1 chr2 chr3 \
#  --output-dir out_bed \
#  --threads 4 \
#  --from-bp 135000000 --to-bp 135200000 --chromosome chr2

# Solo lista operazioni, senza scrivere file
# python rfmix2_plot_pipeline_pro.py \
# --prefix results/sample_ --output-dir out_bed --dry-run

"""
RFMix2 → unified ancestry BED generator (single-file, professionalized)

What it does (end‑to‑end):
1) Discovers and combines per‑chromosome RFMix2 .msp.tsv chunks into a single MSP per individual.
2) Converts MSP into haplotype BEDs (hap1 / hap2) with ancestry labels.
3) Merges hap BEDs into a final BED with rendering hints (geom_rect + colors),
   plus optional feature highlight lines (geom_line).

Key improvements over the original:
- Robust file discovery and individual ID parsing (regex‑based, supports chr and no‑chr tokens).
- Safer I/O (streaming copy for huge files, avoids readlines() on big MSPs).
- Pandas CSV parsing with comment handling and typed columns.
- Natural sorting by chromosome and position.
- Optional parallel processing across individuals (--threads).
- Proper logging (--debug / --quiet) instead of prints.
- Options: --require-all-chroms, --keep-temp, --dry-run, color config file (JSON).
- Resilient ancestry mapping: default to "ancestry{int}", or load from JSON mapping.

Supports up to ancestry0..14 (15 total ancestries) by default.

Dependencies: Python ≥3.8, pandas

Author: Alessandro Lisi
"""



from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
LOGGER = logging.getLogger("rfmix2")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
CHROM_TOKEN_RE = re.compile(r"^(?:chr)?(\w+)$", re.IGNORECASE)


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def strip_chr(token: str) -> str:
    """Remove 'chr' prefix if present (case‑insensitive)."""
    m = CHROM_TOKEN_RE.match(str(token))
    return m.group(1) if m else str(token)


def chrom_sort_key(chrom: str) -> Tuple[int, str]:
    """Natural chromosome sort: 1..22,X,Y,MT (case insensitive)."""
    c = strip_chr(chrom).upper()
    if c in {"X", "XX"}:
        return (1000, "X")
    if c in {"Y"}:
        return (1001, "Y")
    if c in {"MT", "M"}:
        return (1002, "MT")
    try:
        return (int(c), "")
    except ValueError:
        # fallback: non‑numeric contigs at the end, but keep stable order
        return (2000, c)


@dataclass
class Feature:
    chrom: str
    start_bp: int
    end_bp: int
    label: str = "manual_feature"


@dataclass(frozen=True)
class FlareSnpAnnotation:
    chrom: str
    pos: int
    variant_id: str
    snp_label: str
    ref: str
    alt: str
    vcf_sample: str
    gt: str
    an1_code: Optional[int]
    an1_label: str
    an2_code: Optional[int]
    an2_label: str
    anp1: Optional[str] = None
    anp2: Optional[str] = None


@dataclass(frozen=True)
class FlareAnnotationSet:
    by_sample: Dict[str, List[FlareSnpAnnotation]]
    ancestry_map: Dict[int, str]
    matched_variant_count: int
    requested_variant_count: int


# -----------------------------------------------------------------------------
# Discovery & Combination
# -----------------------------------------------------------------------------
MSP_SUFFIX = ".msp.tsv"

DEFAULT_ANCESTRY_COLORS = [
    "#a32e2e",
    "#0a0ae0",
    "#bfa004",
    "#d18311",
    "#22ba9d",
    "#839dfc",
    "#9a5dc1",
    "#26962b",
    "#707070",
    "#00cfff",
    "#790ee0",
    "#ff4d6d",
    "#2d6a4f",
    "#f77f00",
    "#4ea8de",
]




def _extract_individual_from_basename(basename: str, chrom: str) -> str:
    """Infer the 'individual' stem from a basename like '<prefix>_<sample>_<chr>.msp.tsv'.

    Strategy: remove the trailing '<sep>chrom' token and the suffix; then rstrip separators.
    This preserves the original prefix+sample stem, which is what the original script used
    as the individual key.
    """
    name = basename
    if name.endswith(MSP_SUFFIX):
        name = name[: -len(MSP_SUFFIX)]
    # Remove a single trailing chrom token preceded by optional separators
    # e.g., 'foo_bar_chr1' -> 'foo_bar', 'foo1chr1' -> 'foo1'
    chrom_escaped = re.escape(chrom)
    name = re.sub(rf"[._-]?{chrom_escaped}$", "", name)
    return name.rstrip("._-")


def discover_inputs(prefix: str, chroms: Sequence[str]) -> Tuple[Dict[str, List[Path]], List[str]]:
    """Discover .msp.tsv files grouped by individual stem.

    This implementation avoids substring/glob pitfalls (e.g., searching for chrom '1' accidentally
    matching files ending in '21'). We instead:
    - glob all MSPs for the prefix once
    - parse the chromosome token from each filename
    - filter by the requested chromosome set (normalized via strip_chr)

    If prefix is an existing .msp.tsv file, treat it as a pre-combined MSP.
    This supports files such as sample_chr1-22.msp.tsv that already contain all
    requested chromosomes in one table.
    """
    individuals: Dict[str, List[Path]] = {}

    prefix_path = Path(prefix)

    if prefix_path.is_file():
        if prefix_path.name.endswith(MSP_SUFFIX):
            individual = prefix_path.name[: -len(MSP_SUFFIX)]
            return {individual: [prefix_path]}, [individual]
        LOGGER.warning("Prefix points to a file that is not an MSP TSV: %s", prefix_path)
        return {}, []

    search_dir = prefix_path.parent if str(prefix_path.parent) not in ("", ".") else Path(".")
    prefix_base = prefix_path.name

    # Normalize requested chrom set (e.g., {'1','2',...,'22','X'})
    req_norm = {strip_chr(c).upper() for c in chroms}

    # One pass: collect all MSP files matching the prefix
    all_hits = sorted(search_dir.glob(f"{prefix_base}*{MSP_SUFFIX}"))
    if not all_hits:
        return {}, []

    for fp in all_hits:
        if not fp.is_file():
            continue
        chrom_tok = _guess_chrom_from_filename(fp.name)
        chrom_norm = strip_chr(chrom_tok).upper()
        if chrom_norm not in req_norm:
            continue

        ind = _extract_individual_from_basename(fp.name, chrom_tok)
        individuals.setdefault(ind, []).append(fp)

    # Sort each individual's files by chrom natural order then by name
    for ind, files in individuals.items():
        individuals[ind] = sorted(files, key=lambda p: (chrom_sort_key(_guess_chrom_from_filename(p.name)), p.name))

    unique_inds = sorted(individuals.keys())
    return individuals, unique_inds
def filter_headers_for_bed(headers: List[str]) -> List[str]:
    """Keep only informative comment headers for BED outputs.

    We drop the '#chm\t...' header line because in multi-sample MSPs it includes all sample columns
    and is misleading in per-sample BED outputs.
    """
    out: List[str] = []
    for h in headers:
        if h.startswith("#chm"):
            continue
        out.append(h)
    return out


def _guess_chrom_from_filename(basename: str) -> str:
    """Best effort to recover the chrom token from a filename that ends with '<chrom>.msp.tsv'."""
    core = basename
    if core.endswith(MSP_SUFFIX):
        core = core[: -len(MSP_SUFFIX)]
    # Grab last token after separators and hope it's the chrom (works for common RFMix2 names)
    last = re.split(r"[._-]", core)[-1]
    # normalize to chrX if it's purely numeric or already has chr
    if CHROM_TOKEN_RE.match(last):
        return last if last.lower().startswith("chr") else f"chr{last}"
    return last


# -----------------------------------------------------------------------------
# MSP → hap BEDs
# -----------------------------------------------------------------------------
BASE_MSP_COLS = [
    "chm", "spos", "epos", "sgpos", "egpos", "n_snps",
]


def parse_msp_columns_from_header(path: Path) -> List[str]:
    """Return column names from the '#chm\t...' header line.

    RFMix2 MSPs often include the column header as a commented line beginning with '#chm'.
    We normalize spaces to underscores (e.g., 'n snps' -> 'n_snps').
    """
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#chm"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                return [c.strip().replace(" ", "_") for c in cols]
    # Fallback: legacy single-individual format (no explicit header line)
    return BASE_MSP_COLS + ["ind1", "ind2"]


_HAP_COL_RE = re.compile(r"^(?P<sample>.+)\.(?P<hap>[01])$")


def infer_sample_hap_pairs(columns: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """Infer diploid sample -> (hap0_col, hap1_col) from columns.

    For multi-sample MSPs, hap columns are typically named like '<sample>.0' and '<sample>.1'.
    """
    pairs: Dict[str, Dict[str, str]] = {}
    for c in columns:
        m = _HAP_COL_RE.match(str(c))
        if not m:
            continue
        sample = m.group("sample")
        hap = m.group("hap")
        pairs.setdefault(sample, {})[hap] = str(c)

    out: Dict[str, Tuple[str, str]] = {}
    for sample, d in pairs.items():
        if "0" in d and "1" in d:
            out[sample] = (d["0"], d["1"])
    return out


def safe_sample_name(sample: str) -> str:
    """Make a sample name safe for filesystem paths."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", sample)


def read_msp(path: Path, columns: Sequence[str]) -> pd.DataFrame:
    """Load an MSP TSV with comment headers ('#') into a typed DataFrame.

    Supports both legacy single-individual MSPs (ind1/ind2) and multi-sample MSPs with many hap columns.
    """
    cols = list(columns)

    # Build dtype map: first 6 columns typed, remaining columns as nullable integers (ancestry codes).
    # Some MSP exports can contain a final tab-only line; nullable dtypes let us drop it after parsing.
    dtype_map: Dict[str, object] = {
        "chm": str,
        "spos": float,
        "epos": float,
        "sgpos": float,
        "egpos": float,
        "n_snps": "Int64",
    }
    for c in cols[6:]:
        # ancestry codes; allow missing
        dtype_map[str(c)] = "Int16"

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        dtype=dtype_map,
        engine="c",
    )

    df = df.dropna(how="all")
    df = df[df["chm"].notna() & df["spos"].notna() & df["epos"].notna()].copy()
    df = df[df["chm"].astype(str).str.strip().ne("")]

    # Cast positions to integers if they are integral floats
    for col in ("spos", "epos"):
        df[col] = df[col].astype(float).round().astype(int)
    return df


def load_headers_from_msp(path: Path) -> List[str]:
    headers: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
            else:
                break
    return headers


def open_text_auto(path: Path):
    """Open plain text or gzip-compressed text using the file suffix."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def parse_name_code_tokens(payload: str) -> Dict[int, str]:
    """Parse tokens like 'Europeans=0\tEast_Asia=1' into {0: 'Europeans'}.

    Both RFMix/MSP '#Subpopulation order/codes:' and FLARE '##ANCESTRY=<...>'
    use name=integer payloads; separators may be tabs, spaces, or commas.
    """
    out: Dict[int, str] = {}
    payload = payload.strip().strip("<>")
    for item in re.split(r"[,\t ]+", payload):
        item = item.strip()
        if not item or "=" not in item:
            continue
        name, code = item.rsplit("=", 1)
        name = name.strip()
        try:
            out[int(code)] = name
        except ValueError:
            LOGGER.warning("Could not parse ancestry code in token: %s", item)
    return out


def parse_subpopulation_codes(headers: Sequence[str]) -> Dict[int, str]:
    """Parse the RFMix MSP '#Subpopulation order/codes:' header."""
    for header in headers:
        if not header.startswith("#Subpopulation order/codes:"):
            continue
        payload = header.split(":", 1)[1].strip()
        return parse_name_code_tokens(payload)
    return {}


def build_default_ancestry_map() -> Dict[int, str]:
    return {i: f"ancestry{i}" for i in range(0, 128)}


def build_ancestry_map_from_headers(headers: Sequence[str]) -> Dict[int, str]:
    """Prefer the ancestry names declared by the MSP; fallback to ancestry0..N."""
    parsed = parse_subpopulation_codes(headers)
    if parsed:
        return parsed
    return build_default_ancestry_map()


def parse_flare_ancestry_header(line: str) -> Dict[int, str]:
    """Parse FLARE header line: ##ANCESTRY=<Europeans=0,...>."""
    m = re.match(r"^##ANCESTRY=<(.+)>\s*$", line.strip())
    if not m:
        return {}
    return parse_name_code_tokens(m.group(1))


def normalize_snp_chrom(chrom: str) -> str:
    return remove_chr_prefix(str(chrom))


def safe_bed_field(value: object) -> str:
    return str(value).replace("\t", "_").replace("\n", "_").replace("\r", "_")


def parse_snp_list(path: Optional[Path]) -> Tuple[Set[Tuple[str, int]], Set[str], Dict[Tuple[str, int], str], Dict[str, str]]:
    """Read a flexible SNP list.

    Accepted forms include:
    - chr4 46567563
    - chr4:46567563
    - 4:46567563:G:A (also treated as a VCF ID)
    - rsID or VCF variant ID
    - tabular files with headers containing CHROM/POS and optional ID.

    Returns requested positions/IDs plus stable short labels (SNP1, SNP2, ...)
    assigned from the order of the SNP list.
    """
    labels_by_position: Dict[Tuple[str, int], str] = {}
    labels_by_id: Dict[str, str] = {}
    if path is None:
        return set(), set(), labels_by_position, labels_by_id

    header_index: Optional[Dict[str, int]] = None
    chrom_keys = {"chrom", "chromosome", "chr", "#chrom", "#chromosome"}
    pos_keys = {"pos", "position", "bp", "basepair", "base_pair", "start"}
    id_keys = {"id", "snp", "snpid", "rsid", "variant", "variant_id"}
    request_order = 0

    def register_request(position: Optional[Tuple[str, int]] = None, variant_id: Optional[str] = None) -> None:
        nonlocal request_order
        if position is None and not variant_id:
            return
        existing = labels_by_position.get(position) if position is not None else None
        if existing is None and variant_id:
            existing = labels_by_id.get(variant_id)
        if existing is None:
            request_order += 1
            existing = f"SNP{request_order}"
        if position is not None:
            labels_by_position.setdefault(position, existing)
        if variant_id:
            labels_by_id.setdefault(variant_id, existing)

    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#") and not line.lower().startswith("#chrom"):
                continue
            toks = re.split(r"\s+", line)
            toks_lower = [t.lower() for t in toks]

            if header_index is None and (set(toks_lower) & chrom_keys) and (set(toks_lower) & pos_keys):
                header_index = {name: i for i, name in enumerate(toks_lower)}
                continue

            if header_index is not None:
                chrom_i = next((header_index[k] for k in chrom_keys if k in header_index), None)
                pos_i = next((header_index[k] for k in pos_keys if k in header_index), None)
                id_i = next((header_index[k] for k in id_keys if k in header_index), None)
                position: Optional[Tuple[str, int]] = None
                variant_id = toks[id_i] if id_i is not None and id_i < len(toks) else None
                if chrom_i is not None and pos_i is not None and max(chrom_i, pos_i) < len(toks):
                    try:
                        position = (normalize_snp_chrom(toks[chrom_i]), int(float(toks[pos_i])))
                    except ValueError:
                        LOGGER.warning("Skipping SNP-list line %d with non-numeric POS: %s", line_no, line)
                register_request(position, variant_id)
                continue

            if len(toks) == 1:
                token = toks[0]
                m = re.match(r"^(?:chr)?([0-9XYM]+):(\d+)(?::.*)?$", token, flags=re.IGNORECASE)
                position = None
                if m:
                    position = (normalize_snp_chrom(m.group(1)), int(m.group(2)))
                register_request(position, token)
                continue

            # Common two-column / three-column forms: CHROM POS [ID]
            try:
                position = (normalize_snp_chrom(toks[0]), int(float(toks[1])))
                variant_id = toks[2] if len(toks) >= 3 else None
                register_request(position, variant_id)
                continue
            except ValueError:
                pass

            # Fallback: keep the first token as an ID.
            register_request(variant_id=toks[0])

    return set(labels_by_position), set(labels_by_id), labels_by_position, labels_by_id


def load_sample_map(path: Optional[Path]) -> Dict[str, str]:
    """Load LAP-sample to VCF-sample mapping: LAP_sample<TAB>VCF_sample."""
    out: Dict[str, str] = {}
    if path is None:
        return out
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            toks = re.split(r"\s+", line)
            if len(toks) < 2:
                continue
            if toks[0].lower() in {"lap_sample", "lap", "sample"} and toks[1].lower() in {"vcf_sample", "vcf"}:
                continue
            out[toks[0]] = toks[1]
    return out


def load_highlight_regions(path: Optional[Path]) -> List[Feature]:
    """Load manual highlight regions from a text file.

    Expected columns are: chrom start end. Header lines and comments are allowed.
    Optional extra columns are ignored except the fourth column, which is kept as a label.
    """
    features: List[Feature] = []
    if path is None:
        return features

    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            toks = re.split(r"\s+", line)
            if len(toks) < 3:
                LOGGER.warning("Skipping highlight-region line %d with fewer than 3 columns: %s", line_no, line)
                continue
            if toks[0].lower() in {"chrom", "chromosome", "chr", "#chrom"}:
                continue
            try:
                start = int(float(toks[1]))
                end = int(float(toks[2]))
            except ValueError:
                LOGGER.warning("Skipping highlight-region line %d with non-numeric coordinates: %s", line_no, line)
                continue
            if end < start:
                start, end = end, start
            if end == start:
                LOGGER.warning("Skipping zero-width highlight region in line %d: %s", line_no, line)
                continue
            label = toks[3] if len(toks) >= 4 else f"highlight_region_{len(features) + 1}"
            features.append(Feature(chrom=toks[0], start_bp=start, end_bp=end, label=label))

    LOGGER.info("Loaded %d highlight region(s) from %s", len(features), path)
    return features


def load_flare_annotations(
    vcf_path: Path,
    snp_list_path: Optional[Path],
    target_samples: Optional[Set[str]] = None,
) -> FlareAnnotationSet:
    """Parse a FLARE .anc.vcf(.gz) and collect AN1/AN2 calls for requested SNPs."""
    requested_positions, requested_ids, snp_labels_by_position, snp_labels_by_id = parse_snp_list(snp_list_path)
    annotate_all = not requested_positions and not requested_ids
    if annotate_all:
        LOGGER.warning("No --snp-list records were found; annotating all variants present in the FLARE VCF.")

    ancestry_map: Dict[int, str] = {}
    samples: List[str] = []
    selected_indices: List[Tuple[int, str]] = []
    by_sample: Dict[str, List[FlareSnpAnnotation]] = {}
    matched_variants = 0

    with open_text_auto(vcf_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("##ANCESTRY="):
                ancestry_map.update(parse_flare_ancestry_header(line))
                continue
            if line.startswith("#CHROM"):
                cols = line.split("\t")
                samples = cols[9:]
                if target_samples:
                    wanted = set(target_samples)
                    selected_indices = [(i, sample) for i, sample in enumerate(samples) if sample in wanted]
                else:
                    selected_indices = list(enumerate(samples))
                if not selected_indices:
                    LOGGER.warning("No requested VCF samples were found in %s", vcf_path)
                continue
            if line.startswith("#"):
                continue
            if not samples:
                continue

            parts = line.split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_s, variant_id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            chrom_norm = normalize_snp_chrom(chrom)
            variant_matches = annotate_all or (chrom_norm, pos) in requested_positions or variant_id in requested_ids
            if not variant_matches:
                continue
            if annotate_all:
                snp_label = f"SNP{matched_variants + 1}"
            else:
                candidate_labels = [
                    label
                    for label in (
                        snp_labels_by_position.get((chrom_norm, pos)),
                        snp_labels_by_id.get(variant_id),
                    )
                    if label
                ]
                if candidate_labels:
                    snp_label = min(candidate_labels, key=lambda label: int(label[3:]) if label.startswith("SNP") and label[3:].isdigit() else 10**9)
                else:
                    snp_label = f"SNP{matched_variants + 1}"

            fmt = parts[8].split(":")
            fmt_index = {name: i for i, name in enumerate(fmt)}
            try:
                an1_i = fmt_index["AN1"]
                an2_i = fmt_index["AN2"]
            except ValueError:
                raise ValueError(f"FLARE VCF FORMAT lacks AN1/AN2 at {chrom}:{pos}")
            except KeyError:
                raise ValueError(f"FLARE VCF FORMAT lacks AN1/AN2 at {chrom}:{pos}")

            gt_i = fmt_index.get("GT")
            anp1_i = fmt_index.get("ANP1")
            anp2_i = fmt_index.get("ANP2")

            matched_variants += 1
            for sample_idx, sample_name in selected_indices:
                col_idx = 9 + sample_idx
                if col_idx >= len(parts):
                    continue
                values = parts[col_idx].split(":")
                if max(an1_i, an2_i) >= len(values):
                    continue
                an1_raw, an2_raw = values[an1_i], values[an2_i]
                gt = values[gt_i] if gt_i is not None and gt_i < len(values) else "."
                anp1 = values[anp1_i] if anp1_i is not None and anp1_i < len(values) else None
                anp2 = values[anp2_i] if anp2_i is not None and anp2_i < len(values) else None

                def parse_code(raw_code: str) -> Tuple[Optional[int], str]:
                    if raw_code in {"", ".", "NA", "nan"}:
                        return None, "unknown"
                    try:
                        code = int(raw_code)
                    except ValueError:
                        return None, raw_code
                    return code, ancestry_map.get(code, f"ancestry{code}")

                an1_code, an1_label = parse_code(an1_raw)
                an2_code, an2_label = parse_code(an2_raw)
                by_sample.setdefault(sample_name, []).append(
                    FlareSnpAnnotation(
                        chrom=chrom_norm,
                        pos=pos,
                        variant_id=variant_id,
                        snp_label=snp_label,
                        ref=ref,
                        alt=alt,
                        vcf_sample=sample_name,
                        gt=gt,
                        an1_code=an1_code,
                        an1_label=an1_label,
                        an2_code=an2_code,
                        an2_label=an2_label,
                        anp1=anp1,
                        anp2=anp2,
                    )
                )

    requested_count = len(set(snp_labels_by_position.values()) | set(snp_labels_by_id.values()))
    if requested_count and matched_variants == 0:
        LOGGER.warning("No SNPs from %s matched variants in %s", snp_list_path, vcf_path)
    LOGGER.info(
        "Loaded FLARE annotations: %d matching variant(s), %d sample(s) with annotations.",
        matched_variants,
        len(by_sample),
    )
    return FlareAnnotationSet(by_sample, ancestry_map, matched_variants, requested_count)


def resolve_vcf_sample(
    lap_sample: str,
    flare_annotations: Optional[FlareAnnotationSet],
    sample_map: Dict[str, str],
    flare_sample_override: Optional[str],
) -> Optional[str]:
    if flare_annotations is None:
        return None
    if flare_sample_override:
        return flare_sample_override
    if lap_sample in sample_map:
        return sample_map[lap_sample]
    sample_safe = safe_sample_name(lap_sample)
    if sample_safe in sample_map:
        return sample_map[sample_safe]
    if lap_sample in flare_annotations.by_sample:
        return lap_sample
    if sample_safe in flare_annotations.by_sample:
        return sample_safe
    return None


def infer_lap_sample_names(individuals_map: Dict[str, List[Path]]) -> Set[str]:
    """Infer LAP output sample names from MSP headers without loading the full tables."""
    samples: Set[str] = set()
    for individual, files in individuals_map.items():
        if not files:
            continue
        try:
            cols = parse_msp_columns_from_header(files[0])
        except Exception as exc:
            LOGGER.debug("Could not infer sample names from %s: %s", files[0], exc)
            continue
        pairs = infer_sample_hap_pairs(cols)
        if pairs:
            samples.update(pairs.keys())
        else:
            # Legacy two-haplotype MSPs have no explicit sample name; the individual stem is the best available hint.
            samples.add(individual)
    return samples


def build_lap_color_header(code_to_name: Dict[int, str], colors: Dict[str, str]) -> Optional[str]:
    entries: List[str] = []
    seen: Set[str] = set()
    for code in sorted(code_to_name):
        name = code_to_name[code]
        # When no MSP header is available, the fallback map is ancestry0..ancestry127.
        # Keep the emitted palette compact and compatible with the built-in 15 colors.
        if code > 14 and name == f"ancestry{code}":
            continue
        color = colors.get(name, colors.get(f"ancestry{code}"))
        if not color or name in seen:
            continue
        entries.append(f"{safe_bed_field(name)}={safe_bed_field(color)}")
        seen.add(name)
    if not entries:
        return None
    return "#LAP ancestry/colors:\t" + "\t".join(entries)


def msp_to_hap_beds(
    msp_path: Path,
    out_prefix: Path,
    ancestry_map: Dict[int, str],
    keep_headers: Optional[List[str]] = None,
) -> Dict[str, Tuple[Path, Path]]:
    """Convert an MSP to per-sample hap BEDs.

    Returns:
        Dict[sample] -> (hap1_bed_path, hap2_bed_path)

    Supports:
    - Legacy single-individual MSPs with columns ind1/ind2
    - Multi-sample MSPs with hap columns like '<sample>.0' and '<sample>.1'
    """
    cols = parse_msp_columns_from_header(msp_path)
    df = read_msp(msp_path, cols)

    # Convert ancestry integer codes to labels; if not numeric leave as is
    def map_series(s: pd.Series) -> pd.Series:
        codes = pd.to_numeric(s, errors="coerce")
        labels = codes.map(lambda x: ancestry_map.get(int(x), f"ancestry{int(x)}") if pd.notna(x) else None)
        out = labels.astype(object)
        mask_na = codes.isna()
        out.loc[mask_na] = s.loc[mask_na]
        return out

    results: Dict[str, Tuple[Path, Path]] = {}

    # Multi-sample mode: infer sample pairs from column names
    pairs = infer_sample_hap_pairs(cols)
    if pairs:
        for sample, (c0, c1) in sorted(pairs.items()):
            sample_safe = safe_sample_name(sample)
            # hap1 corresponds to .0, hap2 to .1 (consistent within our toolchain)
            hap_outputs: List[Tuple[str, str]] = [("hap1", c0), ("hap2", c1)]
            created: List[Path] = []

            for hap_name, col in hap_outputs:
                sub = df[["chm", "spos", "epos", col]].copy()
                sub[col] = map_series(sub[col])
                sub.columns = ["chm", "spos", "epos", "ancestry"]

                out_path = out_prefix.with_name(out_prefix.name + f"_{sample_safe}_{hap_name}.bed")
                with out_path.open("w", encoding="utf-8") as out:
                    if keep_headers:
                        out.write("\n".join(keep_headers) + "\n")
                    sub.to_csv(
                        out,
                        sep="\t",
                        header=False,
                        index=False,
                        quoting=csv.QUOTE_NONE,
                        escapechar="\\",
                        lineterminator="\n",
                    )
                LOGGER.info("Wrote %s", out_path)
                created.append(out_path)

            results[sample] = (created[0], created[1])

        return results

    # Legacy mode: expect ind1/ind2
    if "ind1" not in df.columns or "ind2" not in df.columns:
        raise ValueError(f"MSP appears to have no sample hap columns and is missing ind1/ind2: {msp_path}")

    hap_outputs_legacy: List[Tuple[str, str]] = [("hap1", "ind1"), ("hap2", "ind2")]
    created_legacy: List[Path] = []
    for hap_name, col in hap_outputs_legacy:
        sub = df[["chm", "spos", "epos", col]].copy()
        sub[col] = map_series(sub[col])
        sub.columns = ["chm", "spos", "epos", "ancestry"]

        out_path = out_prefix.with_name(out_prefix.name + f"_{hap_name}.bed")
        with out_path.open("w", encoding="utf-8") as out:
            if keep_headers:
                out.write("\n".join(keep_headers) + "\n")
            sub.to_csv(
                out,
                sep="\t",
                header=False,
                index=False,
                quoting=csv.QUOTE_NONE,
                escapechar="\\",
                lineterminator="\n",
            )
        LOGGER.info("Wrote %s", out_path)
        created_legacy.append(out_path)

    results["individual"] = (created_legacy[0], created_legacy[1])
    return results


# -----------------------------------------------------------------------------
# hap BEDs → final BED with colors & feature
# -----------------------------------------------------------------------------

def remove_chr_prefix(chrom: str) -> str:
    return strip_chr(str(chrom))


def _iter_bed_lines(bed_path: Path, geom_value: str) -> Iterable[Tuple[str, int, int, str]]:
    with bed_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line_no, line in enumerate(fh, start=1):
            if line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 4:
                continue
            chrom, s, e, ancestry = toks[0], toks[1], toks[2], toks[3]
            try:
                start = int(float(s))
                end = int(float(e))
            except ValueError:
                LOGGER.warning("Skipping malformed BED interval in %s:%d", bed_path, line_no)
                continue
            if end <= start:
                LOGGER.warning("Skipping non-positive BED interval in %s:%d (%s:%s-%s)", bed_path, line_no, chrom, s, e)
                continue
            yield (remove_chr_prefix(chrom), start, end, ancestry)


def load_haplotype_ancestry_intervals(bed_path: Path) -> Dict[str, List[Tuple[int, int, str]]]:
    intervals: Dict[str, List[Tuple[int, int, str]]] = {}
    for chrom, start, end, ancestry in _iter_bed_lines(bed_path, "geom_rect"):
        intervals.setdefault(chrom, []).append((start, end, ancestry))
    for chrom in intervals:
        intervals[chrom].sort(key=lambda item: (item[0], item[1]))
    return intervals


def lookup_ancestry_at_position(intervals: Dict[str, List[Tuple[int, int, str]]], chrom: str, pos: int) -> str:
    chrom_norm = remove_chr_prefix(chrom)
    for start, end, ancestry in intervals.get(chrom_norm, []):
        if start <= pos <= end:
            return ancestry
    return "."


def split_gt(gt: str) -> Tuple[str, str, bool, str]:
    gt = str(gt or ".")
    if "|" in gt:
        toks = gt.split("|")
        sep = "|"
        phased = True
    elif "/" in gt:
        toks = gt.split("/")
        sep = "/"
        phased = False
    else:
        toks = [gt]
        sep = ""
        phased = False
    hap1 = toks[0] if len(toks) > 0 else "."
    hap2 = toks[1] if len(toks) > 1 else "."
    return hap1, hap2, phased, sep


def allele_from_gt_token(token: str, ref: str, alt: str) -> Tuple[str, str, bool]:
    token = str(token or ".")
    if token in {"", "."}:
        return ".", "missing", False
    alt_alleles = str(alt or ".").split(",") if alt not in {"", "."} else []
    try:
        allele_idx = int(token)
    except ValueError:
        return token, "unknown", False
    if allele_idx == 0:
        return ref, "REF", False
    alt_type = "ALT" if len(alt_alleles) <= 1 and allele_idx == 1 else f"ALT{allele_idx}"
    if 1 <= allele_idx <= len(alt_alleles):
        return alt_alleles[allele_idx - 1], alt_type, True
    return ".", alt_type, True


def parse_float_or_none(value: Optional[str]) -> Optional[float]:
    if value in {None, "", ".", "NA", "nan"}:
        return None
    try:
        return float(str(value).split(",")[0])
    except ValueError:
        return None


def concordance_label(flare_ancestry: str, rfmix_ancestry: str) -> str:
    if flare_ancestry in {"", ".", "unknown"} or rfmix_ancestry in {"", "."}:
        return "NA"
    return "match" if flare_ancestry == rfmix_ancestry else "mismatch"


def snp_label_sort_key(label: str) -> Tuple[int, str]:
    if label.startswith("SNP") and label[3:].isdigit():
        return (int(label[3:]), "")
    return (10**9, label)


VARIANT_REPORT_COLUMNS = [
    "sample",
    "vcf_sample",
    "snp_label",
    "chrom",
    "pos",
    "variant_id",
    "REF",
    "ALT",
    "GT",
    "phased",
    "hap1_allele",
    "hap2_allele",
    "hap1_allele_type",
    "hap2_allele_type",
    "hap1_FLARE_ancestry",
    "hap2_FLARE_ancestry",
    "hap1_FLARE_code",
    "hap2_FLARE_code",
    "hap1_RFMix_segment_ancestry",
    "hap2_RFMix_segment_ancestry",
    "FLARE_SNP_ancestry",
    "RFMix_segment_ancestry",
    "ALT_allele_haplotype",
    "ALT_allele",
    "ALT_allele_FLARE_ancestry",
    "ALT_allele_RFMix_ancestry",
    "FLARE_RFMix_concordance",
    "ANP1",
    "ANP2",
    "confidence",
    "interpretation",
]


def write_variant_ancestry_report(
    sample: str,
    vcf_sample: str,
    annotations: Sequence[FlareSnpAnnotation],
    hap1_path: Path,
    hap2_path: Path,
    output_path: Path,
) -> None:
    hap1_intervals = load_haplotype_ancestry_intervals(hap1_path)
    hap2_intervals = load_haplotype_ancestry_intervals(hap2_path)

    sorted_annotations = sorted(
        annotations,
        key=lambda ann: (snp_label_sort_key(ann.snp_label), chrom_sort_key(ann.chrom), ann.pos, ann.variant_id),
    )

    with output_path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(VARIANT_REPORT_COLUMNS)

        for ann in sorted_annotations:
            hap1_token, hap2_token, phased, _sep = split_gt(ann.gt)
            hap1_allele, hap1_type, hap1_is_alt = allele_from_gt_token(hap1_token, ann.ref, ann.alt)
            hap2_allele, hap2_type, hap2_is_alt = allele_from_gt_token(hap2_token, ann.ref, ann.alt)
            hap1_rfmix = lookup_ancestry_at_position(hap1_intervals, ann.chrom, ann.pos)
            hap2_rfmix = lookup_ancestry_at_position(hap2_intervals, ann.chrom, ann.pos)

            hap_data = {
                "hap1": {
                    "allele": hap1_allele,
                    "allele_type": hap1_type,
                    "is_alt": hap1_is_alt,
                    "flare": ann.an1_label,
                    "rfmix": hap1_rfmix,
                    "anp": ann.anp1,
                },
                "hap2": {
                    "allele": hap2_allele,
                    "allele_type": hap2_type,
                    "is_alt": hap2_is_alt,
                    "flare": ann.an2_label,
                    "rfmix": hap2_rfmix,
                    "anp": ann.anp2,
                },
            }
            alt_haps = [hap for hap in ("hap1", "hap2") if hap_data[hap]["is_alt"]]

            if alt_haps and phased:
                alt_haplotype = ";".join(alt_haps)
                alt_allele = ";".join(f"{hap}={hap_data[hap]['allele']}" for hap in alt_haps)
                alt_flare = ";".join(f"{hap}={hap_data[hap]['flare']}" for hap in alt_haps)
                alt_rfmix = ";".join(f"{hap}={hap_data[hap]['rfmix']}" for hap in alt_haps)
                posterior_values = [parse_float_or_none(hap_data[hap]["anp"]) for hap in alt_haps]
                posterior_values = [p for p in posterior_values if p is not None]
                confidence = f"{min(posterior_values):.6g}" if posterior_values else "."
                interpretation = (
                    f"ALT allele is carried on {alt_haplotype}; in this individual, the carrier haplotype local ancestry "
                    f"assigned by FLARE is {alt_flare}. This describes the local ancestry background of the carrier "
                    "haplotype, not necessarily the evolutionary origin of the mutation."
                )
            elif alt_haps:
                alt_haplotype = "unresolved_unphased"
                alt_allele = ";".join(hap_data[hap]["allele"] for hap in alt_haps)
                alt_flare = "."
                alt_rfmix = "."
                confidence = "."
                interpretation = "ALT allele is observed, but GT is not phased; ALT haplotype ancestry cannot be assigned rigorously."
            else:
                alt_haplotype = "."
                alt_allele = "."
                alt_flare = "."
                alt_rfmix = "."
                confidence = "."
                interpretation = "No ALT allele is observed in this sample at this variant."

            concordance = (
                f"hap1={concordance_label(ann.an1_label, hap1_rfmix)};"
                f"hap2={concordance_label(ann.an2_label, hap2_rfmix)}"
            )

            writer.writerow([
                sample,
                vcf_sample,
                ann.snp_label,
                ann.chrom,
                ann.pos,
                ann.variant_id,
                ann.ref,
                ann.alt,
                ann.gt,
                "yes" if phased else "no",
                hap1_allele,
                hap2_allele,
                hap1_type,
                hap2_type,
                ann.an1_label,
                ann.an2_label,
                ann.an1_code if ann.an1_code is not None else ".",
                ann.an2_code if ann.an2_code is not None else ".",
                hap1_rfmix,
                hap2_rfmix,
                f"hap1={ann.an1_label};hap2={ann.an2_label}",
                f"hap1={hap1_rfmix};hap2={hap2_rfmix}",
                alt_haplotype,
                alt_allele,
                alt_flare,
                alt_rfmix,
                concordance,
                ann.anp1 or ".",
                ann.anp2 or ".",
                confidence,
                interpretation,
            ])
    LOGGER.info("Wrote LAP-VAR report %s", output_path)


@dataclass(frozen=True)
class FinalSegment:
    chrom: str
    start: int
    end: int
    feature: str
    color: str
    haplotype: str
    extra: Tuple[str, ...] = ()


def _final_line(seg: FinalSegment) -> str:
    base = [
        seg.chrom,
        str(seg.start),
        str(seg.end),
        seg.feature,
        seg.color,
        str(seg.haplotype),
    ]
    return "\t".join(base + [safe_bed_field(x) for x in seg.extra])


def merge_adjacent_rects(segments: Sequence[FinalSegment], max_gap_bp: int = 0) -> List[FinalSegment]:
    """Merge adjacent same-color intervals within the same chromosome and haplotype.

    RFMix2 can emit many consecutive windows with the same ancestry. Keeping every
    window as a separate SVG rectangle makes antialiasing seams more likely during
    PDF conversion, while merging same-call windows preserves the biological call.
    """
    if max_gap_bp < 0:
        raise ValueError("max_gap_bp must be >= 0")

    rects = [s for s in segments if s.feature == "geom_rect"]
    others = [s for s in segments if s.feature != "geom_rect"]

    rects.sort(key=lambda s: (chrom_sort_key(s.chrom), s.haplotype, s.start, s.end, s.color))
    merged: List[FinalSegment] = []

    for seg in rects:
        if not merged:
            merged.append(seg)
            continue

        last = merged[-1]
        same_track = (
            last.chrom == seg.chrom
            and last.feature == seg.feature
            and last.color == seg.color
            and last.haplotype == seg.haplotype
        )
        touches = seg.start <= last.end + max_gap_bp
        if same_track and touches:
            merged[-1] = FinalSegment(
                last.chrom,
                last.start,
                max(last.end, seg.end),
                last.feature,
                last.color,
                last.haplotype,
                last.extra,
            )
        else:
            merged.append(seg)

    return merged + others


def build_color_map(
    args: argparse.Namespace,
    ancestry_labels: Iterable[str],
    code_to_name: Optional[Dict[int, str]] = None,
) -> Dict[str, str]:
    """Build color lookup keyed by both ancestryN and declared ancestry names.

    This is important when RFMix/MSP and FLARE use different numeric ancestry codes:
    FLARE codes are first translated to names, then names are colored according to the
    MSP/RFMix palette.
    """
    cfg: Dict[str, str] = {}
    if args.color_config:
        cfg_path = Path(args.color_config)
        with cfg_path.open("r", encoding="utf-8") as fh:
            cfg = {str(k): str(v) for k, v in json.load(fh).items()}

    colors: Dict[str, str] = {}
    code_to_name = code_to_name or {}
    max_code = max([14] + list(code_to_name.keys()))
    for i in range(0, max_code + 1):
        default_color = getattr(args, f"ancestry{i}", args.unknown)
        ancestry_key = f"ancestry{i}"
        colors[ancestry_key] = cfg.get(ancestry_key, default_color)
        if i in code_to_name:
            name = code_to_name[i]
            colors[name] = cfg.get(name, cfg.get(ancestry_key, default_color))

    for lab in ancestry_labels:
        colors.setdefault(str(lab), cfg.get(str(lab), args.unknown))

    # Explicit config entries not otherwise covered are still accepted.
    for name, color in cfg.items():
        colors[name] = color
    return colors


def color_for_ancestry_label(label: str, code: Optional[int], colors: Dict[str, str], unknown_color: str) -> str:
    """Return the color for a FLARE ancestry label using names before numeric codes.

    FLARE and RFMix/MSP can use different numeric ancestry orders. Once FLARE has
    provided a real ancestry name, do not silently fall back to ancestry{code};
    that would reintroduce code-order mismatches.
    """
    if label in colors:
        return colors[label]

    label_s = str(label)
    synthetic_label = code is not None and label_s == f"ancestry{code}"
    unknown_label = label_s.strip().lower() in {"", ".", "na", "nan", "unknown"}
    if code is not None and (synthetic_label or unknown_label) and f"ancestry{code}" in colors:
        return colors[f"ancestry{code}"]
    return unknown_color


def process_haps_to_final(
    hap1_path: Path,
    hap2_path: Path,
    colors: Dict[str, str],
    unknown_color: str,
    features: Optional[Sequence[Feature]],
    final_output_path: Path,
    header_from: Optional[List[str]] = None,
    merge_adjacent: bool = True,
    max_merge_gap_bp: int = 0,
    snp_annotations: Optional[Sequence[FlareSnpAnnotation]] = None,
) -> None:
    segments: List[FinalSegment] = []

    for bed_path, gv in ((hap1_path, "1"), (hap2_path, "2")):
        for chrom, start, end, ancestry in _iter_bed_lines(bed_path, gv):
            color = colors.get(str(ancestry), unknown_color)
            segments.append(FinalSegment(chrom, start, end, "geom_rect", color, gv))

    for feature in features or []:
        fchrom = remove_chr_prefix(feature.chrom)
        segments.extend([
            FinalSegment(fchrom, feature.start_bp, feature.end_bp, "geom_line", "#000000", "1", (feature.label,)),
            FinalSegment(fchrom, feature.start_bp, feature.end_bp, "geom_line", "#000000", "2", (feature.label,)),
        ])

    if snp_annotations:
        for ann in snp_annotations:
            color1 = color_for_ancestry_label(ann.an1_label, ann.an1_code, colors, unknown_color)
            color2 = color_for_ancestry_label(ann.an2_label, ann.an2_code, colors, unknown_color)
            common_extra = (ann.variant_id, ann.vcf_sample)
            segments.extend([
                FinalSegment(
                    ann.chrom,
                    ann.pos,
                    ann.pos,
                    "geom_snp",
                    color1,
                    "1",
                    common_extra + (ann.an1_label, str(ann.an1_code) if ann.an1_code is not None else ".", ann.snp_label),
                ),
                FinalSegment(
                    ann.chrom,
                    ann.pos,
                    ann.pos,
                    "geom_snp",
                    color2,
                    "2",
                    common_extra + (ann.an2_label, str(ann.an2_code) if ann.an2_code is not None else ".", ann.snp_label),
                ),
            ])

    raw_count = len([s for s in segments if s.feature == "geom_rect"])
    if merge_adjacent:
        segments = merge_adjacent_rects(segments, max_gap_bp=max_merge_gap_bp)
        merged_count = len([s for s in segments if s.feature == "geom_rect"])
        if merged_count != raw_count:
            LOGGER.info("Merged %d ancestry rectangles into %d for %s", raw_count, merged_count, final_output_path.name)

    def sort_key(seg: FinalSegment) -> Tuple[Tuple[int, str], int, int, str]:
        ckey = chrom_sort_key(seg.chrom)
        hap_order = int(seg.haplotype) if str(seg.haplotype).isdigit() else 99
        return (ckey, seg.start, hap_order, seg.feature)

    segments.sort(key=sort_key)
    lines = [_final_line(seg) for seg in segments]

    with final_output_path.open("w", encoding="utf-8") as out:
        if header_from:
            out.write("\n".join(header_from) + "\n")
        out.write("\n".join(lines) + "\n")
    LOGGER.info("Wrote %s", final_output_path)


# -----------------------------------------------------------------------------
# Combine per-chrom MSPs → one MSP per individual
# -----------------------------------------------------------------------------

def combine_msp_for_individual(files: Sequence[Path], dest_path: Path, require_headers_from_first: bool = True) -> List[str]:
    """Concatenate MSPs skipping comment headers after the first file.

    Returns the header lines captured from the first file.
    """
    headers: List[str] = []
    wrote_header = False

    with dest_path.open("w", encoding="utf-8") as out:
        for i, fp in enumerate(files):
            LOGGER.debug("Concatenating %s (%d/%d)", fp, i + 1, len(files))
            with fp.open("r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if line.startswith("#"):
                        if not wrote_header:
                            headers.append(line.rstrip("\n"))
                            out.write(line)
                        # Comment/header lines are copied only once, from the first chunk.
                        continue
                    out.write(line)
            if not wrote_header:
                wrote_header = True
    return headers


# -----------------------------------------------------------------------------
# Orchestration per individual
# -----------------------------------------------------------------------------

def process_individual(
    individual: str,
    files: Sequence[Path],
    out_dir: Path,
    args: argparse.Namespace,
    features: Sequence[Feature],
    keep_temp: bool,
    merge_adjacent: bool,
    max_merge_gap_bp: int,
    flare_annotations: Optional[FlareAnnotationSet] = None,
    sample_map: Optional[Dict[str, str]] = None,
) -> List[Tuple[str, Path]]:
    LOGGER.info("Processing individual: %s", individual)
    combined_path = out_dir / f"{individual}{MSP_SUFFIX}"

    # Combine MSPs
    headers = combine_msp_for_individual(files, combined_path)
    msp_ancestry_map = build_ancestry_map_from_headers(headers)
    ancestry_labels = set(msp_ancestry_map.values()) | {f"ancestry{i}" for i in range(0, 15)}
    colors = build_color_map(args, ancestry_labels, code_to_name=msp_ancestry_map)

    headers_bed = filter_headers_for_bed(headers)
    color_header = build_lap_color_header(msp_ancestry_map, colors)
    if color_header and color_header not in headers_bed:
        headers_bed.append(color_header)
    if flare_annotations is not None:
        headers_bed.append("#LAP FLARE SNP columns: chrom start end geom_snp color haplotype variant_id vcf_sample ancestry_label ancestry_code snp_label")

    # Convert to hap beds (may yield many samples)
    out_prefix = combined_path.with_suffix("")
    sample_haps = msp_to_hap_beds(combined_path, out_prefix, msp_ancestry_map, keep_headers=headers_bed)

    results: List[Tuple[str, Path]] = []
    temp_paths: List[Path] = [combined_path]
    sample_map = sample_map or {}

    for sample, (hap1, hap2) in sample_haps.items():
        sample_safe = safe_sample_name(sample)
        final_bed = out_dir / f"{sample_safe}.bed"

        sample_snp_annotations: List[FlareSnpAnnotation] = []
        vcf_sample = resolve_vcf_sample(sample, flare_annotations, sample_map, args.flare_sample)
        if flare_annotations is not None:
            if vcf_sample is None:
                LOGGER.warning("No FLARE sample match for LAP sample '%s'; SNP annotation skipped for this sample.", sample)
            else:
                sample_snp_annotations = flare_annotations.by_sample.get(vcf_sample, [])
                LOGGER.info(
                    "Adding %d FLARE SNP annotation(s) to %s using VCF sample %s.",
                    len(sample_snp_annotations),
                    sample,
                    vcf_sample,
                )

        process_haps_to_final(
            hap1,
            hap2,
            colors,
            args.unknown,
            features,
            final_bed,
            header_from=headers_bed,
            merge_adjacent=merge_adjacent,
            max_merge_gap_bp=max_merge_gap_bp,
            snp_annotations=sample_snp_annotations,
        )
        if flare_annotations is not None and vcf_sample is not None:
            report_path = out_dir / f"{sample_safe}.variant_ancestry.txt"
            write_variant_ancestry_report(
                sample,
                vcf_sample,
                sample_snp_annotations,
                hap1,
                hap2,
                report_path,
            )
        results.append((sample, final_bed))
        temp_paths.extend([hap1, hap2])

    # Cleanup
    if not keep_temp:
        for pth in temp_paths:
            try:
                pth.unlink()
                LOGGER.debug("Removed %s", pth)
            except Exception as e:
                LOGGER.warning("Could not remove %s: %s", pth, e)

    return results


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine RFMix2 MSP chunks into final ancestry BEDs with colors and optional feature lines."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Discovery / I/O
    parser.add_argument("--prefix", required=True, help="Prefix for RFMix2 MSP files (may include path).")
    parser.add_argument(
        "--chr",
        nargs="+",
        default=[f"chr{i}" for i in range(1, 23)],
        help="List of chromosome tokens present in filenames (e.g., chr1 chr2 ... OR 1 2 ...)",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory.")

    # Behavior toggles
    parser.add_argument("--require-all-chroms", action="store_true", help="Skip individuals missing any requested chromosome.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep combined MSP and hap BEDs.")
    parser.add_argument("--dry-run", action="store_true", help="Only list planned operations, do not process.")
    parser.add_argument("--threads", type=int, default=1, help="Parallel individuals (I/O-bound, ThreadPool).")
    parser.add_argument(
        "--no-merge-adjacent",
        action="store_true",
        help="Keep consecutive same-ancestry RFMix windows as separate BED records.",
    )
    parser.add_argument(
        "--max-merge-gap-bp",
        type=int,
        default=0,
        help="Allow merging same-ancestry windows separated by up to this many bp.",
    )

    # Colors & mapping
    parser.add_argument("--color-config", type=str, help="JSON file mapping ancestry labels to colors.")
    parser.add_argument("--unknown", default="#808080", help="Color for unknown/unspecified ancestries.")
    for i, default in enumerate(DEFAULT_ANCESTRY_COLORS):
        parser.add_argument(f"--ancestry{i}", default=default, help=f"Color for ancestry{i}")

    # Optional manual feature
    parser.add_argument("--from-bp", type=int, dest="from_bp", help="Feature start (bp)")
    parser.add_argument("--to-bp", type=int, dest="to_bp", help="Feature end (bp)")
    parser.add_argument("-c", "--chromosome", type=str, help="Chromosome for feature (e.g., chr2 or 2)")
    parser.add_argument(
        "--highlight-regions",
        "--regions-to-highlight",
        dest="highlight_regions",
        type=str,
        help="TXT/TSV file of regions to highlight with dashed black lines: chrom start end [label].",
    )

    # FLARE SNP annotations
    parser.add_argument("--flare-vcf", type=str, help="FLARE .anc.vcf or .anc.vcf.gz containing FORMAT/AN1 and FORMAT/AN2.")
    parser.add_argument("--snp-list", type=str, help="TXT/TSV list of SNPs to annotate: CHROM POS, CHROM:POS, or VCF IDs.")
    parser.add_argument("--flare-sample", type=str, help="Use this VCF sample for all LAP output samples; useful for legacy single-sample MSPs.")
    parser.add_argument("--flare-sample-map", type=str, help="Two-column map: LAP_sample<TAB>VCF_sample.")

    # Logging
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--debug", action="store_true", help="Verbose debug logs.")
    g.add_argument("--quiet", action="store_true", help="Only warnings and errors.")

    return parser.parse_args(argv)


def setup_logging(args: argparse.Namespace) -> None:
    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")


def build_ancestry_map() -> Dict[int, str]:
    # Backward-compatible alias retained for external imports.
    return build_default_ancestry_map()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    setup_logging(args)

    out_dir = Path(args.output_dir)
    ensure_dir(out_dir)

    # Normalize chrom tokens to include both '1' and 'chr1' variants when searching
    chroms: List[str] = []
    for tok in args.chr:
        tok = str(tok)
        if tok.lower().startswith("chr"):
            chroms.append(tok)
            chroms.append(strip_chr(tok))
        else:
            chroms.append(tok)
            chroms.append(f"chr{tok}")
    # Deduplicate while preserving order
    seen = set()
    chroms = [x for x in chroms if not (x in seen or seen.add(x))]

    individuals_map, individuals = discover_inputs(args.prefix, chroms)
    if not individuals:
        LOGGER.error("No individuals found for prefix '%s' and chromosomes %s", args.prefix, args.chr)
        return 2

    # Optionally filter individuals missing requested chroms
    if args.require_all_chroms:
        req_set = {c for c in chroms if str(c).lower().startswith("chr")}
        def has_all(files: Sequence[Path]) -> bool:
            present = { _guess_chrom_from_filename(f.name).lower() for f in files }
            return req_set.issubset(present)
        filtered = {ind: fs for ind, fs in individuals_map.items() if has_all(fs)}
        missing = sorted(set(individuals_map) - set(filtered))
        if missing:
            LOGGER.warning("Skipping %d individuals missing chromosomes: %s", len(missing), ", ".join(missing))
        individuals_map = filtered
        individuals = sorted(individuals_map)
        if not individuals:
            LOGGER.error("After filtering, no individuals have all requested chromosomes.")
            return 3

    LOGGER.info("Individuals to process: %s", ", ".join(individuals))

    # Dry run: list planned operations and exit
    if args.dry_run:
        for ind in individuals:
            files = individuals_map[ind]
            LOGGER.info("[DRY‑RUN] %s: %d MSP chunks → (per-sample BEDs in %s)", ind, len(files), out_dir)
        return 0

    features: List[Feature] = []
    if args.from_bp is not None and args.to_bp is not None and args.chromosome is not None:
        features.append(Feature(chrom=args.chromosome, start_bp=args.from_bp, end_bp=args.to_bp, label="manual_feature"))
    elif any(value is not None for value in (args.from_bp, args.to_bp, args.chromosome)):
        LOGGER.warning("Ignoring incomplete manual feature; provide -c/--chromosome, --from-bp, and --to-bp together.")

    if args.highlight_regions:
        features.extend(load_highlight_regions(Path(args.highlight_regions)))

    sample_map = load_sample_map(Path(args.flare_sample_map)) if args.flare_sample_map else {}
    target_samples: Optional[Set[str]] = None
    if args.flare_sample:
        target_samples = {args.flare_sample}
    elif sample_map:
        target_samples = set(sample_map.values())
    elif args.flare_vcf:
        inferred_samples = infer_lap_sample_names(individuals_map)
        target_samples = inferred_samples or None
        if inferred_samples:
            LOGGER.info("Restricting FLARE parsing to %d LAP sample name(s).", len(inferred_samples))

    flare_annotations: Optional[FlareAnnotationSet] = None
    if args.flare_vcf:
        flare_annotations = load_flare_annotations(
            Path(args.flare_vcf),
            Path(args.snp_list) if args.snp_list else None,
            target_samples=target_samples,
        )

    results: List[Tuple[str, Path]] = []

    if args.threads and args.threads > 1:
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = {
                ex.submit(
                    process_individual,
                    ind,
                    individuals_map[ind],
                    out_dir,
                    args,
                    features,
                    args.keep_temp,
                    not args.no_merge_adjacent,
                    args.max_merge_gap_bp,
                    flare_annotations,
                    sample_map,
                ): ind
                for ind in individuals
            }
            for fut in as_completed(futs):
                try:
                    res_list = fut.result()
                    results.extend(res_list)
                except Exception as e:
                    LOGGER.error("Failed individual %s: %s", futs[fut], e)
                    return 4
    else:
        for ind in individuals:
            res = process_individual(
                ind,
                individuals_map[ind],
                out_dir,
                args,
                features,
                args.keep_temp,
                not args.no_merge_adjacent,
                args.max_merge_gap_bp,
                flare_annotations,
                sample_map,
            )
            results.extend(res)

    LOGGER.info("Done. Generated %d BED files in %s", len(results), out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
