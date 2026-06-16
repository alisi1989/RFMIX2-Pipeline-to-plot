#!/usr/bin/env python3
# Authors: Alessandro Lisi & Michael C. Campbell
# Be sure to install librsvg beforehand:
#   macOS: brew install librsvg
#   Linux: install the package that provides rsvg-convert

"""Render LAP BED ancestry segments on the hg38 chromosome SVG template."""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


LOGGER = logging.getLogger("plot_lap")


ET.register_namespace("", "http://www.w3.org/2000/svg")
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
ET.register_namespace("inkscape", "http://www.inkscape.org/namespaces/inkscape")
ET.register_namespace("sodipodi", "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd")


CHROMOSOME_COORDINATES = {
    "chromosome1_hap1": {"x": 217.900, "y": 103.20, "width": 27, "height": 666.300},
    "chromosome1_hap2": {"x": 248.500, "y": 103.20, "width": 27, "height": 666.300},
    "chromosome2_hap1": {"x": 290.900, "y": 121.300, "width": 27, "height": 648.200},
    "chromosome2_hap2": {"x": 321.400, "y": 121.300, "width": 27, "height": 648.200},
    "chromosome3_hap1": {"x": 363.800, "y": 238.600, "width": 27, "height": 530.900},
    "chromosome3_hap2": {"x": 394.000, "y": 238.600, "width": 27, "height": 530.900},
    "chromosome4_hap1": {"x": 436.500, "y": 260.200, "width": 27, "height": 509.300},
    "chromosome4_hap2": {"x": 467.000, "y": 260.200, "width": 27, "height": 509.300},
    "chromosome5_hap1": {"x": 509.300, "y": 283.300, "width": 27, "height": 486.200},
    "chromosome5_hap2": {"x": 540.100, "y": 283.300, "width": 27, "height": 486.200},
    "chromosome6_hap1": {"x": 582.100, "y": 312.000, "width": 27, "height": 457.500},
    "chromosome6_hap2": {"x": 613.000, "y": 312.000, "width": 27, "height": 457.500},
    "chromosome7_hap1": {"x": 654.900, "y": 342.600, "width": 27, "height": 426.900},
    "chromosome7_hap2": {"x": 685.800, "y": 342.600, "width": 27, "height": 426.900},
    "chromosome8_hap1": {"x": 727.600, "y": 380.600, "width": 27, "height": 388.900},
    "chromosome8_hap2": {"x": 758.500, "y": 380.600, "width": 27, "height": 388.900},
    "chromosome9_hap1": {"x": 800.200, "y": 398.600, "width": 27, "height": 370.900},
    "chromosome9_hap2": {"x": 831.100, "y": 398.600, "width": 27, "height": 370.900},
    "chromosome10_hap1": {"x": 873.200, "y": 410.800, "width": 27, "height": 358.600},
    "chromosome10_hap2": {"x": 904.100, "y": 410.800, "width": 27, "height": 358.600},
    "chromosome11_hap1": {"x": 946.000, "y": 407.300, "width": 27, "height": 362.000},
    "chromosome11_hap2": {"x": 976.900, "y": 407.300, "width": 27, "height": 362.000},
    "chromosome12_hap1": {"x": 1018.700, "y": 412.300, "width": 27, "height": 357.200},
    "chromosome12_hap2": {"x": 1049.700, "y": 412.300, "width": 27, "height": 357.200},
    "chromosome13_hap1": {"x": 1092.000, "y": 462.600, "width": 27, "height": 306.700},
    "chromosome13_hap2": {"x": 1123.000, "y": 462.600, "width": 27, "height": 306.700},
    "chromosome14_hap1": {"x": 1164.000, "y": 482.400, "width": 27, "height": 287.100},
    "chromosome14_hap2": {"x": 1195.000, "y": 482.400, "width": 27, "height": 287.100},
    "chromosome15_hap1": {"x": 1237.300, "y": 495.700, "width": 27, "height": 273.600},
    "chromosome15_hap2": {"x": 1268.300, "y": 495.700, "width": 27, "height": 273.600},
    "chromosome16_hap1": {"x": 1309.800, "y": 527.100, "width": 27, "height": 242.400},
    "chromosome16_hap2": {"x": 1340.800, "y": 527.100, "width": 27, "height": 242.400},
    "chromosome17_hap1": {"x": 1382.600, "y": 545.700, "width": 27, "height": 223.500},
    "chromosome17_hap2": {"x": 1413.500, "y": 545.700, "width": 27, "height": 223.500},
    "chromosome18_hap1": {"x": 1456.000, "y": 553.500, "width": 27, "height": 215.800},
    "chromosome18_hap2": {"x": 1486.100, "y": 553.500, "width": 27, "height": 215.800},
    "chromosome19_hap1": {"x": 1528.500, "y": 611.800, "width": 27, "height": 157.700},
    "chromosome19_hap2": {"x": 1559.100, "y": 611.800, "width": 27, "height": 157.700},
    "chromosome20_hap1": {"x": 1601.600, "y": 596.200, "width": 27, "height": 173.300},
    "chromosome20_hap2": {"x": 1631.700, "y": 596.200, "width": 27, "height": 173.300},
    "chromosome21_hap1": {"x": 1674.400, "y": 643.300, "width": 27, "height": 125.900},
    "chromosome21_hap2": {"x": 1704.600, "y": 643.300, "width": 27, "height": 125.900},
    "chromosome22_hap1": {"x": 1747.200, "y": 632.700, "width": 27, "height": 136.800},
    "chromosome22_hap2": {"x": 1777.400, "y": 632.700, "width": 27, "height": 136.800},
}

CHROMOSOME_LENGTHS = OrderedDict(
    [
        ("1", 248956422),
        ("2", 242193529),
        ("3", 198295559),
        ("4", 190214555),
        ("5", 181538259),
        ("6", 170805979),
        ("7", 159345973),
        ("8", 145138636),
        ("9", 138394717),
        ("10", 133797422),
        ("11", 135086622),
        ("12", 133275309),
        ("13", 114364328),
        ("14", 107043718),
        ("15", 101991189),
        ("16", 90338345),
        ("17", 83257441),
        ("18", 80373285),
        ("19", 58617616),
        ("20", 64444167),
        ("21", 46709983),
        ("22", 50818468),
    ]
)


DEFAULT_ANCESTRY_COLORS = [
    ("ancestry0", "#a32e2e"),
    ("ancestry1", "#0a0ae0"),
    ("ancestry2", "#bfa004"),
    ("ancestry3", "#d18311"),
    ("ancestry4", "#22ba9d"),
    ("ancestry5", "#839dfc"),
    ("ancestry6", "#9a5dc1"),
    ("ancestry7", "#26962b"),
    ("ancestry8", "#707070"),
    ("ancestry9", "#00cfff"),
    ("ancestry10", "#790ee0"),
    ("ancestry11", "#ff4d6d"),
    ("ancestry12", "#2d6a4f"),
    ("ancestry13", "#f77f00"),
    ("ancestry14", "#4ea8de"),
]

SVG_VIEWBOX_WIDTH = 1800
SVG_VIEWBOX_HEIGHT = 864
PDF_POINTS_PER_SVG_UNIT = 72.0 / 100.0
NATIVE_CANVAS_LEFT_PADDING_PX = 100.0
NATIVE_CANVAS_RIGHT_PADDING_PX = 48.0
NATIVE_LEGEND_TEXT_WIDTH_FACTOR = 0.62
CENTROMERE_OVERLAY_SEGMENT_WIDTH_PX = 2.4
CENTROMERE_OVERLAY_SEGMENT_GAP_PX = 0.9
LEGEND_SWATCH_X = 1450
LEGEND_TEXT_X = 1486
LEGEND_SWATCH_WIDTH = 30
LEGEND_SWATCH_HEIGHT = 18
LEGEND_FONT_SIZE = 17
LEGEND_ROW_STEP = 32
CENTROMERE_LEGEND_LABEL = "Centromere"
ACROCENTRIC_P_ARM_COLOR = "#d3d3d3"
ACROCENTRIC_P_ARM_LEGEND_LABEL = "Acrocentric p-arms/satellites"
SNP_LABEL_HALO_WIDTH_PX = 1.5
ZOOM_PANEL_WIDTH = 900.0
ZOOM_PANEL_HEIGHT = 760.0
ZOOM_CONTEXT_X = 120.0
ZOOM_CONTEXT_Y = 144.0
ZOOM_CONTEXT_WIDTH = 18.0
ZOOM_CONTEXT_HEIGHT = 470.0
ZOOM_TRACK_TOP = 104.0
ZOOM_TRACK_HEIGHT = 540.0
ZOOM_HAP1_X = 360.0
ZOOM_HAP2_X = 448.0
ZOOM_TRACK_WIDTH = 48.0
ZOOM_AXIS_X = 292.0
ZOOM_LABEL_GAP = 15.0
ZOOM_LABEL_MIN_SEPARATION_PX = 16.0
ZOOM_DEFAULT_OVERLAP_PX = 18.0
ZOOM_DEFAULT_PADDING_BP = 250000
ZOOM_DEFAULT_MIN_WINDOW_BP = 500000
ZOOM_DEFAULT_MAX_PANELS = 12

CENTROMERE_BP = {
    "1": (121700000, 125100000),
    "2": (91800000, 96000000),
    "3": (87800000, 94000000),
    "4": (48200000, 51800000),
    "5": (46100000, 51400000),
    "6": (58500000, 62600000),
    "7": (58100000, 62100000),
    "8": (43200000, 47200000),
    "9": (42200000, 45500000),
    "10": (38000000, 41600000),
    "11": (51000000, 55800000),
    "12": (33200000, 37800000),
    "13": (16500000, 18900000),
    "14": (16100000, 18200000),
    "15": (17500000, 20500000),
    "16": (35300000, 38400000),
    "17": (22700000, 27400000),
    "18": (15400000, 21500000),
    "19": (24200000, 28100000),
    "20": (25700000, 30400000),
    "21": (10900000, 13000000),
    "22": (13700000, 17400000),
}

ACROCENTRIC_CHROMS = {"13", "14", "15", "21", "22"}
LINE_GEOMS = {"geom_line", "geom_snp"}
LEGEND_GEOMS = {"geom_rect", "geom_snp"}


@dataclass(frozen=True)
class BedRecord:
    chrom: str
    start: int
    end: int
    line_type: str
    color: str
    haplotype: int
    source_line: int
    extras: Tuple[str, ...] = ()


@dataclass(frozen=True)
class SnpGroup:
    chrom: str
    pos: int
    variant_id: str
    label: str
    source_line: int
    records: Tuple[BedRecord, ...]


@dataclass(frozen=True)
class SnpOverlapCluster:
    chrom: str
    start: int
    end: int
    groups: Tuple[SnpGroup, ...]


def normalize_chromosome(chromosome: str) -> str:
    chrom = str(chromosome).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def chrom_sort_key(chromosome: str) -> Tuple[int, str]:
    chrom = normalize_chromosome(chromosome).upper()
    try:
        return (int(chrom), "")
    except ValueError:
        return (1000, chrom)


def parse_svg_filename(filename: str) -> str:
    return filename if filename.endswith(".svg") else filename + ".svg"


def parse_output_paths(out_arg: str) -> Tuple[str, str]:
    """Return (output_svg_path, output_pdf_path) from a user-provided -O argument."""
    if out_arg.lower().endswith(".pdf"):
        base = os.path.splitext(out_arg)[0]
        return base + ".svg", out_arg
    if out_arg.lower().endswith(".svg"):
        base = os.path.splitext(out_arg)[0]
        return out_arg, base + ".pdf"
    return out_arg + ".svg", out_arg + ".pdf"


def read_header_lines(bed_file: Path) -> List[str]:
    headers: List[str] = []
    with bed_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            headers.append(line.strip())
    return headers


def parse_subpopulation_codes(headers: Sequence[str]) -> Dict[int, str]:
    for header in headers:
        if not header.startswith("#Subpopulation order/codes:"):
            continue
        payload = header.split(":", 1)[1].strip()
        out: Dict[int, str] = {}
        for item in re.split(r"\s+", payload):
            if "=" not in item:
                continue
            name, code = item.rsplit("=", 1)
            try:
                out[int(code)] = name
            except ValueError:
                LOGGER.warning("Could not parse ancestry code in header token: %s", item)
        return out
    return {}


def load_color_config(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    return {str(k): str(v) for k, v in cfg.items()}


def parse_lap_color_header(headers: Sequence[str]) -> "OrderedDict[str, str]":
    """Parse the optional LAP-generated palette header.

    Format: #LAP ancestry/colors:	AncestryName=#rrggbb	...
    """
    palette: "OrderedDict[str, str]" = OrderedDict()
    for header in headers:
        if not header.startswith("#LAP ancestry/colors:"):
            continue
        payload = header.split(":", 1)[1].strip()
        for item in re.split(r"\s+", payload):
            if "=" not in item:
                continue
            name, color = item.rsplit("=", 1)
            if name and color:
                palette[name] = color
        break
    return palette


def build_ancestry_palette(headers: Sequence[str], color_config: Optional[str]) -> "OrderedDict[str, str]":
    code_to_name = parse_subpopulation_codes(headers)
    cfg = load_color_config(color_config)
    lap_palette = parse_lap_color_header(headers)

    palette: "OrderedDict[str, str]" = OrderedDict()
    for code, (default_name, default_color) in enumerate(DEFAULT_ANCESTRY_COLORS):
        name = code_to_name.get(code, default_name)
        color = lap_palette.get(name, cfg.get(name, cfg.get(default_name, default_color)))
        palette[name] = color

    # Preserve any additional palette entries produced by LAP, then explicit config entries.
    for name, color in lap_palette.items():
        palette.setdefault(name, color)
    for name, color in cfg.items():
        palette.setdefault(name, color)
    return palette


def parse_bed_records(bed_file: Path, strict: bool = False) -> List[BedRecord]:
    records: List[BedRecord] = []
    bad_lines = 0

    with bed_file.open("r", encoding="utf-8", errors="ignore") as bed:
        for line_no, raw_line in enumerate(bed, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                bad_lines += 1
                LOGGER.warning("Skipping short BED line %s:%d", bed_file, line_no)
                continue

            chromosome, start, end, line_type, color, haplotype = parts[:6]
            chrom = normalize_chromosome(chromosome)
            if chrom not in CHROMOSOME_LENGTHS:
                msg = f"Unsupported chromosome in {bed_file}:{line_no}: {chromosome}"
                if strict:
                    raise ValueError(msg)
                bad_lines += 1
                LOGGER.warning("%s", msg)
                continue

            try:
                start_i = int(float(start))
                end_i = int(float(end))
                hap_i = int(haplotype)
            except ValueError:
                bad_lines += 1
                LOGGER.warning("Skipping malformed BED line %s:%d", bed_file, line_no)
                continue

            if hap_i not in (1, 2):
                bad_lines += 1
                LOGGER.warning("Skipping unsupported haplotype in %s:%d: %s", bed_file, line_no, haplotype)
                continue
            if line_type in LINE_GEOMS:
                if end_i < start_i:
                    bad_lines += 1
                    LOGGER.warning("Skipping negative line interval in %s:%d", bed_file, line_no)
                    continue
            elif end_i <= start_i:
                bad_lines += 1
                LOGGER.warning("Skipping non-positive interval in %s:%d", bed_file, line_no)
                continue

            records.append(BedRecord(chrom, start_i, end_i, line_type, color, hap_i, line_no, tuple(parts[6:])))

    if bad_lines:
        LOGGER.warning("Skipped %d malformed/unsupported BED lines", bad_lines)
    return sorted(records, key=lambda r: (chrom_sort_key(r.chrom), r.start, r.haplotype, r.end))


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[Tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def log_coverage_summary(records: Sequence[BedRecord]) -> None:
    chroms_present = sorted({r.chrom for r in records if r.line_type == "geom_rect"}, key=chrom_sort_key)
    if not chroms_present:
        LOGGER.warning("No drawable geom_rect records found in BED.")
        return

    missing = [chrom for chrom in CHROMOSOME_LENGTHS if chrom not in chroms_present]
    LOGGER.info(
        "BED contains drawable ancestry for chromosomes: %s (%d/%d).",
        ", ".join(chroms_present),
        len(chroms_present),
        len(CHROMOSOME_LENGTHS),
    )
    if missing:
        LOGGER.info("Chromosomes with no BED ancestry records will remain empty: %s", ", ".join(missing))

    by_track: Dict[Tuple[str, int], List[Tuple[int, int]]] = {}
    for record in records:
        if record.line_type != "geom_rect":
            continue
        chrom_len = CHROMOSOME_LENGTHS[record.chrom]
        start = max(0, min(record.start, chrom_len))
        end = max(0, min(record.end, chrom_len))
        if end > start:
            by_track.setdefault((record.chrom, record.haplotype), []).append((start, end))

    for (chrom, haplotype), intervals in sorted(by_track.items(), key=lambda kv: (chrom_sort_key(kv[0][0]), kv[0][1])):
        covered = sum(end - start for start, end in merge_intervals(intervals))
        pct = covered / CHROMOSOME_LENGTHS[chrom] * 100
        if pct < 98:
            LOGGER.warning("chr%s hap%d covers %.2f%% of hg38 length; visible blank regions are expected.", chrom, haplotype, pct)


def svg_child_insert_index(root: ET.Element) -> int:
    idx = 0
    for child in list(root):
        tag = child.tag.split("}", 1)[-1]
        if tag in {"style", "namedview", "defs"}:
            idx += 1
            continue
        break
    return idx


def add_rect(
    parent: ET.Element,
    x: float,
    y: float,
    width: float,
    height: float,
    color: str,
    element_id: Optional[str] = None,
    opacity: Optional[float] = None,
) -> None:
    attrs = {
        "x": f"{x:.6f}",
        "y": f"{y:.6f}",
        "width": f"{width:.6f}",
        "height": f"{height:.6f}",
        "fill": color,
        "stroke": "none",
        "shape-rendering": "crispEdges",
    }
    if element_id:
        attrs["id"] = element_id
    if opacity is not None:
        attrs["opacity"] = f"{max(0.0, min(1.0, opacity)):.3f}"
    parent.append(ET.Element("rect", attrs))


def segmented_spans(start: float, length: float) -> Iterable[Tuple[float, float]]:
    end = start + length
    if length <= CENTROMERE_OVERLAY_SEGMENT_WIDTH_PX * 1.4:
        yield start, length
        return

    cursor = start
    while cursor < end:
        seg_end = min(cursor + CENTROMERE_OVERLAY_SEGMENT_WIDTH_PX, end)
        if seg_end > cursor:
            yield cursor, seg_end - cursor
        cursor += CENTROMERE_OVERLAY_SEGMENT_WIDTH_PX + CENTROMERE_OVERLAY_SEGMENT_GAP_PX


def add_segmented_rect(
    parent: ET.Element,
    x: float,
    y: float,
    width: float,
    height: float,
    color: str,
    element_id: Optional[str] = None,
    opacity: Optional[float] = None,
) -> None:
    for idx, (seg_x, seg_width) in enumerate(segmented_spans(x, width), start=1):
        seg_id = f"{element_id}_seg{idx}" if element_id else None
        add_rect(parent, seg_x, y, seg_width, height, color, element_id=seg_id, opacity=opacity)


def add_centromere_legend_swatch(parent: ET.Element, y: float, color: str, opacity: float) -> None:
    parent.append(
        ET.Element(
            "rect",
            {
                "x": f"{LEGEND_SWATCH_X:.6f}",
                "y": f"{y:.6f}",
                "width": f"{LEGEND_SWATCH_WIDTH:.6f}",
                "height": f"{LEGEND_SWATCH_HEIGHT:.6f}",
                "fill": "#ffffff",
                "stroke": "#666666",
                "stroke-width": "0.7",
            },
        )
    )
    add_segmented_rect(
        parent,
        LEGEND_SWATCH_X,
        y,
        LEGEND_SWATCH_WIDTH,
        LEGEND_SWATCH_HEIGHT,
        color,
        element_id="lap_centromere_legend",
        opacity=opacity,
    )


def svg_uses_acrocentric_p_arm_color(root: ET.Element) -> bool:
    svg_text = ET.tostring(root, encoding="unicode").lower()
    return ACROCENTRIC_P_ARM_COLOR.lower() in svg_text


def add_feature_lines(parent: ET.Element, record: BedRecord, start_y: float, end_y: float, x: float, width: float) -> None:
    offset = 15
    for y in (start_y, end_y):
        parent.append(
            ET.Element(
                "line",
                {
                    "x1": f"{x - offset:.6f}",
                    "y1": f"{y:.6f}",
                    "x2": f"{x + width + offset:.6f}",
                    "y2": f"{y:.6f}",
                    "stroke": record.color,
                    "stroke-width": "2",
                    "stroke-dasharray": "10",
                },
            )
        )


def add_snp_line(
    parent: ET.Element,
    record: BedRecord,
    y: float,
    x: float,
    width: float,
    stroke_width: float,
    dasharray: Optional[str],
    overhang_px: float,
    inset_px: float,
    halo_color: Optional[str],
    halo_width: float,
) -> None:
    line_start = x + inset_px - overhang_px
    line_end = x + width - inset_px + overhang_px
    if line_end < line_start:
        line_start = x + width / 2.0
        line_end = line_start
    line_attrs = {
        "x1": f"{line_start:.6f}",
        "y1": f"{y:.6f}",
        "x2": f"{line_end:.6f}",
        "y2": f"{y:.6f}",
        "stroke-linecap": "butt",
    }
    if halo_color and halo_width > 0:
        parent.append(
            ET.Element(
                "line",
                {
                    **line_attrs,
                    "stroke": halo_color,
                    "stroke-width": f"{halo_width:g}",
                },
            )
        )

    attrs = {
        **line_attrs,
        "stroke": record.color,
        "stroke-width": f"{stroke_width:g}",
    }
    if dasharray:
        attrs["stroke-dasharray"] = dasharray
    parent.append(ET.Element("line", attrs))


def add_snp_label(
    parent: ET.Element,
    label: str,
    y: float,
    x: float,
    width: float,
    chrom_y: float,
    font_size: float,
    offset_px: float,
    halo_color: Optional[str],
) -> None:
    if not label or font_size <= 0:
        return
    label_y = max(chrom_y + font_size + 1, y - offset_px)
    attrs = {
        "x": f"{x + width / 2:.6f}",
        "y": f"{label_y:.6f}",
        "font-size": f"{font_size:g}",
        "font-family": "Arial, Helvetica, sans-serif",
        "font-weight": "700",
        "text-anchor": "middle",
        "dominant-baseline": "alphabetic",
    }
    if halo_color:
        halo = ET.Element(
            "text",
            {
                **attrs,
                "fill": halo_color,
                "stroke": halo_color,
                "stroke-width": f"{SNP_LABEL_HALO_WIDTH_PX:g}",
                "stroke-linejoin": "round",
            },
        )
        halo.text = label
        parent.append(halo)
    text = ET.Element("text", {**attrs, "fill": "#111111"})
    text.text = label
    parent.append(text)


def build_snp_fallback_labels(records: Sequence[BedRecord]) -> Dict[Tuple[str, int, str], str]:
    labels: Dict[Tuple[str, int, str], str] = {}
    next_idx = 1
    for record in sorted(records, key=lambda r: r.source_line):
        if record.line_type != "geom_snp":
            continue
        variant_id = snp_variant_id(record)
        key = (record.chrom, record.start, variant_id)
        if key in labels:
            continue
        labels[key] = f"SNP{next_idx}"
        next_idx += 1
    return labels


def snp_label_for_record(record: BedRecord, fallback_labels: Dict[Tuple[str, int, str], str]) -> str:
    if len(record.extras) >= 5 and record.extras[4]:
        return record.extras[4]
    variant_id = snp_variant_id(record)
    return fallback_labels.get((record.chrom, record.start, variant_id), "")


def snp_haplotype_label(record: BedRecord, base_label: str) -> str:
    if not base_label:
        return ""
    return f"{base_label}/AN{record.haplotype}"


def snp_variant_id(record: BedRecord) -> str:
    return record.extras[0] if record.extras else f"chr{record.chrom}:{record.start}"


def record_geometry(record: BedRecord) -> Optional[Tuple[Dict[str, float], float, float]]:
    chromosome_id = f"chromosome{record.chrom}_hap{record.haplotype}"
    if chromosome_id not in CHROMOSOME_COORDINATES:
        LOGGER.warning("No template coordinates for %s; skipping line %d", chromosome_id, record.source_line)
        return None

    data = CHROMOSOME_COORDINATES[chromosome_id]
    chrom_len = CHROMOSOME_LENGTHS[record.chrom]
    clipped_start = max(0, min(record.start, chrom_len))
    clipped_end = max(0, min(record.end, chrom_len))
    if record.line_type in LINE_GEOMS:
        if clipped_end < clipped_start or record.end < 0 or record.start > chrom_len:
            LOGGER.warning("Interval outside chromosome bounds at BED line %d; skipping", record.source_line)
            return None
    elif clipped_end <= clipped_start:
        LOGGER.warning("Interval outside chromosome bounds at BED line %d; skipping", record.source_line)
        return None

    start_y = data["y"] + (clipped_start / chrom_len) * data["height"]
    end_y = data["y"] + (clipped_end / chrom_len) * data["height"]
    return data, start_y, end_y


def interval_geometry(chrom: str, haplotype: int, start_bp: int, end_bp: int) -> Optional[Tuple[Dict[str, float], float, float]]:
    chromosome_id = f"chromosome{chrom}_hap{haplotype}"
    if chromosome_id not in CHROMOSOME_COORDINATES or chrom not in CHROMOSOME_LENGTHS:
        return None
    data = CHROMOSOME_COORDINATES[chromosome_id]
    chrom_len = CHROMOSOME_LENGTHS[chrom]
    clipped_start = max(0, min(start_bp, chrom_len))
    clipped_end = max(0, min(end_bp, chrom_len))
    if clipped_end <= clipped_start:
        return None
    start_y = data["y"] + (clipped_start / chrom_len) * data["height"]
    end_y = data["y"] + (clipped_end / chrom_len) * data["height"]
    return data, start_y, end_y


def chromosomes_with_ancestry(records: Sequence[BedRecord]) -> List[str]:
    chroms = sorted({r.chrom for r in records if r.line_type == "geom_rect"}, key=chrom_sort_key)
    return chroms or list(CHROMOSOME_LENGTHS.keys())


def format_bp(value: int) -> str:
    return f"{int(value):,}"


def collect_snp_groups(records: Sequence[BedRecord], fallback_labels: Dict[Tuple[str, int, str], str]) -> List[SnpGroup]:
    grouped: "OrderedDict[Tuple[str, int, str], List[BedRecord]]" = OrderedDict()
    for record in sorted(records, key=lambda r: r.source_line):
        if record.line_type != "geom_snp":
            continue
        key = (record.chrom, record.start, snp_variant_id(record))
        grouped.setdefault(key, []).append(record)

    groups: List[SnpGroup] = []
    for (chrom, pos, variant_id), group_records in grouped.items():
        label = ""
        for record in group_records:
            label = snp_label_for_record(record, fallback_labels)
            if label:
                break
        groups.append(
            SnpGroup(
                chrom=chrom,
                pos=pos,
                variant_id=variant_id,
                label=label or f"chr{chrom}:{pos}",
                source_line=min(record.source_line for record in group_records),
                records=tuple(sorted(group_records, key=lambda r: (r.haplotype, r.source_line))),
            )
        )
    return sorted(groups, key=lambda g: (chrom_sort_key(g.chrom), g.pos, g.source_line))


def snp_group_template_y(group: SnpGroup) -> Optional[float]:
    chromosome_id = f"chromosome{group.chrom}_hap1"
    if chromosome_id not in CHROMOSOME_COORDINATES or group.chrom not in CHROMOSOME_LENGTHS:
        return None
    data = CHROMOSOME_COORDINATES[chromosome_id]
    chrom_len = CHROMOSOME_LENGTHS[group.chrom]
    pos = max(0, min(group.pos, chrom_len))
    return data["y"] + (pos / chrom_len) * data["height"]


def detect_snp_overlap_clusters(
    records: Sequence[BedRecord],
    fallback_labels: Dict[Tuple[str, int, str], str],
    min_px: float,
    min_snps: int,
) -> List[SnpOverlapCluster]:
    """Find SNP groups whose projected labels are too close in the main karyotype."""
    min_snps = max(2, min_snps)
    by_chrom: Dict[str, List[Tuple[SnpGroup, float]]] = {}
    for group in collect_snp_groups(records, fallback_labels):
        y_value = snp_group_template_y(group)
        if y_value is None:
            continue
        by_chrom.setdefault(group.chrom, []).append((group, y_value))

    clusters: List[SnpOverlapCluster] = []
    for chrom, group_pairs in sorted(by_chrom.items(), key=lambda item: chrom_sort_key(item[0])):
        current: List[Tuple[SnpGroup, float]] = []
        last_y: Optional[float] = None

        def flush_current() -> None:
            if len(current) < min_snps:
                return
            groups = tuple(group for group, _y in current)
            clusters.append(
                SnpOverlapCluster(
                    chrom=chrom,
                    start=min(group.pos for group in groups),
                    end=max(group.pos for group in groups),
                    groups=groups,
                )
            )

        for group, y_value in sorted(group_pairs, key=lambda item: (item[1], item[0].source_line)):
            if not current:
                current = [(group, y_value)]
            elif last_y is not None and abs(y_value - last_y) <= min_px:
                current.append((group, y_value))
            else:
                flush_current()
                current = [(group, y_value)]
            last_y = y_value
        flush_current()
    return clusters


def zoom_window_for_cluster(cluster: SnpOverlapCluster, padding_bp: int, min_window_bp: int) -> Tuple[int, int]:
    chrom_len = CHROMOSOME_LENGTHS[cluster.chrom]
    start = max(0, cluster.start - max(0, padding_bp))
    end = min(chrom_len, cluster.end + max(0, padding_bp))
    min_window = max(1, min(max(1, min_window_bp), chrom_len))
    if end - start < min_window:
        center = (cluster.start + cluster.end) / 2.0
        start = int(round(center - min_window / 2.0))
        end = start + min_window
        if start < 0:
            end -= start
            start = 0
        if end > chrom_len:
            start = max(0, start - (end - chrom_len))
            end = chrom_len
    if end <= start:
        end = min(chrom_len, start + 1)
    return start, end


def zoom_y_for_position(pos: int, window_start: int, window_end: int) -> float:
    if window_end <= window_start:
        return ZOOM_TRACK_TOP + ZOOM_TRACK_HEIGHT / 2.0
    fraction = (pos - window_start) / (window_end - window_start)
    return ZOOM_TRACK_TOP + max(0.0, min(1.0, fraction)) * ZOOM_TRACK_HEIGHT


def overlaps_window(record: BedRecord, window_start: int, window_end: int) -> bool:
    if record.line_type in LINE_GEOMS:
        return window_start <= record.start <= window_end
    return record.end > window_start and record.start < window_end


def clipped_to_window(start: int, end: int, window_start: int, window_end: int) -> Optional[Tuple[int, int]]:
    clipped_start = max(start, window_start)
    clipped_end = min(end, window_end)
    if clipped_end <= clipped_start:
        return None
    return clipped_start, clipped_end


def resolve_zoom_label_positions(
    entries: Sequence[Dict[str, object]],
    min_y: float,
    max_y: float,
    min_separation: float,
) -> List[Dict[str, object]]:
    if not entries:
        return []
    positioned = [dict(entry) for entry in sorted(entries, key=lambda item: (float(item["desired_y"]), int(item["source_line"])))]
    previous_y: Optional[float] = None
    for entry in positioned:
        desired = float(entry["desired_y"])
        label_y = max(min_y, desired if previous_y is None else max(desired, previous_y + min_separation))
        entry["label_y"] = label_y
        previous_y = label_y

    overflow = float(positioned[-1]["label_y"]) - max_y
    if overflow > 0:
        for entry in positioned:
            entry["label_y"] = float(entry["label_y"]) - overflow
    underflow = min_y - float(positioned[0]["label_y"])
    if underflow > 0:
        for entry in positioned:
            entry["label_y"] = float(entry["label_y"]) + underflow
    return positioned


def add_svg_centromere_overlays(
    parent: ET.Element,
    records: Sequence[BedRecord],
    color: str,
    opacity: float,
    overlap_px: float,
) -> None:
    for chrom in chromosomes_with_ancestry(records):
        if chrom not in CENTROMERE_BP:
            continue
        start_bp, end_bp = CENTROMERE_BP[chrom]
        for haplotype in (1, 2):
            geom = interval_geometry(chrom, haplotype, start_bp, end_bp)
            if geom is None:
                continue
            data, start_y, end_y = geom
            top = max(data["y"], start_y - overlap_px)
            bottom = min(data["y"] + data["height"], end_y + overlap_px)
            if bottom > top:
                add_segmented_rect(
                    parent,
                    data["x"],
                    top,
                    data["width"],
                    bottom - top,
                    color,
                    element_id=f"lap_centromere_chr{chrom}_hap{haplotype}",
                    opacity=opacity,
                )


def native_track_data(chromosome: str, haplotype: int) -> Dict[str, float]:
    return dict(CHROMOSOME_COORDINATES[f"chromosome{chromosome}_hap{haplotype}"])


def centromere_fraction(chromosome: str) -> float:
    interval = CENTROMERE_BP.get(chromosome)
    chrom_len = CHROMOSOME_LENGTHS.get(chromosome)
    if not interval or not chrom_len:
        return 0.5
    start, end = interval
    return max(0.06, min(0.94, ((start + end) / 2.0) / chrom_len))


def native_chromosome_path(chromosome: str, data: Dict[str, float], mpl_path_cls):
    """Build the original clean chromosome silhouette for the native PDF backend."""
    x = data["x"]
    y = data["y"]
    width = data["width"]
    height = data["height"]
    cx = x + width / 2.0
    top = y
    bottom = y + height
    cap = min(width * 0.48, height * 0.045)
    centromere_y = y + height * centromere_fraction(chromosome)
    centromere_y = max(top + cap * 1.8, min(bottom - cap * 1.8, centromere_y))

    waist_width = width * (0.38 if chromosome in ACROCENTRIC_CHROMS else 0.48)
    waist_half = waist_width / 2.0
    waist_height = max(width * (1.10 if chromosome in ACROCENTRIC_CHROMS else 1.35), 12.0)
    left = x
    right = x + width
    cy = centromere_y
    wh = waist_height / 2.0

    verts = [
        (cx, top),
        (left, top), (left, top), (left, top + cap),
        (left, cy - wh),
        (left, cy - wh * 0.45), (cx - waist_half, cy - wh * 0.35), (cx - waist_half, cy),
        (cx - waist_half, cy + wh * 0.35), (left, cy + wh * 0.45), (left, cy + wh),
        (left, bottom - cap),
        (left, bottom), (left, bottom), (cx, bottom),
        (right, bottom), (right, bottom), (right, bottom - cap),
        (right, cy + wh),
        (right, cy + wh * 0.45), (cx + waist_half, cy + wh * 0.35), (cx + waist_half, cy),
        (cx + waist_half, cy - wh * 0.35), (right, cy - wh * 0.45), (right, cy - wh),
        (right, top + cap),
        (right, top), (right, top), (cx, top),
        (cx, top),
    ]
    codes = [
        mpl_path_cls.MOVETO,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.LINETO,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.LINETO,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.LINETO,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.LINETO,
        mpl_path_cls.CURVE4, mpl_path_cls.CURVE4, mpl_path_cls.CURVE4,
        mpl_path_cls.CLOSEPOLY,
    ]
    return mpl_path_cls(verts, codes)


def svg_px_to_points(value: float) -> float:
    return value * PDF_POINTS_PER_SVG_UNIT


def parse_native_dasharray(dasharray: Optional[str]) -> Optional[List[float]]:
    if not dasharray:
        return None
    values: List[float] = []
    for token in re.split(r"[,\s]+", dasharray.strip()):
        if not token:
            continue
        try:
            values.append(svg_px_to_points(float(token)))
        except ValueError:
            LOGGER.warning("Ignoring malformed dash token in '%s': %s", dasharray, token)
    return values or None


def build_legend_entries(records: Sequence[BedRecord], palette: "OrderedDict[str, str]") -> List[Tuple[str, str]]:
    found_colors_ordered: List[str] = []
    seen_colors = set()
    for record in records:
        if record.line_type not in LEGEND_GEOMS:
            continue
        color_key = record.color.lower()
        if color_key not in seen_colors:
            found_colors_ordered.append(record.color)
            seen_colors.add(color_key)

    found_color_set = {c.lower() for c in found_colors_ordered}
    entries: List[Tuple[str, str]] = []
    used_colors = set()
    for name, color in palette.items():
        if color.lower() in found_color_set and color.lower() not in used_colors:
            entries.append((name, color))
            used_colors.add(color.lower())

    for color in found_colors_ordered:
        if color.lower() not in used_colors:
            entries.append((color, color))
    return entries


def chromosome_template_x_bounds() -> Tuple[float, float]:
    left = min(data["x"] for data in CHROMOSOME_COORDINATES.values())
    right = max(data["x"] + data["width"] for data in CHROMOSOME_COORDINATES.values())
    return left, right


def estimate_native_text_width(text: str, font_size_px: float) -> float:
    return len(text) * font_size_px * NATIVE_LEGEND_TEXT_WIDTH_FACTOR


def native_canvas_x_limits(
    legend_entries: Sequence[Tuple[str, str]],
    include_centromere_legend: bool,
    include_acrocentric_p_arm_legend: bool,
) -> Tuple[float, float]:
    chrom_left, chrom_right = chromosome_template_x_bounds()
    legend_labels = [name for name, _color in legend_entries]
    if include_centromere_legend:
        legend_labels.append(CENTROMERE_LEGEND_LABEL)
    if include_acrocentric_p_arm_legend:
        legend_labels.append(ACROCENTRIC_P_ARM_LEGEND_LABEL)
    legend_right = LEGEND_TEXT_X
    if legend_labels:
        legend_right += max(estimate_native_text_width(label, LEGEND_FONT_SIZE) for label in legend_labels)

    content_right = max(chrom_right, legend_right)
    x_min = chrom_left - NATIVE_CANVAS_LEFT_PADDING_PX
    x_max = x_min + SVG_VIEWBOX_WIDTH
    required_x_max = content_right + NATIVE_CANVAS_RIGHT_PADDING_PX
    if x_max < required_x_max:
        x_max = required_x_max
    return x_min, x_max


def write_native_pdf(
    bed_file: str,
    output_pdf_file: str,
    individual_name: str,
    color_config: Optional[str],
    overlap_px: float,
    missing_color: Optional[str],
    strict: bool,
    snp_line_width: float,
    snp_line_dasharray: Optional[str],
    snp_line_overhang_px: float,
    snp_line_inset_px: float,
    snp_line_halo_color: Optional[str],
    snp_line_halo_width: float,
    snp_labels: bool,
    snp_label_font_size: float,
    snp_label_offset_px: float,
    snp_label_halo_color: Optional[str],
    centromere_overlay: bool,
    centromere_overlay_color: str,
    centromere_overlay_opacity: float,
) -> None:
    """Render a cleaner direct PDF without using the SVG template or rsvg-convert."""
    try:
        import matplotlib

        matplotlib.use("pdf")
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        logging.getLogger("fontTools").setLevel(logging.ERROR)
        import matplotlib.patheffects as path_effects
        import matplotlib.pyplot as plt
        from matplotlib.patches import PathPatch, Rectangle
        from matplotlib.path import Path as MplPath
    except ImportError as exc:
        raise RuntimeError("Native PDF output requires matplotlib. Install matplotlib or use --pdf-backend rsvg.") from exc

    bed_path = Path(bed_file)
    headers = read_header_lines(bed_path)
    palette = build_ancestry_palette(headers, color_config)
    records = parse_bed_records(bed_path, strict=strict)
    log_coverage_summary(records)
    snp_fallback_labels = build_snp_fallback_labels(records)
    dash_pattern = parse_native_dasharray(snp_line_dasharray)
    legend_entries = build_legend_entries(records, palette)
    include_acrocentric_p_arm_legend = True
    native_x_min, native_x_max = native_canvas_x_limits(
        legend_entries,
        centromere_overlay,
        include_acrocentric_p_arm_legend,
    )
    native_width = native_x_max - native_x_min

    fig = plt.figure(figsize=(native_width / 100.0, SVG_VIEWBOX_HEIGHT / 100.0), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(native_x_min, native_x_max)
    ax.set_ylim(SVG_VIEWBOX_HEIGHT, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    chromosome_paths = {}
    for chrom in CHROMOSOME_LENGTHS:
        for haplotype in (1, 2):
            data = native_track_data(chrom, haplotype)
            path = native_chromosome_path(chrom, data, MplPath)
            chromosome_paths[(chrom, haplotype)] = path
            ax.add_patch(PathPatch(path, facecolor="#ffffff", edgecolor="none", zorder=0))

    def add_clipped_rect(
        chrom: str,
        haplotype: int,
        x: float,
        y: float,
        width: float,
        height: float,
        color: str,
        zorder: float,
        alpha: float = 1.0,
    ) -> None:
        path = chromosome_paths.get((chrom, haplotype))
        if path is None:
            return
        rect = Rectangle(
            (x, y),
            width,
            height,
            facecolor=color,
            edgecolor="none",
            linewidth=0,
            antialiased=False,
            alpha=max(0.0, min(1.0, alpha)),
            zorder=zorder,
        )
        rect.set_clip_path(path, ax.transData)
        ax.add_patch(rect)

    if missing_color:
        chroms_with_data = sorted({r.chrom for r in records if r.line_type == "geom_rect"}, key=chrom_sort_key)
        for chrom in chroms_with_data:
            for haplotype in (1, 2):
                data = native_track_data(chrom, haplotype)
                add_clipped_rect(chrom, haplotype, data["x"], data["y"], data["width"], data["height"], missing_color, 0.5)

    snp_label_regions: Dict[Tuple[str, int, str, str, int], Dict[str, object]] = {}
    for record in records:
        geom = record_geometry(record)
        if geom is None:
            continue
        _template_data, start_y, end_y = geom
        chromosome_data = native_track_data(record.chrom, record.haplotype)
        x = chromosome_data["x"]
        width = chromosome_data["width"]

        if record.line_type == "geom_line":
            offset = 15
            for y_value in (start_y, end_y):
                ax.plot(
                    [x - offset, x + width + offset],
                    [y_value, y_value],
                    color=record.color,
                    linewidth=svg_px_to_points(2.0),
                    dashes=parse_native_dasharray("10 10"),
                    solid_capstyle="butt",
                    zorder=4,
                )
            continue

        if record.line_type == "geom_snp":
            line_start = x + snp_line_inset_px - snp_line_overhang_px
            line_end = x + width - snp_line_inset_px + snp_line_overhang_px
            if line_end < line_start:
                line_start = x + width / 2.0
                line_end = line_start
            line_x = [line_start, line_end]
            if snp_line_halo_color and snp_line_halo_width > 0:
                ax.plot(
                    line_x,
                    [start_y, start_y],
                    color=snp_line_halo_color,
                    linewidth=svg_px_to_points(snp_line_halo_width),
                    solid_capstyle="butt",
                    zorder=5,
                )
            line, = ax.plot(
                line_x,
                [start_y, start_y],
                color=record.color,
                linewidth=svg_px_to_points(snp_line_width),
                solid_capstyle="butt",
                zorder=6,
            )
            if dash_pattern:
                line.set_dashes(dash_pattern)

            if snp_labels:
                label = snp_haplotype_label(record, snp_label_for_record(record, snp_fallback_labels))
                variant_id = snp_variant_id(record)
                key = (record.chrom, record.start, variant_id, label, record.haplotype)
                entry = snp_label_regions.setdefault(
                    key,
                    {
                        "label": label,
                        "y": start_y,
                        "x_min": x,
                        "x_max": x + width,
                        "chrom_y": chromosome_data["y"],
                        "source_line": record.source_line,
                    },
                )
                entry["x_min"] = min(float(entry["x_min"]), x)
                entry["x_max"] = max(float(entry["x_max"]), x + width)
                entry["y"] = min(float(entry["y"]), start_y)
                entry["source_line"] = min(int(entry["source_line"]), record.source_line)
            continue

        top = max(chromosome_data["y"], start_y - overlap_px)
        bottom = min(chromosome_data["y"] + chromosome_data["height"], end_y + overlap_px)
        if bottom > top:
            add_clipped_rect(record.chrom, record.haplotype, x, top, width, bottom - top, record.color, 2)

    for chrom in sorted(ACROCENTRIC_CHROMS, key=chrom_sort_key):
        if chrom not in CENTROMERE_BP:
            continue
        cent_start, cent_end = CENTROMERE_BP[chrom]
        p_arm_end_bp = int(round((cent_start + cent_end) / 2.0))
        for haplotype in (1, 2):
            geom = interval_geometry(chrom, haplotype, 0, p_arm_end_bp)
            if geom is None:
                continue
            data, start_y, end_y = geom
            top = max(data["y"], start_y)
            bottom = min(data["y"] + data["height"], end_y)
            if bottom > top:
                add_clipped_rect(
                    chrom,
                    haplotype,
                    data["x"],
                    top,
                    data["width"],
                    bottom - top,
                    ACROCENTRIC_P_ARM_COLOR,
                    3.2,
                )

    if centromere_overlay:
        for chrom in chromosomes_with_ancestry(records):
            if chrom not in CENTROMERE_BP:
                continue
            start_bp, end_bp = CENTROMERE_BP[chrom]
            for haplotype in (1, 2):
                geom = interval_geometry(chrom, haplotype, start_bp, end_bp)
                if geom is None:
                    continue
                data, start_y, end_y = geom
                top = max(data["y"], start_y - overlap_px)
                bottom = min(data["y"] + data["height"], end_y + overlap_px)
                if bottom > top:
                    for seg_x, seg_width in segmented_spans(data["x"], data["width"]):
                        add_clipped_rect(
                            chrom,
                            haplotype,
                            seg_x,
                            top,
                            seg_width,
                            bottom - top,
                            centromere_overlay_color,
                            3,
                            centromere_overlay_opacity,
                        )

    for chrom in CHROMOSOME_LENGTHS:
        for haplotype in (1, 2):
            path = chromosome_paths.get((chrom, haplotype))
            if path is None:
                continue
            ax.add_patch(
                PathPatch(
                    path,
                    facecolor="none",
                    edgecolor="#000000",
                    linewidth=svg_px_to_points(0.9),
                    joinstyle="round",
                    capstyle="round",
                    zorder=7,
                )
            )

    for chrom in CHROMOSOME_LENGTHS:
        hap1 = native_track_data(chrom, 1)
        hap2 = native_track_data(chrom, 2)
        x_center = (hap1["x"] + hap2["x"] + hap2["width"]) / 2.0
        ax.text(
            x_center,
            792,
            f"CHR{chrom}",
            color="#000000",
            fontsize=svg_px_to_points(15),
            fontweight="bold",
            ha="center",
            va="baseline",
            family="Helvetica",
            zorder=8,
        )

    for entry in sorted(snp_label_regions.values(), key=lambda item: int(item["source_line"])):
        label = str(entry["label"])
        if not label or snp_label_font_size <= 0:
            continue
        x_min = float(entry["x_min"])
        x_max = float(entry["x_max"])
        y_value = max(float(entry["chrom_y"]) + snp_label_font_size + 1, float(entry["y"]) - snp_label_offset_px)
        text = ax.text(
            (x_min + x_max) / 2.0,
            y_value,
            label,
            color="#111111",
            fontsize=svg_px_to_points(snp_label_font_size),
            fontweight="bold",
            ha="center",
            va="baseline",
            family="Helvetica",
            zorder=9,
        )
        if snp_label_halo_color:
            text.set_path_effects(
                [
                    path_effects.withStroke(
                        linewidth=svg_px_to_points(SNP_LABEL_HALO_WIDTH_PX),
                        foreground=snp_label_halo_color,
                    )
                ]
            )

    ax.text(
        (native_x_min + native_x_max) / 2.0,
        30,
        individual_name,
        color="#000000",
        fontsize=svg_px_to_points(32),
        ha="center",
        va="baseline",
        zorder=10,
    )

    def add_native_centromere_legend_swatch(y: float) -> None:
        ax.add_patch(
            Rectangle(
                (LEGEND_SWATCH_X, y),
                LEGEND_SWATCH_WIDTH,
                LEGEND_SWATCH_HEIGHT,
                facecolor="#ffffff",
                edgecolor="#666666",
                linewidth=svg_px_to_points(0.7),
                zorder=10,
            )
        )
        for seg_x, seg_width in segmented_spans(LEGEND_SWATCH_X, LEGEND_SWATCH_WIDTH):
            ax.add_patch(
                Rectangle(
                    (seg_x, y),
                    seg_width,
                    LEGEND_SWATCH_HEIGHT,
                    facecolor=centromere_overlay_color,
                    edgecolor="none",
                    alpha=max(0.0, min(1.0, centromere_overlay_opacity)),
                    zorder=11,
                )
            )

    y_offset = 40
    for name, color in legend_entries:
        ax.add_patch(
            Rectangle(
                (LEGEND_SWATCH_X, y_offset),
                LEGEND_SWATCH_WIDTH,
                LEGEND_SWATCH_HEIGHT,
                facecolor=color,
                edgecolor="none",
                zorder=10,
            )
        )
        ax.text(
            LEGEND_TEXT_X,
            y_offset + 14,
            name,
            color="#000000",
            fontsize=svg_px_to_points(LEGEND_FONT_SIZE),
            ha="left",
            va="baseline",
            zorder=10,
        )
        y_offset += LEGEND_ROW_STEP
    if centromere_overlay:
        add_native_centromere_legend_swatch(y_offset)
        ax.text(
            LEGEND_TEXT_X,
            y_offset + 14,
            CENTROMERE_LEGEND_LABEL,
            color="#000000",
            fontsize=svg_px_to_points(LEGEND_FONT_SIZE),
            ha="left",
            va="baseline",
            zorder=10,
        )
        y_offset += LEGEND_ROW_STEP
    if include_acrocentric_p_arm_legend:
        ax.add_patch(
            Rectangle(
                (LEGEND_SWATCH_X, y_offset),
                LEGEND_SWATCH_WIDTH,
                LEGEND_SWATCH_HEIGHT,
                facecolor=ACROCENTRIC_P_ARM_COLOR,
                edgecolor="none",
                zorder=10,
            )
        )
        ax.text(
            LEGEND_TEXT_X,
            y_offset + 14,
            ACROCENTRIC_P_ARM_LEGEND_LABEL,
            color="#000000",
            fontsize=svg_px_to_points(LEGEND_FONT_SIZE),
            ha="left",
            va="baseline",
            zorder=10,
        )

    os.makedirs(os.path.dirname(output_pdf_file) or ".", exist_ok=True)
    fig.savefig(output_pdf_file, format="pdf", facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close(fig)
    LOGGER.info("Native PDF written: %s", output_pdf_file)


def write_native_zoom_panel(
    output_pdf_file: str,
    individual_name: str,
    records: Sequence[BedRecord],
    palette: "OrderedDict[str, str]",
    snp_fallback_labels: Dict[Tuple[str, int, str], str],
    cluster: SnpOverlapCluster,
    window_start: int,
    window_end: int,
    missing_color: Optional[str],
    snp_line_width: float,
    snp_line_dasharray: Optional[str],
    snp_line_overhang_px: float,
    snp_line_inset_px: float,
    snp_line_halo_color: Optional[str],
    snp_line_halo_width: float,
    snp_labels: bool,
    snp_label_font_size: float,
    snp_label_halo_color: Optional[str],
    centromere_overlay: bool,
    centromere_overlay_color: str,
    centromere_overlay_opacity: float,
) -> None:
    try:
        import matplotlib

        matplotlib.use("pdf")
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        logging.getLogger("fontTools").setLevel(logging.ERROR)
        import matplotlib.patheffects as path_effects
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError as exc:
        raise RuntimeError("Zoom PDF output requires matplotlib. Install matplotlib before using --auto-zoom-overlaps.") from exc

    dash_pattern = parse_native_dasharray(snp_line_dasharray)
    chrom = cluster.chrom
    chrom_len = CHROMOSOME_LENGTHS[chrom]
    panel_records = [record for record in records if record.chrom == chrom and overlaps_window(record, window_start, window_end)]
    legend_entries = build_legend_entries(panel_records, palette)
    track_x = {1: ZOOM_HAP1_X, 2: ZOOM_HAP2_X}
    zoom_label_font = max(snp_label_font_size + 2.0, 8.0)

    fig = plt.figure(figsize=(ZOOM_PANEL_WIDTH / 100.0, ZOOM_PANEL_HEIGHT / 100.0), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, ZOOM_PANEL_WIDTH)
    ax.set_ylim(ZOOM_PANEL_HEIGHT, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    region_label = f"chr{chrom}:{format_bp(window_start)}-{format_bp(window_end)}"
    title = f"{individual_name} - zoom {region_label}"
    ax.text(
        ZOOM_PANEL_WIDTH / 2.0,
        34,
        title,
        color="#000000",
        fontsize=svg_px_to_points(24),
        ha="center",
        va="baseline",
        family="Helvetica",
        zorder=20,
    )
    ax.text(
        ZOOM_PANEL_WIDTH / 2.0,
        60,
        f"{len(cluster.groups)} clustered SNPs: " + ", ".join(group.label for group in cluster.groups),
        color="#333333",
        fontsize=svg_px_to_points(13),
        ha="center",
        va="baseline",
        family="Helvetica",
        zorder=20,
    )

    ax.add_patch(
        Rectangle(
            (ZOOM_CONTEXT_X, ZOOM_CONTEXT_Y),
            ZOOM_CONTEXT_WIDTH,
            ZOOM_CONTEXT_HEIGHT,
            facecolor="#f7f7f7",
            edgecolor="#222222",
            linewidth=svg_px_to_points(0.9),
            zorder=1,
        )
    )
    context_start_y = ZOOM_CONTEXT_Y + (window_start / chrom_len) * ZOOM_CONTEXT_HEIGHT
    context_end_y = ZOOM_CONTEXT_Y + (window_end / chrom_len) * ZOOM_CONTEXT_HEIGHT
    ax.add_patch(
        Rectangle(
            (ZOOM_CONTEXT_X - 5, context_start_y),
            ZOOM_CONTEXT_WIDTH + 10,
            max(3.0, context_end_y - context_start_y),
            facecolor="#f3c96b",
            edgecolor="#c97800",
            linewidth=svg_px_to_points(0.7),
            alpha=0.68,
            zorder=2,
        )
    )
    for group in cluster.groups:
        y_value = ZOOM_CONTEXT_Y + (max(0, min(group.pos, chrom_len)) / chrom_len) * ZOOM_CONTEXT_HEIGHT
        ax.plot(
            [ZOOM_CONTEXT_X - 7, ZOOM_CONTEXT_X + ZOOM_CONTEXT_WIDTH + 7],
            [y_value, y_value],
            color="#111111",
            linewidth=svg_px_to_points(1.0),
            solid_capstyle="butt",
            zorder=3,
        )
    ax.text(
        ZOOM_CONTEXT_X + ZOOM_CONTEXT_WIDTH / 2.0,
        ZOOM_CONTEXT_Y - 14,
        f"chr{chrom}",
        color="#000000",
        fontsize=svg_px_to_points(14),
        fontweight="bold",
        ha="center",
        va="baseline",
        family="Helvetica",
        zorder=4,
    )

    ax.plot(
        [ZOOM_AXIS_X, ZOOM_AXIS_X],
        [ZOOM_TRACK_TOP, ZOOM_TRACK_TOP + ZOOM_TRACK_HEIGHT],
        color="#333333",
        linewidth=svg_px_to_points(0.8),
        zorder=2,
    )
    tick_values = [window_start, int(round((window_start + window_end) / 2.0)), window_end]
    for tick in tick_values:
        y_value = zoom_y_for_position(tick, window_start, window_end)
        ax.plot([ZOOM_AXIS_X - 6, ZOOM_AXIS_X], [y_value, y_value], color="#333333", linewidth=svg_px_to_points(0.8), zorder=2)
        ax.text(
            ZOOM_AXIS_X - 10,
            y_value + 4,
            format_bp(tick),
            color="#333333",
            fontsize=svg_px_to_points(10),
            ha="right",
            va="baseline",
            family="Helvetica",
            zorder=2,
        )

    for haplotype, x in track_x.items():
        ax.add_patch(
            Rectangle(
                (x, ZOOM_TRACK_TOP),
                ZOOM_TRACK_WIDTH,
                ZOOM_TRACK_HEIGHT,
                facecolor=missing_color or "#ffffff",
                edgecolor="none",
                zorder=0,
            )
        )
        ax.text(
            x + ZOOM_TRACK_WIDTH / 2.0,
            ZOOM_TRACK_TOP - 18,
            f"hap{haplotype}",
            color="#000000",
            fontsize=svg_px_to_points(14),
            fontweight="bold",
            ha="center",
            va="baseline",
            family="Helvetica",
            zorder=8,
        )

    for record in panel_records:
        if record.line_type != "geom_rect":
            continue
        clipped = clipped_to_window(record.start, record.end, window_start, window_end)
        if clipped is None:
            continue
        y1 = zoom_y_for_position(clipped[0], window_start, window_end)
        y2 = zoom_y_for_position(clipped[1], window_start, window_end)
        ax.add_patch(
            Rectangle(
                (track_x[record.haplotype], y1),
                ZOOM_TRACK_WIDTH,
                max(0.8, y2 - y1),
                facecolor=record.color,
                edgecolor="none",
                antialiased=False,
                zorder=2,
            )
        )

    if chrom in ACROCENTRIC_CHROMS and chrom in CENTROMERE_BP:
        cent_start, cent_end = CENTROMERE_BP[chrom]
        p_arm_end_bp = int(round((cent_start + cent_end) / 2.0))
        clipped = clipped_to_window(0, p_arm_end_bp, window_start, window_end)
        if clipped is not None:
            y1 = zoom_y_for_position(clipped[0], window_start, window_end)
            y2 = zoom_y_for_position(clipped[1], window_start, window_end)
            for x in track_x.values():
                ax.add_patch(
                    Rectangle(
                        (x, y1),
                        ZOOM_TRACK_WIDTH,
                        y2 - y1,
                        facecolor=ACROCENTRIC_P_ARM_COLOR,
                        edgecolor="none",
                        zorder=3.2,
                    )
                )

    centromere_visible = False
    if centromere_overlay and chrom in CENTROMERE_BP:
        clipped = clipped_to_window(CENTROMERE_BP[chrom][0], CENTROMERE_BP[chrom][1], window_start, window_end)
        if clipped is not None:
            centromere_visible = True
            y1 = zoom_y_for_position(clipped[0], window_start, window_end)
            y2 = zoom_y_for_position(clipped[1], window_start, window_end)
            for x in track_x.values():
                for seg_x, seg_width in segmented_spans(x, ZOOM_TRACK_WIDTH):
                    ax.add_patch(
                        Rectangle(
                            (seg_x, y1),
                            seg_width,
                            y2 - y1,
                            facecolor=centromere_overlay_color,
                            edgecolor="none",
                            alpha=max(0.0, min(1.0, centromere_overlay_opacity)),
                            zorder=4,
                        )
                    )

    label_entries: Dict[int, List[Dict[str, object]]] = {1: [], 2: []}
    for record in panel_records:
        if record.line_type != "geom_snp":
            continue
        x = track_x[record.haplotype]
        y_value = zoom_y_for_position(record.start, window_start, window_end)
        line_start = x + snp_line_inset_px - snp_line_overhang_px
        line_end = x + ZOOM_TRACK_WIDTH - snp_line_inset_px + snp_line_overhang_px
        if snp_line_halo_color and snp_line_halo_width > 0:
            ax.plot(
                [line_start, line_end],
                [y_value, y_value],
                color=snp_line_halo_color,
                linewidth=svg_px_to_points(snp_line_halo_width),
                solid_capstyle="butt",
                zorder=6,
            )
        line, = ax.plot(
            [line_start, line_end],
            [y_value, y_value],
            color=record.color,
            linewidth=svg_px_to_points(snp_line_width),
            solid_capstyle="butt",
            zorder=7,
        )
        if dash_pattern:
            line.set_dashes(dash_pattern)

        if snp_labels:
            label = snp_haplotype_label(record, snp_label_for_record(record, snp_fallback_labels))
            if label:
                label_entries[record.haplotype].append(
                    {
                        "label": label,
                        "desired_y": y_value,
                        "line_y": y_value,
                        "source_line": record.source_line,
                        "record": record,
                    }
                )

    for haplotype, x in track_x.items():
        ax.add_patch(
            Rectangle(
                (x, ZOOM_TRACK_TOP),
                ZOOM_TRACK_WIDTH,
                ZOOM_TRACK_HEIGHT,
                facecolor="none",
                edgecolor="#000000",
                linewidth=svg_px_to_points(1.1),
                zorder=8,
            )
        )
        side_entries = resolve_zoom_label_positions(
            label_entries[haplotype],
            ZOOM_TRACK_TOP + zoom_label_font,
            ZOOM_TRACK_TOP + ZOOM_TRACK_HEIGHT - 3,
            ZOOM_LABEL_MIN_SEPARATION_PX,
        )
        for entry in side_entries:
            line_y = float(entry["line_y"])
            label_y = float(entry["label_y"])
            if haplotype == 1:
                text_x = x - ZOOM_LABEL_GAP
                anchor_x = x
                ha = "right"
                leader_text_x = text_x + 3
            else:
                text_x = x + ZOOM_TRACK_WIDTH + ZOOM_LABEL_GAP
                anchor_x = x + ZOOM_TRACK_WIDTH
                ha = "left"
                leader_text_x = text_x - 3
            ax.plot(
                [anchor_x, leader_text_x],
                [line_y, label_y],
                color="#555555",
                linewidth=svg_px_to_points(0.6),
                zorder=8.5,
            )
            text = ax.text(
                text_x,
                label_y + zoom_label_font * 0.35,
                str(entry["label"]),
                color="#111111",
                fontsize=svg_px_to_points(zoom_label_font),
                fontweight="bold",
                ha=ha,
                va="baseline",
                family="Helvetica",
                zorder=9,
            )
            if snp_label_halo_color:
                text.set_path_effects(
                    [
                        path_effects.withStroke(
                            linewidth=svg_px_to_points(SNP_LABEL_HALO_WIDTH_PX),
                            foreground=snp_label_halo_color,
                        )
                    ]
                )

    legend_x = 618.0
    legend_y = 112.0
    ax.text(
        legend_x,
        legend_y - 22,
        "Legend",
        color="#000000",
        fontsize=svg_px_to_points(15),
        fontweight="bold",
        ha="left",
        va="baseline",
        family="Helvetica",
        zorder=20,
    )
    for name, color in legend_entries:
        ax.add_patch(Rectangle((legend_x, legend_y), 22, 14, facecolor=color, edgecolor="none", zorder=20))
        ax.text(
            legend_x + 30,
            legend_y + 12,
            name,
            color="#000000",
            fontsize=svg_px_to_points(12),
            ha="left",
            va="baseline",
            family="Helvetica",
            zorder=20,
        )
        legend_y += 24
    if centromere_visible:
        ax.add_patch(
            Rectangle(
                (legend_x, legend_y),
                22,
                14,
                facecolor="#ffffff",
                edgecolor="#666666",
                linewidth=svg_px_to_points(0.7),
                zorder=20,
            )
        )
        for seg_x, seg_width in segmented_spans(legend_x, 22):
            ax.add_patch(
                Rectangle(
                    (seg_x, legend_y),
                    seg_width,
                    14,
                    facecolor=centromere_overlay_color,
                    edgecolor="none",
                    alpha=max(0.0, min(1.0, centromere_overlay_opacity)),
                    zorder=21,
                )
            )
        ax.text(
            legend_x + 30,
            legend_y + 12,
            CENTROMERE_LEGEND_LABEL,
            color="#000000",
            fontsize=svg_px_to_points(12),
            ha="left",
            va="baseline",
            family="Helvetica",
            zorder=20,
        )
        legend_y += 24
    if chrom in ACROCENTRIC_CHROMS:
        ax.add_patch(Rectangle((legend_x, legend_y), 22, 14, facecolor=ACROCENTRIC_P_ARM_COLOR, edgecolor="none", zorder=20))
        ax.text(
            legend_x + 30,
            legend_y + 12,
            ACROCENTRIC_P_ARM_LEGEND_LABEL,
            color="#000000",
            fontsize=svg_px_to_points(12),
            ha="left",
            va="baseline",
            family="Helvetica",
            zorder=20,
        )

    os.makedirs(os.path.dirname(output_pdf_file) or ".", exist_ok=True)
    fig.savefig(output_pdf_file, format="pdf", facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close(fig)
    LOGGER.info("Native zoom PDF written: %s", output_pdf_file)


def write_native_zoom_overlap_plots(
    bed_file: str,
    output_dir: str,
    output_stem: str,
    individual_name: str,
    color_config: Optional[str],
    missing_color: Optional[str],
    strict: bool,
    snp_line_width: float,
    snp_line_dasharray: Optional[str],
    snp_line_overhang_px: float,
    snp_line_inset_px: float,
    snp_line_halo_color: Optional[str],
    snp_line_halo_width: float,
    snp_labels: bool,
    snp_label_font_size: float,
    snp_label_halo_color: Optional[str],
    centromere_overlay: bool,
    centromere_overlay_color: str,
    centromere_overlay_opacity: float,
    zoom_overlap_px: float,
    zoom_padding_bp: int,
    zoom_min_window_bp: int,
    zoom_min_snps: int,
    zoom_max_panels: int,
) -> List[str]:
    bed_path = Path(bed_file)
    headers = read_header_lines(bed_path)
    palette = build_ancestry_palette(headers, color_config)
    records = parse_bed_records(bed_path, strict=strict)
    snp_fallback_labels = build_snp_fallback_labels(records)
    clusters = detect_snp_overlap_clusters(records, snp_fallback_labels, zoom_overlap_px, zoom_min_snps)
    if not clusters:
        LOGGER.info("No overlapping SNP clusters detected; no zoom panels written.")
        return []

    os.makedirs(output_dir, exist_ok=True)
    written: List[str] = []
    for idx, cluster in enumerate(clusters[: max(0, zoom_max_panels)], start=1):
        window_start, window_end = zoom_window_for_cluster(cluster, zoom_padding_bp, zoom_min_window_bp)
        out_path = Path(output_dir) / f"{output_stem}.zoom{idx:02d}.chr{cluster.chrom}_{window_start}_{window_end}.pdf"
        write_native_zoom_panel(
            str(out_path),
            individual_name,
            records,
            palette,
            snp_fallback_labels,
            cluster,
            window_start,
            window_end,
            missing_color,
            snp_line_width,
            snp_line_dasharray,
            snp_line_overhang_px,
            snp_line_inset_px,
            snp_line_halo_color,
            snp_line_halo_width,
            snp_labels,
            snp_label_font_size,
            snp_label_halo_color,
            centromere_overlay,
            centromere_overlay_color,
            centromere_overlay_opacity,
        )
        written.append(str(out_path))
    if len(clusters) > len(written):
        LOGGER.warning(
            "Detected %d overlapping SNP clusters, wrote %d panels because --zoom-max-panels is %d.",
            len(clusters),
            len(written),
            zoom_max_panels,
        )
    return written


def insert_colored_regions(
    svg_file: str,
    bed_file: str,
    output_svg_file: str,
    output_pdf_file: str,
    individual_name: str,
    color_config: Optional[str],
    overlap_px: float,
    missing_color: Optional[str],
    strict: bool,
    snp_line_width: float,
    snp_line_dasharray: Optional[str],
    snp_line_overhang_px: float,
    snp_line_inset_px: float,
    snp_line_halo_color: Optional[str],
    snp_line_halo_width: float,
    snp_labels: bool,
    snp_label_font_size: float,
    snp_label_offset_px: float,
    snp_label_halo_color: Optional[str],
    centromere_overlay: bool,
    centromere_overlay_color: str,
    centromere_overlay_opacity: float,
    export_pdf: bool,
) -> None:
    bed_path = Path(bed_file)
    headers = read_header_lines(bed_path)
    palette = build_ancestry_palette(headers, color_config)
    records = parse_bed_records(bed_path, strict=strict)
    log_coverage_summary(records)
    snp_fallback_labels = build_snp_fallback_labels(records)

    original_svg_tree = ET.parse(svg_file)
    original_svg_root = original_svg_tree.getroot()
    show_acrocentric_p_arm_legend = svg_uses_acrocentric_p_arm_color(original_svg_root)

    region_group = ET.Element("g", {"id": "lap_colored_regions"})
    annotation_group = ET.Element("g", {"id": "lap_annotations"})

    if missing_color:
        chroms_with_data = sorted({r.chrom for r in records if r.line_type == "geom_rect"}, key=chrom_sort_key)
        for chrom in chroms_with_data:
            for haplotype in (1, 2):
                data = CHROMOSOME_COORDINATES[f"chromosome{chrom}_hap{haplotype}"]
                add_rect(region_group, data["x"], data["y"], data["width"], data["height"], missing_color)

    snp_label_regions: Dict[Tuple[str, int, str, str, int], Dict[str, object]] = {}
    for record in records:
        geom = record_geometry(record)
        if geom is None:
            continue
        chromosome_data, start_y, end_y = geom
        x = chromosome_data["x"]
        width = chromosome_data["width"]

        if record.line_type == "geom_line":
            add_feature_lines(annotation_group, record, start_y, end_y, x, width)
            continue
        if record.line_type == "geom_snp":
            add_snp_line(
                annotation_group,
                record,
                start_y,
                x,
                width,
                snp_line_width,
                snp_line_dasharray,
                snp_line_overhang_px,
                snp_line_inset_px,
                snp_line_halo_color,
                snp_line_halo_width,
            )
            if snp_labels:
                label = snp_haplotype_label(record, snp_label_for_record(record, snp_fallback_labels))
                variant_id = snp_variant_id(record)
                key = (record.chrom, record.start, variant_id, label, record.haplotype)
                entry = snp_label_regions.setdefault(
                    key,
                    {
                        "label": label,
                        "y": start_y,
                        "x_min": x,
                        "x_max": x + width,
                        "chrom_y": chromosome_data["y"],
                        "source_line": record.source_line,
                    },
                )
                entry["x_min"] = min(float(entry["x_min"]), x)
                entry["x_max"] = max(float(entry["x_max"]), x + width)
                entry["y"] = min(float(entry["y"]), start_y)
                entry["source_line"] = min(int(entry["source_line"]), record.source_line)
            continue

        top = max(chromosome_data["y"], start_y - overlap_px)
        bottom = min(chromosome_data["y"] + chromosome_data["height"], end_y + overlap_px)
        if bottom > top:
            add_rect(region_group, x, top, width, bottom - top, record.color)

    if centromere_overlay:
        add_svg_centromere_overlays(
            region_group,
            records,
            centromere_overlay_color,
            centromere_overlay_opacity,
            overlap_px,
        )

    for entry in sorted(snp_label_regions.values(), key=lambda item: int(item["source_line"])):
        x_min = float(entry["x_min"])
        x_max = float(entry["x_max"])
        add_snp_label(
            annotation_group,
            str(entry["label"]),
            float(entry["y"]),
            x_min,
            x_max - x_min,
            float(entry["chrom_y"]),
            snp_label_font_size,
            snp_label_offset_px,
            snp_label_halo_color,
        )

    title_element = ET.Element("text", {"x": "800", "y": "30", "fill": "black", "font-size": "32"})
    title_element.text = individual_name
    annotation_group.append(title_element)

    legend_entries = build_legend_entries(records, palette)
    y_offset = 40
    for name, color in legend_entries:
        add_rect(annotation_group, LEGEND_SWATCH_X, y_offset, LEGEND_SWATCH_WIDTH, LEGEND_SWATCH_HEIGHT, color)
        text_element = ET.Element(
            "text",
            {
                "x": str(LEGEND_TEXT_X),
                "y": str(y_offset + 14),
                "fill": "black",
                "font-size": str(LEGEND_FONT_SIZE),
            },
        )
        text_element.text = name
        annotation_group.append(text_element)
        y_offset += LEGEND_ROW_STEP
    if centromere_overlay:
        add_centromere_legend_swatch(annotation_group, y_offset, centromere_overlay_color, centromere_overlay_opacity)
        text_element = ET.Element(
            "text",
            {
                "x": str(LEGEND_TEXT_X),
                "y": str(y_offset + 14),
                "fill": "black",
                "font-size": str(LEGEND_FONT_SIZE),
            },
        )
        text_element.text = CENTROMERE_LEGEND_LABEL
        annotation_group.append(text_element)
        y_offset += LEGEND_ROW_STEP
    if show_acrocentric_p_arm_legend:
        add_rect(
            annotation_group,
            LEGEND_SWATCH_X,
            y_offset,
            LEGEND_SWATCH_WIDTH,
            LEGEND_SWATCH_HEIGHT,
            ACROCENTRIC_P_ARM_COLOR,
        )
        text_element = ET.Element(
            "text",
            {
                "x": str(LEGEND_TEXT_X),
                "y": str(y_offset + 14),
                "fill": "black",
                "font-size": str(LEGEND_FONT_SIZE),
            },
        )
        text_element.text = ACROCENTRIC_P_ARM_LEGEND_LABEL
        annotation_group.append(text_element)

    original_svg_root.insert(svg_child_insert_index(original_svg_root), region_group)
    original_svg_root.append(annotation_group)

    os.makedirs(os.path.dirname(output_svg_file) or ".", exist_ok=True)
    original_svg_tree.write(output_svg_file, encoding="utf-8", xml_declaration=True)

    if not export_pdf:
        LOGGER.info("SVG written: %s", output_svg_file)
        return

    if not shutil.which("rsvg-convert"):
        raise RuntimeError("rsvg-convert not found. Install librsvg before exporting PDF, or use --svg-only.")

    os.makedirs(os.path.dirname(output_pdf_file) or ".", exist_ok=True)
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", output_pdf_file, output_svg_file], check=True)
    LOGGER.info("File converted to PDF: %s", output_pdf_file)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Insert colored LAP BED regions into the hg38 SVG template.")
    parser.add_argument("-B", type=parse_svg_filename, default="hg38", help='Input build38 "hg38" SVG file without extension')
    parser.add_argument("-I", type=str, required=True, help="Input LAP BED file")
    parser.add_argument("-O", type=str, required=True, help="Output path: .pdf, .svg, or basename/prefix")
    parser.add_argument("--color-config", type=str, help="JSON file mapping ancestry labels to colors.")
    parser.add_argument(
        "--overlap-px",
        type=float,
        default=0.25,
        help="Vertical SVG-pixel overlap between ancestry rectangles to avoid antialiasing gaps.",
    )
    parser.add_argument(
        "--missing-color",
        type=str,
        help="Optional fill color for chromosomes present in the BED before ancestry segments are drawn.",
    )
    parser.add_argument(
        "--snp-line-width",
        type=float,
        default=3.0,
        help="Stroke width for FLARE SNP ancestry lines.",
    )
    parser.add_argument(
        "--snp-line-dasharray",
        default="2 1.5",
        help="SVG stroke-dasharray for FLARE SNP lines; use an empty string for solid lines.",
    )
    parser.add_argument(
        "--snp-line-overhang-px",
        type=float,
        default=0.0,
        help="Horizontal overhang in SVG pixels for FLARE SNP lines beyond each haplotype bar.",
    )
    parser.add_argument(
        "--snp-line-inset-px",
        type=float,
        default=1.0,
        help="Horizontal inset in SVG pixels for FLARE SNP lines inside each haplotype bar.",
    )
    parser.add_argument(
        "--snp-line-halo-color",
        default="#ffffff",
        help="Contrast stroke drawn under FLARE SNP lines; use an empty string to disable.",
    )
    parser.add_argument(
        "--snp-line-halo-width",
        type=float,
        default=4.0,
        help="Stroke width for the contrast line drawn under FLARE SNP lines.",
    )
    parser.add_argument("--no-snp-labels", action="store_true", help="Do not draw short SNP labels above FLARE SNP lines.")
    parser.add_argument(
        "--snp-label-font-size",
        type=float,
        default=5.5,
        help="Font size for haplotype-specific SNP labels drawn above FLARE SNP lines.",
    )
    parser.add_argument(
        "--snp-label-offset-px",
        type=float,
        default=3.0,
        help="Vertical offset in SVG pixels between each SNP line and its short label.",
    )
    parser.add_argument(
        "--snp-label-halo-color",
        default="#ffffff",
        help="Contrast outline color for short SNP labels; use an empty string to disable.",
    )
    parser.add_argument(
        "--pdf-backend",
        choices=("rsvg", "native"),
        default="rsvg",
        help="PDF renderer: rsvg uses the legacy SVG template; native draws a cleaner vector PDF directly from the BED.",
    )
    parser.add_argument(
        "--no-centromere-overlay",
        action="store_true",
        help="Do not draw the semi-transparent segmented centromere overlay.",
    )
    parser.add_argument(
        "--centromere-overlay-color",
        default="#707070",
        help="Fill color for the semi-transparent segmented centromere overlay.",
    )
    parser.add_argument(
        "--centromere-overlay-opacity",
        type=float,
        default=0.42,
        help="Opacity for the segmented centromere overlay, from 0 to 1.",
    )
    parser.add_argument(
        "--auto-zoom-overlaps",
        action="store_true",
        help="Write secondary native PDF zoom panels for clusters of FLARE SNP labels that overlap in the main karyotype.",
    )
    parser.add_argument(
        "--zoom-output-dir",
        type=str,
        help="Directory for auto-generated zoom PDFs. Default: <output-basename>_zoom.",
    )
    parser.add_argument(
        "--zoom-overlap-px",
        type=float,
        default=ZOOM_DEFAULT_OVERLAP_PX,
        help="Projected SVG-pixel distance below which neighboring SNPs are grouped into one zoom panel.",
    )
    parser.add_argument(
        "--zoom-padding-bp",
        type=int,
        default=ZOOM_DEFAULT_PADDING_BP,
        help="Base-pair padding added to each detected SNP cluster in zoom panels.",
    )
    parser.add_argument(
        "--zoom-min-window-bp",
        type=int,
        default=ZOOM_DEFAULT_MIN_WINDOW_BP,
        help="Minimum genomic span shown in each zoom panel.",
    )
    parser.add_argument(
        "--zoom-min-snps",
        type=int,
        default=2,
        help="Minimum number of close SNPs required to generate a zoom panel.",
    )
    parser.add_argument(
        "--zoom-max-panels",
        type=int,
        default=ZOOM_DEFAULT_MAX_PANELS,
        help="Maximum number of automatic zoom panels to write per sample.",
    )
    parser.add_argument("--svg-only", action="store_true", help="Write SVG only and skip PDF conversion.")
    parser.add_argument("--strict", action="store_true", help="Fail instead of warning on unsupported chromosomes.")
    parser.add_argument("--quiet", action="store_true", help="Only warnings and errors.")
    return parser.parse_args(argv)


def setup_logging(args: argparse.Namespace) -> None:
    level = logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    setup_logging(args)

    input_name = os.path.basename(args.I)
    individual_name = input_name[:-4] if input_name.lower().endswith(".bed") else os.path.splitext(input_name)[0]
    out_svg, out_pdf = parse_output_paths(args.O)

    def maybe_write_zoom_panels() -> None:
        if not args.auto_zoom_overlaps:
            return
        out_pdf_path = Path(out_pdf)
        zoom_output_dir = args.zoom_output_dir or f"{out_pdf_path.with_suffix('')}_zoom"
        write_native_zoom_overlap_plots(
            args.I,
            zoom_output_dir,
            out_pdf_path.stem,
            individual_name,
            args.color_config,
            args.missing_color,
            args.strict,
            args.snp_line_width,
            args.snp_line_dasharray or None,
            args.snp_line_overhang_px,
            args.snp_line_inset_px,
            args.snp_line_halo_color or None,
            args.snp_line_halo_width,
            not args.no_snp_labels,
            args.snp_label_font_size,
            args.snp_label_halo_color or None,
            not args.no_centromere_overlay,
            args.centromere_overlay_color,
            args.centromere_overlay_opacity,
            args.zoom_overlap_px,
            args.zoom_padding_bp,
            args.zoom_min_window_bp,
            args.zoom_min_snps,
            args.zoom_max_panels,
        )

    if args.pdf_backend == "native" and not args.svg_only:
        write_native_pdf(
            args.I,
            out_pdf,
            individual_name,
            args.color_config,
            args.overlap_px,
            args.missing_color,
            args.strict,
            args.snp_line_width,
            args.snp_line_dasharray or None,
            args.snp_line_overhang_px,
            args.snp_line_inset_px,
            args.snp_line_halo_color or None,
            args.snp_line_halo_width,
            not args.no_snp_labels,
            args.snp_label_font_size,
            args.snp_label_offset_px,
            args.snp_label_halo_color or None,
            not args.no_centromere_overlay,
            args.centromere_overlay_color,
            args.centromere_overlay_opacity,
        )
        maybe_write_zoom_panels()
        return 0

    insert_colored_regions(
        args.B,
        args.I,
        out_svg,
        out_pdf,
        individual_name,
        args.color_config,
        args.overlap_px,
        args.missing_color,
        args.strict,
        args.snp_line_width,
        args.snp_line_dasharray or None,
        args.snp_line_overhang_px,
        args.snp_line_inset_px,
        args.snp_line_halo_color or None,
        args.snp_line_halo_width,
        not args.no_snp_labels,
        args.snp_label_font_size,
        args.snp_label_offset_px,
        args.snp_label_halo_color or None,
        not args.no_centromere_overlay,
        args.centromere_overlay_color,
        args.centromere_overlay_opacity,
        not args.svg_only,
    )
    maybe_write_zoom_panels()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
