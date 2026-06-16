# LAP SNP Zoom Prototype

This experimental folder adds optional native PDF zoom panels for FLARE SNP markers that are too close to read on the full chromosome karyotype.

Example:

```bash
python3 Plot_LAP_v4.py \
  -I Out-Test/mock_zoom_overlap.bed \
  -B hg38.svg \
  -O Out-Test/mock_zoom_overlap_native.pdf \
  --pdf-backend native \
  --auto-zoom-overlaps
```

The main plot is written normally. When close SNP clusters are detected, extra PDFs are written to:

```text
<output-basename>_zoom/
```

Useful options:

```text
--auto-zoom-overlaps      enable automatic secondary zoom panels
--zoom-output-dir DIR     write zoom PDFs to a custom directory
--zoom-overlap-px 18      projected main-plot distance used to group close SNPs
--zoom-padding-bp 250000  genomic padding around each cluster
--zoom-min-window-bp 500000
--zoom-min-snps 2
--zoom-max-panels 12
```

The zoom panels are native-vector PDFs and do not depend on editing the SVG template.
