# AncestryGrapher Toolkit

Authors: Alessandro Lisi and Michael C. Campbell  
Human Evolutionary Genomics Lab, Department of Biological Sciences, Human and Evolutionary Biology section, University of Southern California

AncestryGrapher is a command-line toolkit for converting RFMix v2 local and global ancestry output into publication-ready ancestry visualizations. The toolkit contains two companion workflows:

- **GAP, Global Ancestry Painting:** summarizes and plots genome-wide ancestry proportions from RFMix `.rfmix.Q` files.
- **LAP, Local Ancestry Painting:** paints local ancestry tracts along human autosomes from RFMix `.msp.tsv` files.

The current LAP workflow also supports **FLARE ancestry VCFs** (`.anc.vcf` / `.anc.vcf.gz`) for SNP-aware visualization and variant-level reporting. When a user supplies a list of variants, LAP can draw haplotype-resolved SNP markers on the karyotype and write a LAP-VAR report describing the ancestry of the haplotype carrying each allele.

The repository keeps the original RFMix-oriented structure:

```text
RFMIX2-Pipeline-to-plot/
├── GAP/
│   └── Scripts/
│       ├── Gap_v2.py
│       ├── GAP_Plot.py
│       └── GAP.py                 # compatibility wrapper for GAP_Plot.py
├── LAP/
│   ├── hg38.svg
│   ├── hg37.svg
│   ├── README_LAP_FLARE.md
│   ├── README_ZOOM_PANELS.md
│   ├── Examples/
│   │   └── mock_zoom_overlap.bed
│   └── Scripts/
│       ├── LAP_v4.py
│       ├── Plot_LAP_v4.py
│       └── LAP.py                 # compatibility wrapper for Plot_LAP_v4.py
├── Genetic_Map/
├── Sample_map_File/
├── Example_Dataset/
└── README.md
```

The historical helper scripts are still present for compatibility, but the recommended current workflows are `Gap_v2.py` for GAP and `LAP_v4.py` plus `Plot_LAP_v4.py` for LAP.

---

## 1. Installation

AncestryGrapher is written in Python and runs on macOS and Linux. Windows users can run it through Anaconda, WSL, or another Unix-like Python environment.

Clone the repository:

```bash
git clone https://github.com/alisi1989/RFMIX2-Pipeline-to-plot.git
cd RFMIX2-Pipeline-to-plot
```

Install Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

The core Python dependencies are:

```text
pandas
numpy
matplotlib
```

For the historical SVG-to-PDF backend, install `librsvg`:

```bash
brew install librsvg
```

On Debian/Ubuntu:

```bash
sudo apt-get install -y librsvg2-bin
```

If `rsvg-convert` is not available, LAP can still write SVG files with `--svg-only`, or can write direct vector PDFs with `--pdf-backend native`.

---

## 2. Input Data Before AncestryGrapher

AncestryGrapher expects RFMix v2 output. The usual upstream workflow is:

1. Phase target and reference VCFs.
2. Run RFMix v2 chromosome by chromosome.
3. Use the RFMix `.rfmix.Q` files for GAP.
4. Use the RFMix `.msp.tsv` files for LAP.
5. Optionally, use FLARE `.anc.vcf.gz` files plus a SNP list for LAP-VAR and SNP-aware plotting.

### Target VCFs

The target VCF should be phased. RFMix is commonly run chromosome by chromosome. The exact input layout can vary, but the final RFMix output names should make the chromosome recoverable from the filename.

Recommended naming:

```text
target_prefix_chr1.rfmix.Q
target_prefix_chr1.msp.tsv
target_prefix_chr2.rfmix.Q
target_prefix_chr2.msp.tsv
...
target_prefix_chr22.rfmix.Q
target_prefix_chr22.msp.tsv
```

For single-sample legacy runs, filenames may include the sample name:

```text
Mozabite1_chr4.rfmix.Q
Mozabite1_chr4.msp.tsv
```

For multi-sample runs, one file per chromosome can contain many target individuals:

```text
target_covid_kinPass_chr4.rfmix.Q
target_covid_kinPass_chr4.msp.tsv
```

The current GAP and LAP scripts handle the multi-sample RFMix layout.

### Reference Sample Map

RFMix requires a sample map for the reference individuals. The file is tab-delimited and contains no target individuals:

```text
ind1    Africa
ind2    Africa
ind3    Africa
ind4    Europe
ind5    Europe
ind6    Middle_East
ind7    Middle_East
ind8    Middle_East
```

The ancestry labels in this file become the ancestry names used by RFMix and then by AncestryGrapher. Keep labels stable and avoid accidental spelling differences such as `MiddleEast` versus `Middle_East`.

### Genetic Map

The `Genetic_Map/` directory contains hg38 genetic maps for chromosomes 1 through 22. RFMix expects genetic map files with columns similar to:

```text
chr    pos     cM
chr1   55550   0
chr1   82571   0.080572
chr1   88169   0.092229
chr1   285245  0.439456
```

### Example RFMix v2 Command

The exact command depends on how the user's data are organized. A typical loop is:

```bash
for chr in {1..22}; do
  rfmix \
    -f Example_Dataset/Target/Mozabite_1.vcf.gz \
    -r Example_Dataset/Reference/Reference_Phased_chr${chr}.vcf.gz \
    -m Sample_map_File/Sample_Reference.txt \
    -g Genetic_Map/chr${chr}.b38.txt \
    -o Example_Dataset/RFMIX2_Output/Mozabite1_chr${chr} \
    --chromosome=${chr}
done
```

Important RFMix outputs:

- `*.rfmix.Q`: global ancestry proportions by chromosome, used by GAP.
- `*.msp.tsv`: local ancestry tracts, used by LAP.
- `*.fb.tsv` and `*.sis.tsv`: additional RFMix outputs, not required by the current AncestryGrapher workflows.

---

## 3. GAP: Global Ancestry Painting

GAP converts RFMix `.rfmix.Q` files into a `.gap` table and then plots stacked global ancestry bars.

Current recommended scripts:

```text
GAP/Scripts/Gap_v2.py
GAP/Scripts/GAP_Plot.py
```

`GAP/Scripts/GAP.py` is a compatibility wrapper for the plotter.

### GAP Step 1: Build a `.gap` Table

From the repository root:

```bash
cd GAP
python3 Scripts/Gap_v2.py \
  --prefix ../Example_Dataset/RFMIX2_Output/Mozabite \
  --chr 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 \
  --output-dir Output_GAP \
  --final-output Output_GAP/Mozabite.gap
```

You can also use shell expansion:

```bash
python3 Scripts/Gap_v2.py \
  --prefix ../Example_Dataset/RFMIX2_Output/Mozabite \
  --chr {1..22} \
  --output-dir Output_GAP \
  --final-output Output_GAP/Mozabite.gap
```

If `--chr` is omitted, GAP uses autosomes 1 through 22.

### GAP File Discovery

`--prefix` is the shared prefix before the chromosome token. For example, if files are:

```text
../Example_Dataset/RFMIX2_Output/Mozabite1_chr1.rfmix.Q
../Example_Dataset/RFMIX2_Output/Mozabite1_chr2.rfmix.Q
```

then a prefix like this is appropriate:

```bash
--prefix ../Example_Dataset/RFMIX2_Output/Mozabite1
```

If files are multi-sample per chromosome:

```text
target_covid_kinPass_chr1.rfmix.Q
target_covid_kinPass_chr2.rfmix.Q
```

then use:

```bash
--prefix target_covid_kinPass
```

### GAP Weighting

By default, GAP combines chromosome-level RFMix Q values using hg38 chromosome-length weighting:

```bash
--weighting chrom-length
```

This is usually preferable for genome-wide global ancestry estimates because longer chromosomes contribute more genomic material. If users want every chromosome to contribute equally, use:

```bash
--weighting equal
```

### Missing Chromosomes

By default, if a sample is missing one requested chromosome, GAP computes that sample's ancestry from the chromosomes available for that sample and logs a warning. To drop incomplete samples:

```bash
--require-all-chroms
```

### Sorting

To sort individuals by a specific ancestry column:

```bash
python3 Scripts/Gap_v2.py \
  --prefix ../Example_Dataset/RFMIX2_Output/Mozabite \
  --output-dir Output_GAP \
  --final-output Output_GAP/Mozabite.sorted.gap \
  --sort-ancestry Middle_East
```

The ancestry name must match the RFMix Q column exactly.

### Custom GAP Colors

Colors can be provided as JSON:

```json
{
  "Africa": "#a38905",
  "Europe": "#a30d05",
  "Middle_East": "#0e6b05"
}
```

Use it with:

```bash
python3 Scripts/Gap_v2.py \
  --prefix ../Example_Dataset/RFMIX2_Output/Mozabite \
  --output-dir Output_GAP \
  --final-output Output_GAP/Mozabite.gap \
  --color-config colors.json
```

### GAP Step 2: Plot Global Ancestry

```bash
python3 Scripts/GAP_Plot.py \
  --input Output_GAP/Mozabite.gap \
  --output Output_GAP/Mozabite_GAP.pdf \
  --title "Mozabite Global Ancestry"
```

SVG output is also supported:

```bash
python3 Scripts/GAP_Plot.py \
  --input Output_GAP/Mozabite.gap \
  --output Output_GAP/Mozabite_GAP.svg
```

If a `.gap` file contains values that do not sum exactly to 1 per sample, use:

```bash
--normalize
```

---

## 4. LAP: Local Ancestry Painting

LAP converts RFMix `.msp.tsv` local ancestry files into one final BED-like file per sample and then draws each sample as a karyotype-like ancestry painting.

Current recommended scripts:

```text
LAP/Scripts/LAP_v4.py
LAP/Scripts/Plot_LAP_v4.py
```

`LAP/Scripts/LAP.py` is a compatibility wrapper for `Plot_LAP_v4.py`.

### LAP Step 1: Build Final LAP BED Files

From the repository root:

```bash
cd LAP
python3 Scripts/LAP_v4.py \
  --prefix ../Example_Dataset/RFMIX2_Output/Mozabite \
  --chr {1..22} \
  --output-dir Output_LAP
```

If `--chr` is omitted, LAP uses `chr1` through `chr22`.

The output directory receives one final BED-like file per sample:

```text
Output_LAP/SAMPLE.bed
```

Each final BED contains ancestry rectangles for haplotype 1 and haplotype 2. The plotter reads these records and paints each chromosome pair.

### Keeping Intermediate Files

By default, LAP v4 writes only the final per-sample BEDs. To keep intermediate combined MSP and haplotype BED files:

```bash
--keep-temp
```

### Parallel Processing

For larger cohorts:

```bash
python3 Scripts/LAP_v4.py \
  --prefix Input/target_chr \
  --output-dir Output_LAP \
  --threads 4
```

The processing is I/O-bound, so moderate thread counts are usually enough.

### Merging Adjacent Windows

By default, LAP merges consecutive same-ancestry RFMix windows to make the final BED cleaner and plotting faster. To keep every RFMix window:

```bash
--no-merge-adjacent
```

To merge same-ancestry windows separated by a small gap:

```bash
--max-merge-gap-bp 1000
```

### Custom LAP Colors

Colors can be assigned by ancestry name:

```bash
python3 Scripts/LAP_v4.py \
  --prefix Input/target_chr \
  --output-dir Output_LAP \
  --color-config colors.json
```

The JSON file should map ancestry labels to hex colors:

```json
{
  "Europeans": "#d18311",
  "Middle_East": "#22ba9d",
  "North_Africa": "#839dfc",
  "SubSaharan_Africa": "#9a5dc1"
}
```

Legacy positional color flags are also supported:

```bash
--ancestry0 "#a32e2e" --ancestry1 "#0a0ae0"
```

### Highlighting Regions

LAP can add dashed black lines to mark genes, candidate loci, GWAS regions, or any user-defined genomic interval.

For one region:

```bash
python3 Scripts/LAP_v4.py \
  --prefix Input/target_chr \
  --output-dir Output_LAP \
  -c chr4 \
  --from-bp 46500000 \
  --to-bp 46600000
```

For many regions, use `--highlight-regions` or its alias `--regions-to-highlight`:

```bash
python3 Scripts/LAP_v4.py \
  --prefix Input/target_chr \
  --output-dir Output_LAP \
  --highlight-regions regions_to_highlight.tsv
```

Format:

```text
chrom  start      end        label
chr4   46500000   46600000   locus_A
chr6   25000000   25500000   locus_B
17     43000000   43100000   locus_C
```

The `label` column is optional. Highlighted regions do not alter FLARE SNP annotation or LAP-VAR reporting.

---

## 5. LAP Step 2: Plot Local Ancestry

Historical RSVG backend:

```bash
python3 Scripts/Plot_LAP_v4.py \
  -I Output_LAP/SAMPLE.bed \
  -B hg38.svg \
  -O Output_LAP/SAMPLE.pdf \
  --pdf-backend rsvg
```

Compatibility wrapper:

```bash
python3 Scripts/LAP.py \
  -I Output_LAP/SAMPLE.bed \
  -B hg38.svg \
  -O Output_LAP/SAMPLE.pdf \
  --pdf-backend rsvg
```

Direct native PDF backend:

```bash
python3 Scripts/Plot_LAP_v4.py \
  -I Output_LAP/SAMPLE.bed \
  -O Output_LAP/SAMPLE.native.pdf \
  --pdf-backend native
```

The RSVG backend uses the editable SVG template (`hg38.svg`) and converts it to PDF with `rsvg-convert`. This is useful when users want the exact template shapes.

The native backend draws a vector PDF directly from Python. It avoids SVG conversion and produces a simpler layer structure for Illustrator, Inkscape, or Affinity editing. It uses the same hg38 chromosome lengths and coordinate mapping as the template-based plot.

### SVG-Only Output

If `rsvg-convert` is not installed:

```bash
python3 Scripts/Plot_LAP_v4.py \
  -I Output_LAP/SAMPLE.bed \
  -B hg38.svg \
  -O Output_LAP/SAMPLE.svg \
  --svg-only
```

### Centromere Overlay

By default, `Plot_LAP_v4.py` draws a semi-transparent, segmented grey overlay on centromere coordinates. This is intended to make centromeric regions visually explicit because local ancestry inference can be less reliable in these regions.

Disable it:

```bash
--no-centromere-overlay
```

Adjust color and opacity:

```bash
--centromere-overlay-color "#707070" \
--centromere-overlay-opacity 0.42
```

The legend includes a `Centromere` entry. The native backend also draws and labels grey acrocentric p-arm/satellite regions for chromosomes 13, 14, 15, 21, and 22. In the RSVG backend, the legend includes `Acrocentric p-arms/satellites` when the SVG template contains those grey regions.

---

## 6. FLARE-Aware SNP Annotation

LAP v4 can integrate FLARE ancestry calls at specific variants. This allows the plot to display SNP markers colored by the FLARE haplotype ancestry and allows LAP to write per-sample variant ancestry reports.

Required inputs:

- RFMix `.msp.tsv` files for local ancestry tracts.
- FLARE `.anc.vcf` or `.anc.vcf.gz` with phased `GT` and ancestry fields such as `AN1` and `AN2`.
- A SNP list with variants of interest.

Example:

```bash
python3 Scripts/LAP_v4.py \
  --prefix Input/target_covid_kinPass_chr \
  --output-dir Output_LAP \
  --flare-vcf Flare_output/flare_chr4_covidanc.vcf.gz \
  --snp-list Flare_output/snps.txt
```

Then plot:

```bash
python3 Scripts/Plot_LAP_v4.py \
  -I Output_LAP/SAMPLE.bed \
  -B hg38.svg \
  -O Output_LAP/SAMPLE.pdf \
  --pdf-backend rsvg
```

### FLARE Code Mapping

RFMix and FLARE can encode ancestries in different numeric orders. LAP does not assume that RFMix ancestry code `1` means the same thing as FLARE ancestry code `1`.

The mapping is:

```text
FLARE AN1/AN2 code -> FLARE ancestry name -> LAP/RFMix color for the same ancestry name
```

Therefore, as long as the ancestry names match, SNP markers are colored correctly even if numeric codes differ.

### SNP List Formats

One SNP per line. Examples:

```text
chr4 46567563
chr4:46567563
4:46567563:G:A
```

Headered tabular format is also supported:

```text
CHROM POS ID
chr4 46567563 4:46567563:G:A
```

SNP labels in plots are assigned by the order of the SNP list:

```text
SNP1, SNP2, SNP3, ...
```

On the plot, labels are haplotype-specific:

```text
SNP1/AN1
SNP1/AN2
```

### Sample Matching Between LAP and FLARE

By default, LAP tries to match the sample name from RFMix/MSP to the sample name in the FLARE VCF.

For a single FLARE sample:

```bash
--flare-sample COV.COV190_111
```

For a mapping file:

```bash
--flare-sample-map sample_map.tsv
```

Format:

```text
LAP_sample        VCF_sample
sampleA           COV.COV190_111
sampleB           COV.COV524_286
```

---

## 7. LAP-VAR: Variant-Aware Local Ancestry Reporting

When `--flare-vcf` and `--snp-list` are supplied, LAP writes a tab-delimited report for every sample:

```text
Output_LAP/SAMPLE.variant_ancestry.txt
```

The report includes:

```text
sample
vcf_sample
snp_label
chrom
pos
variant_id
REF
ALT
GT
phased
hap1_allele
hap2_allele
hap1_FLARE_ancestry
hap2_FLARE_ancestry
hap1_RFMix_segment_ancestry
hap2_RFMix_segment_ancestry
ALT_allele_haplotype
ALT_allele_FLARE_ancestry
ALT_allele_RFMix_ancestry
FLARE_RFMix_concordance
ANP1
ANP2
confidence
interpretation
```

The biologically important point is haplotype-resolved interpretation. For example:

```text
GT = 0|1
AN1 = Europeans
AN2 = Americas
```

The ALT allele is carried on haplotype 2, so the report states that the ALT allele is embedded in an Americas local ancestry background within that analyzed individual.

This does **not** imply that the mutation itself originated in that ancestry group. It reports the local ancestry of the haplotype carrying the allele in the individual being analyzed.

This is useful for candidate variants, GWAS loci, rare variant follow-up, medical genetics, pharmacogenomics, and population-genetic interpretation of admixed haplotypes.

---

## 8. SNP Marker Plot Styling

The plotter draws FLARE SNP markers as short dashed colored lines on each haplotype. Defaults are tuned to be readable on the full karyotype:

```text
--snp-line-inset-px 1.0
--snp-line-dasharray "2 1.5"
--snp-line-halo-width 4.0
--snp-label-font-size 5.5
```

Useful options:

```bash
--snp-line-width 3.0
--snp-line-dasharray "2 1.5"
--snp-line-overhang-px 0.0
--snp-line-inset-px 1.0
--snp-line-halo-color "#ffffff"
--snp-line-halo-width 4.0
--no-snp-labels
--snp-label-font-size 5.5
--snp-label-offset-px 3.0
--snp-label-halo-color "#ffffff"
```

Use an empty dasharray for solid SNP markers:

```bash
--snp-line-dasharray ""
```

---

## 9. Automatic Zoom Panels for Overlapping SNPs

When variants are very close, labels and dashed SNP markers can overlap in the full chromosome plot. The current plotter can automatically create secondary zoom panels:

```bash
python3 Scripts/Plot_LAP_v4.py \
  -I Output_LAP/SAMPLE.bed \
  -B hg38.svg \
  -O Output_LAP/SAMPLE.pdf \
  --pdf-backend rsvg \
  --auto-zoom-overlaps
```

The main plot is written normally. If close SNP clusters are detected, extra native-vector PDFs are written to:

```text
Output_LAP/SAMPLE_zoom/
```

The zoom panels show a local haplotype-resolved view with:

- hap1 and hap2 as enlarged local tracks
- ancestry painting across the zoomed interval
- dashed SNP lines colored by FLARE `AN1`/`AN2`
- readable `SNP#/AN#` labels with leader lines
- genomic coordinates for the zoom window
- a small chromosome context indicator

Useful options:

```bash
--zoom-output-dir DIR
--zoom-overlap-px 18
--zoom-padding-bp 250000
--zoom-min-window-bp 500000
--zoom-min-snps 2
--zoom-max-panels 12
```

The file `LAP/Examples/mock_zoom_overlap.bed` is a small synthetic example for testing:

```bash
cd LAP
python3 Scripts/Plot_LAP_v4.py \
  -I Examples/mock_zoom_overlap.bed \
  -B hg38.svg \
  -O Output_LAP/mock_zoom_overlap.pdf \
  --pdf-backend native \
  --auto-zoom-overlaps
```

---

## 10. Final BED Format Used by LAP

The final LAP BED is a tab-delimited file. Core ancestry records use:

```text
chrom  start  end  geom_rect  color  haplotype
```

Highlighted region records use:

```text
chrom  start  end  geom_line  color  haplotype  label
```

FLARE SNP records use:

```text
chrom  start  end  geom_snp  color  haplotype  variant_id  vcf_sample  ancestry_label  ancestry_code  snp_label
```

Example:

```text
4  46567563  46567563  geom_snp  #22ba9d  1  4:46567563:G:A  COV.COV190_111  Middle_East  1  SNP1
4  46567563  46567563  geom_snp  #9a5dc1  2  4:46567563:G:A  COV.COV190_111  SubSaharan_Africa  3  SNP1
```

---

## 11. Troubleshooting

### `rsvg-convert not found`

Install `librsvg` or use native PDF output:

```bash
python3 Scripts/Plot_LAP_v4.py \
  -I Output_LAP/SAMPLE.bed \
  -O Output_LAP/SAMPLE.native.pdf \
  --pdf-backend native
```

### SNP colors do not match expectations

Check that ancestry names match between RFMix/MSP and FLARE. Numeric codes can differ; names are what matter. For example, `Europeans` in FLARE must correspond to `Europeans` in the LAP/RFMix palette.

### FLARE sample names do not match LAP sample names

Use:

```bash
--flare-sample SAMPLE_NAME
```

or:

```bash
--flare-sample-map sample_map.tsv
```

### Some chromosomes are blank

This usually means the BED has no ancestry records for those chromosomes, or the requested chromosome list did not match the file naming convention. Check `--prefix`, `--chr`, and chromosome tokens such as `1` versus `chr1`.

### GAP ancestry sums are slightly off

Use `GAP_Plot.py --normalize` for plotting if needed. Also inspect missing chromosome warnings from `Gap_v2.py`.

---

## 12. Citation and Contact

For questions about the pipeline, please contact:

- Alessandro Lisi: alisi@usc.edu
- Michael C. Campbell: mc44680@usc.edu

Please cite the AncestryGrapher toolkit and the underlying local ancestry inference software used in your analysis.
