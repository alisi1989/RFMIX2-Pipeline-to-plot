# LAP + FLARE SNP annotations

Questa versione aggiunge a LAP il supporto agli output FLARE `.anc.vcf` / `.anc.vcf.gz` con campi `FORMAT/AN1` e `FORMAT/AN2`.

## Logica implementata

1. LAP legge il mapping ancestry dell'MSP RFMix dal commento:

   `#Subpopulation order/codes: Americas=0 ... Europeans=3 ...`

2. LAP legge il mapping ancestry del VCF FLARE dal commento:

   `##ANCESTRY=<Europeans=0,East_Asia=1,Americas=2,...>`

3. I codici FLARE `AN1` e `AN2` vengono tradotti prima nei nomi ancestry del VCF e poi colorati usando la palette ancestry dell'MSP/LAP. Questo evita errori quando i codici numerici FLARE e RFMix non hanno lo stesso ordine.

4. Per ogni SNP annotato vengono scritte due righe nel BED finale:

   - `geom_snp` haplotype `1`: colore da `AN1`
   - `geom_snp` haplotype `2`: colore da `AN2`

   Nel plot, le due linee sono disegnate sui due cromosomi omologhi: la metà/hap1 prende il colore di `AN1`, la metà/hap2 prende il colore di `AN2`.

5. Per ogni campione LAP scrive anche un report LAP-VAR:

   ```text
   SAMPLE.variant_ancestry.txt
   ```

   Il report collega genotipo phased, allele su hap1/hap2, ancestry FLARE e ancestry del segmento RFMix che contiene la variante.

## Uso base

```bash
python LAP_v4.py \
  --prefix /path/to/rfmix/prefix_ \
  --chr chr4 \
  --output-dir out_bed \
  --flare-vcf flare_chr4.anc.vcf.gz \
  --snp-list snps.txt
```

Poi:

```bash
python Plot_LAP_v4.py \
  -B hg38.svg \
  -I out_bed/SAMPLE.bed \
  -O plots/SAMPLE.pdf
```

Se `rsvg-convert` non è installato, puoi produrre solo SVG:

```bash
python Plot_LAP_v4.py \
  -B hg38.svg \
  -I out_bed/SAMPLE.bed \
  -O plots/SAMPLE.svg \
  --svg-only
```

## Formati accettati per `--snp-list`

Una riga per SNP. Sono accettati, per esempio:

```text
chr4 46567563
chr4:46567563
4:46567563:G:A
```

Sono supportati anche file tabulari con header tipo:

```text
CHROM POS ID
chr4 46567563 4:46567563:G:A
```

## Sample matching

Per default LAP prova a usare lo stesso nome sample presente nelle colonne MSP (`sample.0`, `sample.1`) e nel VCF FLARE.

Se il nome non coincide:

```bash
python LAP_v4.py ... --flare-sample-map sample_map.tsv
```

con `sample_map.tsv` a due colonne:

```text
LAP_sample  VCF_sample
sampleA     COV.COV190_111
```

Per un singolo individuo puoi forzare il sample VCF:

```bash
python LAP_v4.py ... --flare-sample COV.COV190_111
```

## Nuove opzioni principali

- `--flare-vcf`: VCF/VCF.GZ FLARE con `GT:AN1:AN2` o comunque `AN1`/`AN2` nel FORMAT.
- `--snp-list`: TXT/TSV con SNP da annotare.
- `--flare-sample`: usa un unico sample VCF per tutti gli output LAP.
- `--flare-sample-map`: mappa `LAP_sample -> VCF_sample`.
- `--highlight-regions` / `--regions-to-highlight`: TXT/TSV con regioni da evidenziare con linee nere tratteggiate (`chrom start end [label]`).
- `-c/--chromosome`, `--from-bp`, `--to-bp`: evidenzia una singola regione manuale, mantenuto per compatibilita.
- `Plot_LAP_v4.py --snp-line-width`: spessore delle linee SNP.
- `Plot_LAP_v4.py --snp-line-dasharray`: tratteggio SVG; stringa vuota per linea continua. Il default usa segmenti brevi e fitti.
- `Plot_LAP_v4.py --snp-line-overhang-px`: estensione orizzontale oltre ciascuna barra haplotipica.
- `Plot_LAP_v4.py --snp-line-inset-px`: accorcia leggermente le linee SNP dentro ciascuna barra, evitando che coprano il bordo nero.
- `Plot_LAP_v4.py --snp-line-halo-color`: colore di contrasto sotto le linee SNP; stringa vuota per disattivarlo.
- `Plot_LAP_v4.py --snp-line-halo-width`: spessore della linea di contrasto sotto le linee SNP.
- `Plot_LAP_v4.py --no-snp-labels`: nasconde le etichette haplotype-specific `SNP1/AN1`, `SNP1/AN2`, ... sopra le linee SNP.
- `Plot_LAP_v4.py --snp-label-font-size`: dimensione delle etichette SNP haplotype-specific sopra le linee SNP.
- `Plot_LAP_v4.py --svg-only`: salta la conversione PDF.
- `Plot_LAP_v4.py --pdf-backend native`: scrive un PDF vettoriale diretto da Python, senza passare da `hg38.svg` / `rsvg-convert`.
- `Plot_LAP_v4.py --pdf-backend rsvg`: comportamento storico, usa `hg38.svg` e `rsvg-convert`.
- `Plot_LAP_v4.py --no-centromere-overlay`: disattiva il layer grigio segmentato sulle regioni centromeriche.
- `Plot_LAP_v4.py --centromere-overlay-color`: colore del layer centromerico.
- `Plot_LAP_v4.py --centromere-overlay-opacity`: trasparenza del layer centromerico, da `0` a `1`.

## Report LAP-VAR

Quando `--flare-vcf` e `--snp-list` sono usati, LAP scrive un file tab-delimited per campione:

```text
out_dir/SAMPLE.variant_ancestry.txt
```

Le colonne principali includono:

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

Interpretazione importante: `ALT_allele_FLARE_ancestry` indica l'ancestry locale dell'haplotype che porta l'allele ALT nell'individuo analizzato. Non implica necessariamente l'origine evolutiva della mutazione.

## Evidenziare Piu Regioni

Per aggiungere linee nere tratteggiate su piu regioni del cariogramma puoi usare:

```bash
python LAP_v4.py \
  --prefix /path/to/rfmix/prefix_ \
  --output-dir out_bed \
  --highlight-regions regions_to_highlight.tsv
```

Formato del file:

```text
chrom  start  end  label
chr4   46500000   46600000   locus_A
chr6   25000000   25500000   locus_B
17     43000000   43100000   locus_C
```

La colonna `label` e opzionale. Le regioni FLARE/SNP e i report LAP-VAR restano invariati.

## PDF Nativo Pulito

Per ottenere un PDF piu semplice da modificare in Illustrator/Inkscape/Affinity, puoi evitare la conversione SVG:

```bash
python Plot_LAP_v4.py \
  -I out_bed/SAMPLE.bed \
  -O plots/SAMPLE.native.pdf \
  --pdf-backend native
```

Il backend nativo usa le stesse coordinate e lunghezze hg38 del plot storico, ma disegna direttamente cromosomi, segmenti ancestry, linee SNP, regioni evidenziate, titolo e legenda. La forma esterna riprende la prima sagoma semplice del backend nativo: cromosomi arrotondati con centromero morbido, senza antenne, bande G, effetti 3D o shape sperimentali. Le regioni acrocentriche superiori di chr13/14/15/21/22 sono comunque evidenziate in grigio e riportate in legenda come `Acrocentric p-arms/satellites`. I testi sono esportati come font TrueType editabili quando possibile. Il vecchio backend resta disponibile con `--pdf-backend rsvg`.

## Overlay Centromerico

Per default `Plot_LAP_v4.py` disegna un layer grigio semitrasparente e segmentato verticalmente sulle coordinate centromeriche definite in `CENTROMERE_BP`. Il layer viene disegnato sopra i segmenti ancestry RFMix, ma sotto linee SNP/regioni e label, usando la stessa trasformazione bp-to-pixel del resto del plot. La legenda include una voce `Centromere` con lo stesso pattern segmentato. Nel backend `rsvg`, se lo SVG di base contiene le regioni grigie degli acrocentrici, la legenda include anche `Acrocentric p-arms/satellites`; nel backend `native` questa voce e l'overlay grigio sono disegnati direttamente.

Per disattivarlo:

```bash
python Plot_LAP_v4.py ... --no-centromere-overlay
```

Per regolare colore e trasparenza:

```bash
python Plot_LAP_v4.py ... \
  --centromere-overlay-color "#707070" \
  --centromere-overlay-opacity 0.42
```

## Note sul BED finale

Le righe FLARE aggiunte al BED hanno questo formato:

```text
chrom  start  end  geom_snp  color  haplotype  variant_id  vcf_sample  ancestry_label  ancestry_code  snp_label
```

Esempio:

```text
4  46567563  46567563  geom_snp  #bfa004  1  4:46567563:G:A  COV.COV190_111  East_Asia  1  SNP1
4  46567563  46567563  geom_snp  #a32e2e  2  4:46567563:G:A  COV.COV190_111  Americas   2  SNP1
```
