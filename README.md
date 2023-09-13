

Author: Alessandro Lisi, Michael Campbell
(Campbell Computational Genomics LAB @ USC)

In this guide, we present a method for plotting the outputs of RFmix version 2. The process begins with a Pipeline.txt file and follows these steps:\

Prepare Input Dataset:

The target file must be in VCF format and phased (not gzipped) and should contain all chromosomes together, with one individual per VCF.
The reference file must also be in VCF format and phased (not gzipped), but it should be split into chromosomes 1 to 22.
Sample Map File:

Create a sample map file in .ref/.txt format. It should include all individuals from the reference dataset but should not contain the target individuals. The first column should list the individuals in the same order as in the VCF file, and the second column should specify the ancestry based on your own classification (e.g., "Africa").

Example (tab-delimited):

ind1	Africa \
ind2	Africa \
ind3	Africa \
ind4	Europe \
ind5	Europe \
ind6	MiddleEast \
ind7	MiddleEast \
ind8	MiddleEast \
ind9	MiddleEast 


RFMIX2 Execution:

We recommend running RFMIX2 for one individual at a time. You can use a loop to include all target individuals in the analysis. This will generate more output files, one for each individual.

The genetic map file should contain all chromosomes together, from 1 to 22. The columns should represent chromosome, position, and cM, tab-delimited.

Example:

chr	pos	cM \
chr1	55550	0 \
chr1	82571	0.080572 \
chr1	88169	0.092229 \
chr1	285245	0.439456 \
chr1	629218	1.478148 \
chr1	629241	1.478214 \
. \
. \
. \
chr22	45679	0.086453

Summary Needed for RFMIX2:

-f target VCF/BCF file\
-r reference VCF/BCF file\
-m target map file\
-g genetic map file\
-o output basename\
--chromosome=chromosome to analyse
	
	
Example Script for Dataset:

Make sure you have RFmix2 installed. RFmix version 2 directly works with phased VCF input files and is significantly faster than RFmix 1.5.4.


for i in {1..22}; \
do \
rfmix -e 2 -w 0.5 -f Example_dataset/Target/Mozabite1_ind1_allchr.vcf \
-r Example_dataset/Reference/Reference_chr$i.vcf \
-m Sample_map_File/sample_file.txt -g map/all_chr.txt\
-o Output/Mozabite1_ind1_chr$i --chromosome=$i; \
done

Based on your dataset, this process may take from 3 to 30 minutes per chromosome.

Combining Output Files:

Four different output files are generated from RFMIX2 for each chromosome, but we are interested in the .msv.tsp file. 
To combine these files for all chromosomes, you can use the following command:

for i in {1..22}; do tail -n +3 "Mozabite1_ind1_chr$i.msp.tsv"; done > Mozabite1_ind1_allchr.msp.tsv

in Mozabite_ind1_chr2.msp.tsv output, you can see, on the top, 
that each population has a specific number corresponding to ancestry. e.g Africa=0 Europe=1 MiddleEast=2

Running the R Script:

Now, run the R script using the command:

Rscript rfmix2tobed.R

The script accepts the input file "Mozabite_ind1_allchr.msp.tsv," but you can change the input name and path. It will produce two .bed files (e.g., Mozabite1_ind1_hap1.bed and Mozabite1_ind1_hap2.bed), where Africa is labeled as ANC0, Europe as ANC1, and Middle East as ANC2. You can adjust this order or ancestry labels according to your project's needs. The ancestry labels must match the order in the input msp.tsv file.

Using rfmix2bedtotagore.py:

Run the rfmix2bedtotagore.py script with python rfmix2bedtotagore.py --help to see the instructions. You can choose colors, represented in code notation (e.g., #0b1b56), for each ancestry determined in RFMIX.

Example:

python rfmix2bedtotagore.py -1 Output/Mozabite1_ind1_hap1.bed -2 Output/Mozabite1_ind1_hap2.bed \
--anc0 #80cdc1 --anc1 #dfc27d --anc2 #075716 -o Tagore/Mozabite1/Mozabite1_ind1_tagore.bed

These colors indicate, for example, that African ancestry is represented by light blue, European ancestry by light brown, and Middle Eastern ancestry by green.

Plotting with Tagore:

Install Tagore using pip3 install tagore and ensure you have rsvg installed with pip3 install rsvg or brew install rsvg.

Run tagore --help to see how to plot the results. For instance:

run tagore --help. to see how to plot 
e.g

tagore -i Tagore/Mozabite1/Mozabite1_ind1_tagore.bed -p Tagore/Mozabite1/Mozabite1_ind1_tagore -b hg38 -ofmt png

You can choose between PNG or PDF output. An SVG file will also be generated automatically, which can be edited with Illustrator. Note that Tagore does not provide a legend.

End

This pipeline guides you through the process of using RFmix2 to analyze phased VCF files, convert the results to BED format, and visualize the local ancestry using Tagore.



