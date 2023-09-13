# RFmix2-Pipeline-to-plot
Here we present a method to plot the outputs of RFmix version 2

starts from Pipeline.txt and follow up

Author "Alessandro Lisi, Michael Campbell"
"Campbell Computational Genomics LAB @USC")

#Pipeline for RFmix 2 and plot with Tagore all thr chrs

#Prepare input dataset

#The target file must to be in vcf and phased(not gzip) all chrs together. one individual per vcf

#The reference file must to be in vcf and phased (not gzip) split in chrs 1 to 22. 

In this script Mozabite is used as Target. 
African, Europe and MiddleEast is used as Reference. 
For this simulation analysis, we chose the Mozabites because 
they are a population highly amixed with European, Middle Eastern and Africans

Prepare the samplefile 
the sample map file in .ref/.txt must contain all the individuals of the dataset used as reference. 
it must not contain the target individuals. the first column are the individual in the same order 
of the vcf file. the second column must contain the ancestry based on your own. 

e.g tab delimited

ind1	Africa \
ind2	Africa \
ind3	Africa \
ind4	Europe \ 
ind5	Europe \
ind6	MiddleEast \
ind7	MiddleEast \
ind8	MiddleEast \
ind9	MiddleEast \


We suggest to run RFMIX2 just on individual at a time. you can use loop to include all individuals
of target into the analysis. you will obtein more output for each individual.

The genetic map file must to be contain all chromosomes together. from 1 to 22. the first column is chr, 
second column is pos and the tirtdh column is cM. is tab delimited 

chr	pos	cM
chr1	55550	0
chr1	82571	0.080572
chr1	88169	0.092229
chr1	285245	0.439456
chr1	629218	1.478148
chr1	629241	1.478214
.
.
.
chr22	45679	0.086453

summary needed for RFMIX2


-f <target VCF/BCF file>
	-r <reference VCF/BCF file>
	-m <target map file>
	-g <genetic map file> 
	-o <output basename>
	--chromosome=<chromosome to analyze>
	
	
script for example dataset. is in Example_Dataset folder

make sure you have RFmix2 installed. is crucial to use RFmix version 2 
( "conda install rfmix" ) or "https://github.com/slowkoni/rfmix" 
Why RFmix2? compare to Rfmix 1.5.4 Rfmix2 works directly with vcf phased input file and is 20x faster than rfmix1.5.4

for i in {1..22}; \
do \
rfmix -e 2 -w 0.5 -f Example_dataset/Target/Mozabite1_ind1_allchr.vcf \
-r Example_dataset/Reference/Reference_chr$i.vcf \
-m Sample_map_File/sample_file.txt -g map/all_chr.txt\
-o Output/Mozabite1_ind1_chr$i --chromosome=$i; \
done

based on your dataset takes from 3 to 30 minutes per chromosomes

4 different output are generated from RFMIX2 for each chromosome. we are looking just on .msv.tsp one.

condense all the msp.tsv together without the first 2 line for the other files
e.g: 

for i in {1..22}; do tail -n +3 "Mozabite1_ind1_chr$i.msp.tsv"; done > Mozabite1_ind1_allchr.msp.tsv

in Mozabite_ind1_chr2.msp.tsv output, you can see, on the top, 
that each population has a specific number corresponding to ancestry. e.g Africa=0 Europe=1 MiddleEast=2

now run the R script. with 

Rscript rfmix2tobed.R

the script accepts the input file " Mozabite_ind1_allchr.msp.tsv" you can change the input name and path
you will obetin 2 .bed file. e.g Mozabite1_ind1_hap1.bed Mozabite1_ind1_hap2.bed 
(where Africa is ANC0 e Europe is ANC1 e MiddleEast is ANC2). you can change 
this order or ancestry based on your project. the ancestry needs to be in the same order of input msp.tsv

now run rfmix2bedtotagore.py with "python rfmix2bedtotagore.py --help". see the istructions
for all ancestry that you decided in RFMIX you can choose a color with code notation. e.g (#0b1b56)
e.g

python rfmix2bedtotagore.py -1 Output/Mozabite1_ind1_hap1.bed -2 Output/Mozabite1_ind1_hap2.bed \
--anc0 #80cdc1 --anc1 #dfc27d --anc2 #075716 -o Tagore/Mozabite1/Mozabite1_ind1_tagore.bed

using these color meens that African ancestry is light blue, light brown is Europe and green is Middle_east

Plot with Tagore
https://github.com/jordanlab/tagore#installation
now install tagore with : pip3 install tagore
make sure you have rsvg installed : pip3 install rsvg or brew install rsvg

run tagore --help. to see how to plot 
e.g

tagore -i Tagore/Mozabite1/Mozabite1_ind1_tagore.bed -p Tagore/Mozabite1/Mozabite1_ind1_tagore -b hg38 -ofmt png
you can choose png or pdf. automaticaly .svg is generated. is editable with Illustrator. 
Legend is not provided by Tagore

end



