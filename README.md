

Author: Alessandro Lisi and Michael C. Campbell
(Computational Genomics Lab @ USC, Department of Biological Sciences, Human and Evolutionary Biology section)

In this README file, we present a method for plotting the output of RFMix version 2 (RFMix2). This pipeline will guide you through the process of using RFMix2 to analyze phased VCF files, convert the results to BED format, and visualize the local ancestry using LAP.

1. Prepare Input Dataset

The target file must be in VCF with one individual per VCF. This VCF should contain all chromosomes together and then phased (not gzipped). 

Furthermore, the reference file also should be in VCF and phased (not gzipped). However, in this case, the reference file should be split into individual chromosomes (e.g., 1 through 22 if you are interested in autosomal DNA).


2. Create a Sample Map File

Create a sample map file in *.ref or *.txt tab-delimited format. The map file should include all individuals from the reference dataset (and none of the target individuals). The first column should list the individuals in the same order as in the VCF file. The second column should specify the ancestry based on your own classification. Please see an example of a sample map file in tab-delimited format below:

<pre>
ind1	Africa 
ind2	Africa 
ind3	Africa 
ind4	Europe 
ind5	Europe 
ind6	MiddleEast
ind7	MiddleEast 
ind8	MiddleEast 
ind9	MiddleEast 
</pre>


3. Execute RFMix2

First, we recommend running RFMix2 on one individual at a time (RFmix2 works on one individual at a time). You can use a for loop that includes all target individuals in the analysis, which will generate multiple output files for each individual.

Secondly, the genetic map file should contain all chromosomes together (e.g., 1 through 22 for autosomal DNA). In this file, the columns should be chromosome, position, and cM, respectively, in tab-delimited format. Please see an example of the genetic map file below:

<pre>
chr    pos    cM
chr1   55550  0
chr1   82571  0.080572
chr1   88169  0.092229
chr1   285245 0.439456
chr1   629218 1.478148
chr1   629241 1.478214
</pre>


To run RFMix2, you will need the following files:

-f target VCF/BCF file\
-r reference VCF/BCF file\
-m target map file\
-g genetic map file\
-o output basename\
--chromosome=chromosome to analyze
	
Example Script for Dataset:

Make sure you have the RFmix2 software installed. RFMix2 directly works with phased VCF input files and is significantly faster than RFMix 1.5.4.

To run the RFMix2 software, you can use the command below. If you would like to change the window-size of local ancestry, please refer to the RFMix2 manual.

<pre>
<code>
for i in {1..22}; \
do \
rfmix -e 2 -w 0.5 -f Example_dataset/Target/Mozabite1_allchr.vcf \
-r Example_dataset/Reference/Reference_chr$i.vcf \
-m Sample_map_File/sample_file.txt -g map/all_chr.txt\
-o Output/Mozabite1_chr$i --chromosome=$i; \
done
</code>	
</pre>



Based on your dataset, this process may take from 3 to 30 minutes per chromosome.


4. Combine Output Files

After running RFMix2, four different types of output file are generated for each chromosome; 1) *.Q (global ancestry); 2) *.tsv (marginal probability); 3) *.sis.tsv (condensed information from *.msp.tsv); and 4) *.msp.tsv (crf point). Of these different outputs, we are interested in the *.msp.tsv files. 

If you look at our example output file, Mozabite1_chr2.msp.tsv, you can see at the top that donor populations have a specific number that corresponds to ancestry (e.g., Africa=0, Europe=1, MiddleEast=2).

To combine the *.msp.tsv files for all chromosomes, you can use the following command:


<pre>
<code>
head -n 2 Mozabite1_chr1.msp.tsv > Mozabite1_allchr.msp.tsv 
</code>
</pre>

This above command will preserve the first two rows, which contain header information from the RFMix v.2 output. In addition, the user-specified output filename in the above command must be used in the next command:

<pre>
<code>
for i in {1..22}; do tail -n +3 "Mozabite1_chr$i.msp.tsv" >> Mozabite1_allchr.msp.tsv;done
</code>
</pre>



5. Convert RFMix v. 2 output files to *.bed files.".

In this step, we recommend users acquaint themselves with the usage of this Python script by typing the following command:

python RFMix2ToBed.py --help

The basic usage of RFMix2ToBed.py is:

python RFMix2ToBed.py --input [input filename] –output [output filename]

Here, the input filename must include the *.msp.tsv extension from Step 1. Furthermore, only a prefix name is required in the output filename (i.e. it is not necessary to add the *.bed extension) as follows:

<pre>
<code>
python RFMix2ToBed.py --input Moazbite1_allchr.msp.tsv --output Mozabite1_allchr
</code>
</pre>


This python script accepts the "Mozabite1_allchr.msp.tsv" file as input (please note that you can change the input name and path). This script will generate two *.bed files (in this case, Mozabite1_hap1.bed and Mozabite1_hap2.bed) in which Africa is labeled as ancestry0, Europe as ancestry1, and Middle East as ancestry2. The ancestry labels must match the order in the *.msp.tsv input file.

The RFMix2ToBed.py script will generate two output files, namely Mozabite1_hap1.bed and Mozabite1_hap2.bed.


6. Create the color scheme for ancestry painting along chromosomes.

The two output files from Step 2 will serve as the input files for Step 3.
Prior to running the BedToLap.py.py script, type "python BedToLap.py --help" to see the instructions. You can choose your colors, represented by code notation (e.g., "#0b1b56"), for each ancestry determined by RFMix2. For example,


<pre>
<code>
python BedToLAP.py --bed1 Mozabite1_allchr_hap1.bed --bed2 Mozabite1_allchr_hap2.bed \
--ancestry0 "#0000ff" --ancestry1 "#850b39" --ancestry2 "#F4A500" --out Mozabite1.bed
</code>
</pre>



These color codes indicate that African ancestry is represented by Blue, European ancestry by dark red, and Middle Eastern ancestry by gold.

The output file generated is "Output/Mozabite1.bed"


7. Generate the ancestry plot for each chromosome.


In this step, we suggest users acquaint themselves with the usage of the LAP.py script by typing: python LAP.py --help



<pre>
<code>
python LAP.py -I Mozabite1.bed -O Mozabite1_LAP -B hg38 
</code>
</pre>


where the -O flag specifies the output filename without an extension and -B indicates the genomic build (either “hg37” or “hg38”). The resulting output file (e.g., Mozabite1_LAP.svg and Mozabite1_LAP.pdf) can be edited in Adobe Illustrator or Inkscape. Furthermore, the color quality of the image is 4k resolution (4210 x 1663). 

An SVG file will be generated automatically, which can be edited with the Illustrator software. 

For any questions about this pipeline, please contact Alessandro Lisi by email, alisi@usc.edu

End





