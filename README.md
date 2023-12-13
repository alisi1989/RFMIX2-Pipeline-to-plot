

Author: Alessandro Lisi and Michael C. Campbell
(Computational Genomics Lab @ USC, Department of Biological Sciences, Human and Evolutionary Biology section)

In this README file, we present a method for plotting the output of RFMix version 2 (RFMix2). This pipeline will guide you through the process of using RFMix2 to analyze phased VCF files, convert the results to BED format, and visualize the local ancestry using Tagore.

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
rfmix -e 2 -w 0.5 -f Example_dataset/Target/Mozabite1_ind1_allchr.vcf \
-r Example_dataset/Reference/Reference_chr$i.vcf \
-m Sample_map_File/sample_file.txt -g map/all_chr.txt\
-o Output/Mozabite1_ind1_chr$i --chromosome=$i; \
done
</code>	
</pre>



Based on your dataset, this process may take from 3 to 30 minutes per chromosome.


4. Combine Output Files

After running RFMix2, four different types of output file are generated for each chromosome; 1) *.Q (global ancestry); 2) *.tsv (marginal probability); 3) *.sis.tsv (condensed information from *.msp.tsv); and 4) *.msp.tsv (crf point). Of these different outputs, we are interested in the *.msp.tsv files. 

If you look at our example output file, Mozabite1_ind1_chr2.msp.tsv, you can see at the top that donor populations have a specific number that corresponds to ancestry (e.g., Africa=0, Europe=1, MiddleEast=2).

To combine the *.msp.tsv files for all chromosomes, you can use the following command:


<pre>
<code>
for i in {1..22}; do tail -n +3 "Mozabite1_ind1_chr$i.msp.tsv"; done > Mozabite1_ind1_allchr.msp.tsv
</code>
</pre>



5. type "python RFMix2ToBed.py --help to see the opstions".

Next, run the script using the following command:

<pre>
<code>
python scripts/RFMix2ToBed.py --input Output/Mozabite1_ind1_allchr.msp.tsv --output Output/Mozabite1_ind1.bed
</code>
</pre>


This python script accepts the "Mozabite_ind1_allchr.msp.tsv" file as input (please note that you can change the input name and path). This script will generate two *.bed files (in this case, Mozabite1_ind1_hap1.bed and Mozabite1_ind1_hap2.bed) in which Africa is labeled as ancestry0, Europe as ancestry1, and Middle East as ancestry2. You can modify this order or the ancestry labels according to your needs. However, the ancestry labels must match the order in the *.msp.tsv input file.


6. Apply the BedToLap.py script

Prior to running the BedToLap.py.py script, type "python BedToLap.py --help" to see the instructions. You can choose your colors, represented by code notation (e.g., "#0b1b56"), for each ancestry determined by RFMix2. For example,


<pre>
<code>
python scripts/BedToLap.py -bed1 Output/Mozabite1_ind1_hap1.bed -bed2 Output/Mozabite1_ind1_hap2.bed \
--ancestry0 #80cdc1 --ancestry1 #dfc27d --ancestry2 #075716 -o Output/Mozabite1_ind1_LAP.bed
</code>
</pre>



These color codes indicate that African ancestry is represented by light blue, European ancestry by light brown, and Middle Eastern ancestry by green.

The output file generated is "Output/Mozabite1_ind1_LAP.bed"


7. Plot with LAP


Run the command "python LAP.py --help" to see how to plot the results. For instance, you can run the following command to generate a plot in png format:

<pre>
<code>
python scripts/LAP.py -I Output/Mozabite1_ind1_LAP.bed -O plot/Mozabite1_ind1_tagore -B hg38 
</code>
</pre>



Here, "plot/Mozabite1_ind1_LAP.bed" is the file generated in Step 6 above; -O is the prefix of the output file; -B is the genome reference build

An SVG file will be generated automatically, which can be edited with the Illustrator software. Please note that LAP.py does not provide a legend for plots. Use LAP_legend.py instead.

For any questions about this pipeline, please contact Alessandro Lisi by email, alisi@usc.edu

End





