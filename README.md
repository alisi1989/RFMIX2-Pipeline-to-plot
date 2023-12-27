

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
After running RFMix2, four different types of output file are generated for each chromosome; 1) *.Q (global ancestry); 2) *.tsv (marginal probability); 3) *.sis.tsv (condensed information from *.msp.tsv); and 4) *.msp.tsv (crf point). Of these different outputs, we are interested in the *.msp.tsv files. 

Overview of the AncestryGrapher toolkit 

Unequivocally, inferences of genetic ancestry are informative for mapping the population origins of genetic risk alleles associated with complex diseases (Cheng et al. 2009; Daya et al. 2014; Freedman et al. 2006) and for understanding the genetic history of admixed populations, including the timing of admixture events (Browning, Waples, and Browning 2023; Daya et al. 2014; Uren, Hoal, and Möller 2020). 
The AncestryGrapher toolkit enables users to visualize global and local ancestry with two distinct pipelines, Global Ancestry Painting (GAP) and Local Ancestry Painting (LAP; Figure 1), that run in a command-line Terminal window on Mac OS X and Linux machines. To execute these pipelines on a Microsoft Windows computer, users will need to install the Anaconda command-line environment on the host machine. 
The AncestryGrapher toolkit can be downloaded to users’ local computers by pressing the “code” button shaded in green on the Github page (https://github.com/alisi1989/RFmix2-Pipeline-to-plot.git). Using the command-line interface in the Terminal window, users will change the current working directory to the directory where the downloaded “RFmix2-Pipeline-to-plot-main.zip” folder is located. To unzip this folder, type “unzip RFmix2-Pipeline-to-plot-main.zip” at the command line prompt (typically indicated by a “$” sign) and the uncompressed “RFmix2-Pipeline-to-plot-main” folder will appear. To demonstrate the utility of our method, we applied the AncestryGrapher toolkit to the output files from an RFMix v.2 analysis of the Mozabite Berber population from North Africa. 


Global Ancestry Painting (GAP)

The GAP pipeline consists of four separate Python scripts: 1) individuals_collapse.py; 2) RFMix2ToBed4GAP.py; 3) BedToGAP.py (this script creates the input file for GAP); 4) GAP.py. Users will need to change their working directory to “GlobalAncestryPaint” inside of the “RFmix2-Pipeline-to-plot-main” directory.

Step 1: Combine the output files for all the chromosomes per individual into a single file.

The RFMix v.2 software generates a *.rfmix.Q (Global ancestry information) output file for each chromosome. To combine chromosomes per individual into a single file, users will execute the following the Python script:  


<pre>
<code>
python individuals_collapse.py
</code>
</pre>

However, before running this script, users will need to manually change the sample name, the name of the input file, and the name of the output file. More specifically, users should modify the output file name in line 4, the name of the input file in line 14, and then the sample name in lines 17 through 20. This process should be repeated for each individual in the dataset. The script also calculates the average of each ancestry component across all chromosomes in each individual. For example, if there are ten individuals in a dataset, this script will generate ten separate output files with ancestry information for a given individual. Alternatively, if users wish to analyze a single chromosome only, they can proceed directly to Step 2.

Step 2: Merge all the individuals or all the chromosomes into one file. 

To combine the separate files from Step 1 into a single file, users will run RFMix2ToBed4GAP.py as follows: 


<pre>
<code>
python RFMix2ToBed4GAP.py --input Mozabite{1..27}_allchr.rfmix.Q --output Mozabite_allind_allchr.bed
</code>
</pre>

where {1..27} signifies the number of individuals in a given dataset; specifically, {1..27} indicates 1 to 27 individuals to combine into a single file. Likewise, users will provide a range for the number of individuals they wish to combine beginning with “1”. If users only have a single individual in their dataset, a range is not required. Finally, users must add the .bed extension to the end of the output filename (e.g., Mozabite_allind_allchr.bed).

Step 3: Create the input file for GAP to visualize the global ancestry proportions. 

The output file generated in Step 2 will serve as the input file for Step 3. However, there are additional options that users should consider, which can be accessed by typing the following command:

python BedToGAP.py --help

Notably, users can incorporate up to five distinct colors, one for each ancestry component, with the --ancestry flag:


<pre>
<code>
python BedToGAP.py --input Mozabite_allind_allchr.bed --ancestry1 Africa "#0000ff" --ancestry2 Europe "#850b39" --ancestry3 Middle_East "#F4A500" --out Mozabite_GAP.bed
</code>
</pre>

where the --ancestry flag accepts two arguments: population ancestry name and the hex color code, which can be found on any website on the internet (e.g., computerhope.com). It is important to note, however, that the population ancestry name provided after the --ancestry flag must be identical to the population ancestry name in the header of the output file from Step 2. In addition, the *.bed extension must be specified in the output filename (e.g., Mozabite1_GAP.bed). 

Step 4: Generate the plot for global ancestry proportion.

In this step, GAP.py accepts a single input file as follows:

<pre>
<code>
python GAP.py --input Mozabite_GAP.bed --output Mozabite_GAP.pdf
</code>
</pre>


Users also must specify the output filename, adding either a “.pdf” or a “.svg” extension to the end. This script will generate an output file with ancestry proportions for each individual in a bar plot in either “pdf” or “svg” format (Figure 1). 


Local Ancestry Painting (LAP)

The LAP pipeline consists of three separate Python scripts: 1) RFMIX2ToBed.py; 2) BedToLAP.py; 3) LAP.py. 

The working directory for this pipeline will be “LocalAncestryPaint”, and all output files should be saved within this directory.

Step 1: Combine the output files from RFMix v.2 into a single file.

RFMix v.2 generates the *.msp.tsv output file for each chromosome. To combine these chromosomes into a single file, use the combine command line in Github as follows: 

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


The resulting output file will contain the header and the ancestry information for all the autosomal chromosomes.

Alternatively, if users require a subset of chromosomes in a single file, they can run the following command: 

head -n 2 Mozabite1_chr1.msp.tsv > Mozabite1_chr5_7_12.msp.tsv

for i in 5 7 12; do tail -n +3 "Mozabite1_chr${i}.msp.tsv" >> Mozabite1_chr5_7_12.msp.tsv;done

where 5 7 12 correspond to chromosomes 5, 7, and 12, respectively.


Step 2: Convert RFMix v.2 output files to *.bed input files.

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


The RFMix2ToBed.py script will generate two output files, namely Mozabite1_hap1.bed and Mozabite1_hap2.bed.

Step 3: Create the color scheme for ancestry painting along chromosomes.

The two output files from Step 2 will serve as the input files for Step 3:


<pre>
<code>
python BedToLAP.py –bed1 Mozabite1_allchr_hap1.bed –bed2 Mozabite1_allchr_hap2.bed --ancestry0 "#0000ff" --ancestry1 "#850b39" --ancestry2 "#F4A500" --out Mozabite1.bed
</code>
</pre>


In addition, users can enter up to five distinct colors (color codes), one for each ancestry component using the –ancestry flag. This flag requires a single argument, specifically a user-defined hex color code, which can be found on any website on the internet (e.g., computerhope.com). If no color is specified, default colors—for a maximum of four ancestry components—will be selected. More explicitly, the default color for ancestry0 is blue; magenta for ancestry1; gold for ancestry2 and dark gold for ancestry3. Finally, the *.bed extension must be added to the output file otherwise the output file cannot be used in the next step.

Alternatively, to highlight a specific gene or genomic region on a chromosome, users can run the same command as above with additional parameters. Specifically,


<pre>
<code>
python BedToLAP.py --bed Mozabite1_hap1.bed --bed Mozabite1_hap2.bed –ancestry0 "#0000ff"--ancestry1 "#850b39" --ancestry2 "#F4A500" --chr 2 --from-bp 135787850 --to-bp 135837184 –feature-type triangle --out Mozabite1.bed
</code>
</pre>



where --chr specifies the chromosome on a given the gene or genomic region of interest is located; --from-bp and --to-bp indicate the start and end positions in base pairs of the gene or genomic region of interest; the --feature-type flag specifies either a dashed line or arrow to indicate a gene or genomic region of interest.

For more detailed information about the different options in this script, users can type the following command in the Terminal window:

python BedToLAP.py --help

Step 4: Generate the ancestry plot for each chromosome.

In this step, we suggest users acquaint themselves with the usage of the LAP.py script by typing: 

python LAP.py --help

An example of the basic usage of the LAP.py script is as follows:

python LAP.py -I Mozabite1.bed -O Mozabite1_LAP -B hg38 --ancestry0 [legend name] ["color code"] --output [output filename] 

where the -I specifies the input file from Step 3; -O flag specifies the output filename without an extension; and -B indicates the genomic build (either “hg37” or “hg38”). Furthermore, users can enter up to five different legend names and corresponding colors with the --ancestry flag; however, the user-defined colors must be identical to the colors selected in Step 3. The output files will be in high-quality “svg” and “pdf” formats. For example,



<pre>
<code>
python LAP.py -I Mozabite1.bed -O Mozabite1_LAP -B hg38 --ancestry0 Africa "#0000ff" --ancestry1 Europe "#850b39" --ancestry2 Middle_East "#F4A500" 
</code>
</pre>

The resulting output files (e.g., Mozabite1_LAP.svg and Mozabite1_LAP.pdf) will contain ancestry-informative karyotypes along with a legend of ancestry origin and the name of the individual. Moreover, the images in the output files will have 4k resolution (4210 x 1663) and can be edited in Adobe Illustrator or Inkscape. 


For any questions about this pipeline, please contact Alessandro Lisi by email, alisi@usc.edu

End





