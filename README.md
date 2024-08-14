

Authors: Alessandro Lisi and Michael C. Campbell
(Human Evolutionary Genomics Lab @ USC, Department of Biological Sciences, Human and Evolutionary Biology section)

In this README file, we present a method for plotting the output of RFMIX version 2 (RFMIX2). This document will guide users through the process of applying RFMIX2 to analyze phased VCF files, converting the results to a *.bed format, and visualizing global and local ancestry using GAP and LAP, respectively.

1. Prepare Input Dataset

The target file must be in VCF and phased with one individual per VCF. This VCF should contain all chromosomes together. 

Furthermore, the reference file also should be in VCF and phased. However, in this case, the reference file should be split into individual chromosomes (e.g., 1 through 22 if users are interested in autosomal DNA).

2. Create a Sample Map File

Create a sample map file in *.ref or *.txt tab-delimited format. The map file should only contain individuals from the reference dataset (and none of the target individuals). The first column should list the individuals in the same order as they appear in the VCF file. The second column should specify the ancestry name (or classification), which is determined by the user. Please see an example of a sample map file in tab-delimited format below:

<pre>
<code>
ind1	Africa 
ind2	Africa 
ind3	Africa 
ind4	Europe 
ind5	Europe 
ind6	MiddleEast
ind7	MiddleEast 
ind8	MiddleEast 
ind9	MiddleEast 
<code>
<pre>



3. Execute RFMIX2

RFMIX2 analyzes one individual at a time. Consequently, users might consider employing a for loop that includes all target individuals in a given dataset, which will generate multiple output files for each individual.

Furthermore, the genetic map file should contain all chromosomes together (e.g., 1 through 22 for autosomal DNA). In this file, the columns should be chromosome, position, and cM, respectively, in tab-delimited format. Please see an example of the genetic map file below:

<pre>
<code>
chr    pos    cM
chr1   55550  0
chr1   82571  0.080572
chr1   88169  0.092229
chr1   285245 0.439456
chr1   629218 1.478148
chr1   629241 1.478214
<code>
<pre>



To run RFMIX2, users will need the following:

-f target VCF/BCF file \ # vcf/bcf file containing data from the target population(s)
-r reference VCF/BCF file \ # vcf file containing data from the reference population(s)
-m target map file \ # file containing genetic map for SNP loci in the target population
-g genetic map file \ # file containing genetic map for SNP loci in the reference population
-o output basename \ # the prefix of an output file name (without an extension)
--chromosome #chromosome to analyze

The genetic map for each chromosome is provided in the “Genetic_Map” folder in the “RFMIX2-Pipeline-to-plot-main” directory.

Example of basic usage:

To run the RFMIX2 software, users can specify the command below. If users wish to change the window size of local ancestry, we recommend they refer to the RFMIX2 manual.

<pre>
<code>
for i in {1..22}; do
    for j in {1..27}; do
        rfmix -f Example_Dataset/Target/Mozabite_${j}.vcf.gz \
              -r Example_Dataset/Reference/Reference_Phased_chr${i}.vcf.gz \
              -m Sample_map_File/Sample_Reference.txt \
              -g Genetic_Map/chr${i}.b38.txt \
              -o Example_Dataset/RFMIX2_Output/Mozabite${j}_chr${i} \
              --chromosome=${i}
    done
done
<code>
<pre>

where variable i in a for loop refers to chromosome number (in this case, chromosomes 1 through 22); variable j in a for loop refers to the individuals in the dataset (in this case, 1 through 27). 

It is important to note that the output file name (in this case, Mozabite${j}_chr${i} ) must contain the individual name (Mozabite${j}) followed by “_chr” and then the chromosome number (_chr${i}). 

After running RFMIX2, four different types of output files are generated for each chromosome: 1) *.Q (global ancestry); 2) *.tsv (marginal probability); 3) *.sis.tsv (condensed information from *.msp.tsv); and 4) *.msp.tsv (crf point). Of these different output files, the *.rfmix.Q  and the *.msp.tsv will be the input files in the AncestryGrapher toolkit pipelines. 

Overview of the AncestryGrapher toolkit 

Indeed, inferences of genetic ancestry are informative for mapping the population origins of genetic risk alleles associated with complex diseases (Cheng et al. 2009; Daya et al. 2014; Freedman et al. 2006) and for understanding the genetic history of admixed populations, including the timing of admixture events (Browning, Waples, and Browning 2023; Daya et al. 2014; Uren, Hoal, and Möller 2020). 

The AncestryGrapher toolkit enables users to visualize global and local ancestry with two distinct pipelines, Global Ancestry Painting (GAP) and Local Ancestry Painting (LAP), that run in a command-line Terminal window on Mac OS X and Linux machines. To execute these pipelines on a Microsoft Windows computer, users will need to install the Anaconda command-line environment on the host machine. 

The AncestryGrapher toolkit can be downloaded to users’ local computers by pressing the “code” button shaded in green on the Github page (https://github.com/alisi1989/RFMIX2-Pipeline-to-plot.git). Using the command-line interface in the Terminal window, users will change the current working directory to the directory where the downloaded “RFMIX2-Pipeline-to-plot-main.zip” folder is located. To unzip this folder, type “unzip RFMIX2-Pipeline-to-plot-main.zip” at the command line prompt (typically indicated by a “$” sign), and the uncompressed “RFMIX2-Pipeline-to-plot-main” folder will appear. To demonstrate the utility of our method, we applied the AncestryGrapher toolkit to the output files from an RFMIX2 analysis of the Finnish (European), Mozabite Berber (North African), and Bedouin (Middle Eastern) populations.  


##################Global Ancestry Painting (GAP)########################

Before proceeding with the pipelines for LAP, users must ensure they have the following Python packages installed:
1.	"argparse"
2.	"pandas"
3.	"matplotlib"
4.	"king"
5.	"glob"
6.	"numpy"
7.	"os"
8.	"click"

These packages can be installed with pip or pip3

e.g., pip3 install arparse

The GAP pipeline consists of three separate Python scripts: 1) RFMIX2ToBed4GAP.py; 2) BedToGAP.py (this script creates the input file for GAP); 3) GAP.py. Users will need to change their working directory to “GAP” in the “RFMIX2-Pipeline-to-plot-main” directory (e.g., cd RFMIX2-Pipeline-to-plot-main/GAP/)

Step 1: Combine the RFMIX2 output files for all the chromosomes per individual into a single file and merge all the individuals together.

The RFMIX2 software generates a *.rfmix.Q (global ancestry information) output file for each chromosome per individual.  Users must ensure that the output file names from RFMIX2 include the name of the individual (e.g., Mozabite1) followed by “_chr” and then the chromosome number (e.g., _chr2 for chromosome 2). This entire name or prefix must appear before the *.rfmix.Q extension (e.g., Mozabite1_chr2.rfmix.Q)._

To combine chromosomes per individual into a single file, users will execute the Python script below:

Basic command line:

<pre>
<code>
python RFMIX2ToBed4GAP.py --prefix [argument] --chr [argument] --output [argument] --sort-ancestry [argument] 
</code>
</pre>

where users need to enter: 1) the prefix of the RFMIX2 output filename (without the file extension) after the "--prefix" flag; 2) chromosome number placed between curly braces, {}, after the "--chr" flag; and 3) the prefix of an output filename (without a file extension) after the "--output" flag. This script will automatically generate a *.bed file.

“--sort-ancestry” is an optional flag that can be used to sort by ancestry. Importantly, the population ancestry name provided after the “--sort-ancestry” flag must be identical to the population ancestry name specified in the header of the *.rfmix.Q output file. The script will then sort this ancestry from the largest ancestry proportion to the smallest.
 
Example of basic usage:

<pre>
<code>
python Scripts/RFMIX2ToBed4GAP.py --prefix ../Example_Dataset/RFMIX2_Output/Mozabite --chr {1..22} --output Output_GAP
</code>
</pre>

Example of basic usage for target chromosomes:

<pre>
<code>
python Scripts/RFMIX2ToBed4GAP.py --prefix ../Example_Dataset/RFMIX2_Output/Mozabite --chr 2 5 7 --output Output_GAP
</code>
</pre>

Example of usage with “--sort-ancestry” flag:

<pre>
<code>
python Scripts/RFMIX2ToBed4GAP.py --prefix ../Example_Dataset/RFMIX2_Output/Mozabite --chr {1..22} --output Output_GAP/ --sort-ancestry Middle_East
</code>
</pre>




Step 2: Create the input file for GAP to visualize the global ancestry proportions. 

The output file generated in Step 1 will serve as the input file for Step 2. However, there are additional options that users may wish to consider. These options can be accessed by typing the following command:

python BedToGAP.py --help

Basic command line:

<pre>
<code>
python BedToGAP.py --input [argument] --ancestry [argument] --out [argument] 
</code>
</pre>

where the “--input” flag accepts the file name from Step 1 (*.bed), and the “--out” flag accepts the user-specified output file name with the *.bed extension. If specified, the --ancestry flag requires two arguments: 1) population ancestry name, and 2) the hex color code, which can be found on any website on the internet (e.g., computerhope.com). It is important to note, however, that the population ancestry name provided after the “--ancestry” flag must be identical to the population ancestry name in the header of the output file from Step 1. By default, this script can assign up to ten distinct colors, one for each ancestry component. Alternatively, users can assign their own colors to ancestry components using the “--ancestry” flags.

Example of basic usage (with default ancestry colors):

<pre>
<code>
python Scripts/BedToGAP.py --input Output_GAP/Mozabite.bed --out Output_GAP/Mozabite_GAP.bed
</code>
</pre>

Example of usage with “--ancestry” flag  (for user-specified colors):

<pre>
<code>
python Scripts/BedToGAP.py --input Output_GAP/Mozabite.bed --ancestry0 Africa #a38905 --ancestry1 Europe #a30d05 --ancestry2 Middle_East #0e6b05 --out Output_GAP/Mozabite_GAP.bed
</code>
</pre>

If specified, the --ancestry flag requires two arguments: 1) population ancestry name, and 2) the hex color code. Again, the population ancestry must be identical to the population ancestry name present in the header of the output file from RFMIX2. 

Step 3: Generate the plot for global ancestry proportion with GAP.py.

In this step, GAP.py requires a single input file name and a user-specified output file name as arguments:

python GAP.py --input [argument] --output [argument] 

where users will enter the output file name from Step 2 as an argument for the “--input" flag.  Furthermore, users also must specify the output filename, adding either a “.pdf” or a “.svg” extension to the end. This script will generate an output file with ancestry proportions for each individual in a bar plot in either “pdf” or “svg” format. 

Examples of usage:

<pre>
<code>
python Scripts/GAP.py --input Output_GAP/Mozabite_GAP.bed --output Output_GAP/Mozabite_GAP.pdf
</code>
</pre>

<pre>
<code>
python Scripts/GAP.py --input Output_GAP/Mozabite_GAP.bed --output Output_GAP/Mozabite_GAP.svg
</code>
</pre>


###############Local Ancestry Painting (LAP)######################

Before proceeding with the pipelines for LAP, users must ensure they have the following Python packages installed:
9.	"argparse"
10.	"pandas"
11.	"matplotlib"
12.	"king"
13.	"glob"
14.	"numpy"
15.	"os"
16.	"click"

These packages can be installed with pip or pip3

e.g., pip3 install arparse

Users will also need to install a library (“librsvg”) to create, edit, and convert pdf and svg files. This library can be installed on a MacOS machine with the following command typed in a Terminal window: brew install rsvg).

Alternatively, this library can be downloaded manually from:  https://download.gnome.org/sources/librsvg. 

In addition, the “librsvg” library can be installed on a Linux/Ubuntu/Debian/ machine with the following command typed in a Terminal window:

sudo apt-get install -y librsvg2-dev 

Alternatively, this library can be downloaded manually from https://manpages.ubuntu.com/manpages/trusty/man1/rsvg-convert.1.html

The LAP pipeline consists of three separate Python scripts: 1) RFMIX2ToBed.py; 2) BedToLAP.py; 3) LAP.py. Users will need to change their working directory to “LAP” in the “RFMIX2-Pipeline-to-plot-main” folder (e.g., cd RFMIX2-Pipeline-to-plot-main/LAP/).

Step 1: Combine the output files from RFMIX2 into a single file and generate *.bed input files with RFMIX2ToBed.py.

RFMIX2 generates a *.msp.tsv (local ancestry information) output file for each chromosome for a given individual. Users must ensure that the output file names from RFMIX2 include the name of the individual (e.g., Mozabite1) followed by “_chr” and then the chromosome number (e.g., _chr2 for chromosome 2). This entire name or prefix must precede the *. msp.tsv extension (e.g., Mozabite1_chr2.msp.tsv). 

In this step, we recommend that users acquaint themselves with the usage of this Python script by typing the following command:

<pre>
<code>
python Scripts/RFMIX2ToBed.py --help
</code>
</pre>

As stated above, RFMIX2 generates a *.msp.tsv (local ancestry information) output file for each chromosome for a given individual. Users will combine these chromosomes into a single file for each individual, which will contain header and ancestry information for all chromosomes for each individual in the dataset.

Basic command line:

<pre>
<code>
python RFMIX2ToBed.py --prefix [argument] --chr [argument ] --output [argument]
</code>
</pre>

where users are required to: 1) enter the prefix of the RFMIX2 output filename (without the file extension) after the "--prefix" flag; 2) enter a range of chromosome numbers placed between curly braces, {}, after the "--chr" flag; and 3) provide the pathway to where the output files will be saved. File names for the output files are NOT required.  

The RFMIX2ToBed.py script will automatically generate two *.bed files (*_hap1.bed and *_hap2.bed) for each individual with the same prefix. The *_hap1.bed and *_hap2.bed files correspond to diploid chromosomes (i.e., the maternal and paternal copies of chromosomes)

Example of usage:

<pre>
<code>
python Scripts/RFMIX2ToBed.py --prefix ../Example_Dataset/RFMIX2_Output/Mozabite --chr {1..22} --output Output_LAP/
</code>
</pre>

Alternatively, if users require a subset of chromosomes in a single file, they can run the following command: 

<pre>
<code>
python Scripts/RFMIX2ToBed.py --prefix ../Example_Dataset/RFMIX2_Output/Mozabite --chr 2 3 5 7 --output Output_LAP/
</code>
</pre>

In this scenario, specific chromosome numbers, separated by spaces, will appear after the “--chr" flag.

In either example,  the RFMIX2ToBed.py script will generate two output files for each individual, namely Mozabite1_hap1.bed and Mozabite1_hap2.bed.

Step 2: Create the color scheme for ancestry painting along chromosomes with BedToLAP.py.

In this step, we recommend that users acquaint themselves with the usage of this Python script by typing the following command in the Terminal window:

<pre>
<code>
python BedToLAP.py --help
</code>
</pre>

Basic command line:

<pre>
<code>
python BedToLAP.py --bed1 [argument] --bed2 [argument] --out [argument]
</code>
</pre>

where users will provide; 1) *_hap1.bed filename after the “–-bed1” flag (from Step 1); 2) the *_hap2.bed filename after the “--bed2” flag (from Step 1); and 3) a user-specified output filename with the *.bed extension after the “--out" flag. The *.bed extension must be added to the output file name; otherwise, the output file cannot be used in the next step (Step 3).

Example of usage:

<pre>
<code>
for i in {1..27}; do python Scripts/BedToLAP.py --bed1 Output_LAP/Mozabite${i}_hap1.bed --bed2 Output_LAP/Mozabite${i}_hap2.bed --out Output_LAP/Mozabite${i}.bed; done
</code>
</pre>

where variable i in a for loop refers to the individuals in the dataset (1 through 27). This script will automatically assign a default color to each ancestry component (default colors for a maximum of ten ancestry components will be generated). However, users can also choose up to ten distinct colors, one for each ancestry component with the “--ancestry" flag (please use "python BedToLap.py --help" to see the options). The “--ancestry" flag requires two single arguments: 1) population ancestry name, and 2) the hex color code, which can be found on any website on the internet (e.g., computerhope.com). Again, the population ancestry must be identical to the population ancestry name present in the header of the output file from RFMIX2. 

In addition, to highlight a specific gene or genomic region on a chromosome, users can run the same command as above with additional parameters. Specifically,

<pre>
<code>
for i in {1..27}; do python Scripts/BedToLAP.py --bed1 Output_LAP/Mozabite${i}_hap1.bed --bed2 Output_LAP/Mozabite${i}_hap2.bed --out Output_LAP/Mozabite${i}.bed --chr 2 --from-bp 135787850 --to-bp 155837184; done
</code>
</pre>

where variable i in a for loop refers to the individuals (1 through 27) in the dataset; “--chr" accepts a chromosome number as an argument; --from-bp and --to-bp flags require the start and end positions of a gene or genomic region of interest in base pairs, respectively; a dashed black line, indicating a gene or genomic region of interest, will be added to the final plot generated in Step 3.

Step 3: Generate the ancestry plot for each chromosome with LAP.py.

In this step, we suggest that users acquaint themselves with the usage of the LAP.py script by typing: 

python LAP.py --help

Basic command line:

<pre>
<code>
python LAP.py -I [argument] -O [argument] -B [argument] 
</code>
</pre>

Example of usage:

<pre>
<code>
for i in {1..27}; do python Scripts/LAP.py -I Output_LAP/Mozabite${i}.bed -B hg38 -O Output_LAP/Mozabite${i}.pdf; done
</code>
</pre>

where variable i in the for loop refers to the individuals (1 through 27); the -I flag specifies the input file from Step 3; the -O flag requires a user-specified output filename with either a .pdf or .svg extension; and -B indicates the genomic build (either “hg37” or “hg38”). The output files will be in a high-quality editable “.pdf” format.

Regardless, the resulting output file (e.g., Mozabite1_LAP.pdf) will contain ancestry-informative karyograms along with a legend of ancestry origin and the names of individuals in the dataset. Furthermore, the images in the output files will have 4k resolution (4210 x 1663) and can be edited in Adobe Illustrator or Inkscape. 

For any questions about this pipeline, please contact Alessandro Lisi (alisi@usc.edu) or Michael Campbell (mc44680@usc.edu) by email.

End





