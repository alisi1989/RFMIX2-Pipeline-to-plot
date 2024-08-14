Authors: Alessandro Lisi and Michael C. Campbell  
(Human Evolutionary Genomics Lab @ USC, Department of Biological Sciences, Human and Evolutionary Biology section)

In this README file, we present a method for plotting the output of RFMIX version 2 (RFMIX2). This document will guide users through the process of applying RFMIX2 to analyze phased VCF files, converting the results to a *.bed format, and visualizing global and local ancestry using GAP and LAP, respectively.

1. **Prepare Input Dataset**

   The target file must be in VCF and phased with one individual per VCF. This VCF should contain all chromosomes together. 

   Furthermore, the reference file also should be in VCF and phased. However, in this case, the reference file should be split into individual chromosomes (e.g., 1 through 22 if users are interested in autosomal DNA).

2. **Create a Sample Map File**

   Create a sample map file in *.ref or *.txt tab-delimited format. The map file should only contain individuals from the reference dataset (and none of the target individuals). The first column should list the individuals in the same order as they appear in the VCF file. The second column should specify the ancestry name (or classification), which is determined by the user. Please see an example of a sample map file in tab-delimited format below:

   <pre>
   ind1    Africa 
   ind2    Africa 
   ind3    Africa 
   ind4    Europe 
   ind5    Europe 
   ind6    MiddleEast
   ind7    MiddleEast 
   ind8    MiddleEast 
   ind9    MiddleEast 
   </pre>

3. **Execute RFMIX2**

   RFMIX2 analyzes one individual at a time. Consequently, users might consider employing a for loop that includes all target individuals in a given dataset, which will generate multiple output files for each individual.

   Furthermore, the genetic map file should contain all chromosomes together (e.g., 1 through 22 for autosomal DNA). In this file, the columns should be chromosome, position, and cM, respectively, in tab-delimited format. Please see an example of the genetic map file below:

   <pre>
   chr    pos    cM
   chr1   55550  0
   chr1   82571  0.080572
   chr1   88169  0.092229
   chr1   285245 0.439456
   chr1   629218 1.478148
   chr1   629241 1.478214
   </pre>

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

   <pre><code>
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
   </code></pre>

   where variable `i` in a for loop refers to chromosome number (in this case, chromosomes 1 through 22); variable `j` in a for loop refers to the individuals in the dataset (in this case, 1 through 27). 

   It is important to note that the output file name (in this case, `Mozabite${j}_chr${i}` ) must contain the individual name (`Mozabite${j}`) followed by “_chr” and then the chromosome number (`_chr${i}`). 

   After running RFMIX2, four different types of output files are generated for each chromosome: 1) *.Q (global ancestry); 2) *.tsv (marginal probability); 3) *.sis.tsv (condensed information from *.msp.tsv); and 4) *.msp.tsv (crf point). Of these different output files, the *.rfmix.Q and the *.msp.tsv will be the input files in the AncestryGrapher toolkit pipelines. 

**Overview of the AncestryGrapher toolkit**

Indeed, inferences of genetic ancestry are informative for mapping the population origins of genetic risk alleles associated with complex diseases (Cheng et al. 2009; Daya et al. 2014; Freedman et al. 2006) and for understanding the genetic history of admixed populations, including the timing of admixture events (Browning, Waples, and Browning 2023; Daya et al. 2014; Uren, Hoal, and Möller 2020). 

The AncestryGrapher toolkit enables users to visualize global and local ancestry with two distinct pipelines, Global Ancestry Painting (GAP) and Local Ancestry Painting (LAP), that run in a command-line Terminal window on Mac OS X and Linux machines. To execute these pipelines on a Microsoft Windows computer, users will need to install the Anaconda command-line environment on the host machine. 

The AncestryGrapher toolkit can be downloaded to users’ local computers by pressing the “code” button shaded in green on the Github page (https://github.com/alisi1989/RFMIX2-Pipeline-to-plot.git). Using the command-line interface in the Terminal window, users will change the current working directory to the directory where the downloaded “RFMIX2-Pipeline-to-plot-main.zip” folder is located. To unzip this folder, type “unzip RFMIX2-Pipeline-to-plot-main.zip” at the command line prompt (typically indicated by a “$” sign), and the uncompressed “RFMIX2-Pipeline-to-plot-main” folder will appear. To demonstrate the utility of our method, we applied the AncestryGrapher toolkit to the output files from an RFMIX2 analysis of the Finnish (European), Mozabite Berber (North African), and Bedouin (Middle Eastern) populations.  

---

**Global Ancestry Painting (GAP)**

Before proceeding with the pipelines for LAP, users must ensure they have the following Python packages installed:

1. "argparse"
2. "pandas"
3. "matplotlib"
4. "king"
5. "glob"
6. "numpy"
7. "os"
8. "click"

These packages can be installed with pip or pip3:


pip3 install argparse
