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
