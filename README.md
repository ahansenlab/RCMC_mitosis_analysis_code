# RCMC_mitosis_analysis_code
This repository contains source code for the article "Dynamics of microcompartment formation at the mitosis-to-G1 transition", used in the analysis of RCMC data.

Code is provided either in the form of Python scripts or as Jupyter notebooks to be run in conda environments containing the required packages. Additionally, genomic positions of microcompartments identified in the paper are included in bedpe format. Many of the analyses described here are reproduced or derived from those published with the 2023 Nature Genetics RCMC paper, "Region Capture Micro-C reveals coalescence of enhancers and promoters into nested microcompartments", for which the source code repository is https://github.com/ahansenlab/RCMC_analysis_code.

## Code summary
### Micro-C alignment (microcprocessing.py)
Required packages:
-	bwamem2
-	samtools
-	sambamba
-	pairtools
-	cooler
-	pairix

Python script used to align reads in .fastq format from paired-end sequencing of Micro-C experiments and produces as output .pairs, .cool and .mcool files compatible with downstream applications such as HiGlass.

Example usage:

```
python /path/to/script/microcbowtie2.py --file_1 pair1.fastq --file_2 pair2.fasq -g mm39 -t 36 -o exampleoutput
```
