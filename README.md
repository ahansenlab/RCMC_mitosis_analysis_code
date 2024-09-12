# RCMC_mitosis_analysis_code
This repository contains source code for the article "Dynamics of microcompartment formation at the mitosis-to-G1 transition", used in the analysis of RCMC data.

Code is provided either in the form of Python scripts or as Jupyter notebooks to be run in conda environments containing the required packages. Additionally, genomic positions of microcompartments identified in the paper are included in bedpe format. Many of the analyses described here are reproduced or derived from those published with the 2023 Nature Genetics RCMC paper, "Region Capture Micro-C reveals coalescence of enhancers and promoters into nested microcompartments", for which the source code repository is https://github.com/ahansenlab/RCMC_analysis_code.

## Code summary
### Genome-wide sequencing alignment (microcprocessing.py)
Required packages:
-	bwa or bowtie2
-	samtools
-	sambamba
-	pairtools
-	cooler
-	pairix

Python script used to align reads in .fastq format from paired-end sequencing of Micro-C & RCMC experiments; produces as output .pairs, .cool and .mcool files compatible with downstream applications such as HiGlass. Sequencing data is aligned genome-wide, after which the unique pairs can be filtered (e.g., using pairtools select) to retain only the reads where both paired ends are within a captured ROI and then merged (e.g., using pairtools merge) to create ROI-only RCMC .pairs and .mcool files.

Use "python /path/to/script/microcprocessing.py -h" to see a description of script functionality. Example usage:

```
python /path/to/script/microcprocessing.py --file_1 pair1.fastq --file_2 pair2.fasq -g mm39 -a bwa -w expand -t 36 -o exampleoutput
```
