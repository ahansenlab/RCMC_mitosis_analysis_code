# RCMC_mitosis_analysis_code
This repository contains source code for the article "Dynamics of microcompartment formation at the mitosis-to-G1 transition", used in the analysis of RCMC data.

Code is provided either in the form of Python / R scripts or as Jupyter notebooks to be run in conda environments containing the required packages. Additionally, genomic positions of microcompartments identified in the paper are included in bedpe format. Many of the analyses described here are reproduced or derived from those published with the 2023 Nature Genetics RCMC paper, "Region Capture Micro-C reveals coalescence of enhancers and promoters into nested microcompartments", for which the source code repository is https://github.com/ahansenlab/RCMC_analysis_code.

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

### Visualization of contact maps and genomic tracks (ContactMaps&GenomicTracks.ipynb)
Required packages:
-	cooltools
-	cooler
-	coolbox
-	matplotlib

Jupyter notebook used to generate visualizations of contact maps and genomics tracks for figures. Contact map visualization is accomplished using cooltools and requires a .mcool file of contacts from RCMC or a comparable method. Genomic track visualization is accomplished using coolbox and requires a .mcool file of contacts, gene annotations (.gtf format or similar), and ChIP-seq, RNA-seq, and ATAC-seq datasets (.bw format).

### Lists of probes used for capturing regions of interest (MitosisRCMC_CaptureProbes_80mers_mm39.bed)

BED format file listing the genomic locations of all probes used for capturing the _Id1_ (chr2), _Klf1_ (chr8), _Cdt1_ (chr8), _Dag1_ (chr9), and _Myc_ (chr15) regions used in the RCMC experiments. Coordinates are provided for the mm39 reference genome. Columns in the file are as follows: the first is the chromosome the region is located on, the second is the start coordinate of the probe range, and the third is the end coordinate of the probe range, with probe ranges being in multiples of 80 representing one or more consecutive 80mer probes.

### List of manually-annotated microcompartment loops (MitosisRCMC_LoopCalls_mm39.bed & MitosisRCMC_LoopCalls_PlusMin1kb_mm39.bedpe) and anchors (MitosisRCMC_AnchorCalls_PlusMin1kb_mm39.bed)

BEDPE & BED format files listing the superset of 3350 manually-annotated microcompartment loops and 363 loop anchors across the five M-to-G1 conditions and for the five captured regions, used in the analysis scripts above. Coordinates are provided for the mm39 reference genome.
For the **BED format anchors file**, loop anchors are listed as plus-and-minus 1kb from the anchor's point coordinate, with the first column of the file being the chromosome for the anchor, the second column being the coordinate of the loop anchor minus 1 kb, the third being the coordinate of the loop anchor plus 1 kb.
For the **BED format loops file**, loop anchors are provided as point coordinates, with the first column of the file being the chromosome for both anchors, the second column being the coordinate of the left loop anchor, and the third column being the coordinate of the right loop anchor.
For the **BEDPE format loops file**, loop anchors are listed as plus-and-minus 1kb from each anchor’s point coordinate. Columns in the BEDPE file are as follows: the first is the chromosome of the left loop anchor, the second is the coordinate of the left loop anchor minus 1 kb, the third is the coordinate of the left loop anchor plus 1 kb, and the remaining three columns are the same for the right loop anchor.

### Calculation of average read counts and/or read-containing bin fractions by contact distance (AverageReadCount&FillFractionByDistPlots.ipynb)
Required packages:
-	cooltools
-	cooler
-	matplotlib

Jupyter notebook used to calculate the fraction of bins in .mcool-derived contact maps which contain at least one read pair at a given resolution and contact distance from the diagonal. The notebook takes an unbalanced .mcool of contacts from RCMC or a comparable method, tabulates the occupied contact bin fraction at specified contact distances, and generates a plot of occupied bin fraction by contact distance.

### Generation of P(s) curves (Ps_CurveAnalysis.ipynb)
Required packages:
-	cooltools
-	cooler
-	matplotlib

Jupyter notebook used to calculate the P(s) curves for provided mcools. The notebook calculates interaction frequencies across contact distances, truncates calculated distances to match the size of RCMC regions, and plots both P(s) curves and derivative plots.

### Finding chromatin features overlapping microcompartment anchors (loopFeatureOverlap.R)
Required packages:
-	plyr
-	dplyr
-	reshape2
-	purrr
-	grid
-	IRanges
-	GenomicRanges
-	arrangements
-	foreach

R script used to classify microcompartment interactions by finding overlap between identified microcompartments (.bedpe) and chromatin features (.bed) such as promoters, enhancers, CTCF binding sites, etc. It outputs individual .bedpe files of interactions according to combinatorial classification of chromatin features (e.g for enhancer (E) and promoter (P): P-P, E-E, E-P, E-null, P-null, null-null), including interactions which have no overlap (null category). Classification can be mutually exclusive (E-P cannot also be P-P) or inclusive.

Example usage:
```
Rscript /path/to/script/loopFeatureOverlap.R -l interactions.bedpe -b promoter.bed,enhancer.bed -i P,E -o outputdirectory/
```

### Plotting loop contact distances and loops per anchor (LoopHistograms&SwarmPlot.ipynb)

Jupyter notebook used to calculate histograms of loop contact distances and the number of loops formed by each loop anchor, as well as swarm plots of the number of loops formed by each anchor separated by anchor identity. The notebook takes a bed file of all loop anchors (see MitosisRCMC_AnchorCalls_PlusMin1kb_mm39.bed) and a bedpe file of all loops (see MitosisRCMC_LoopCalls_PlusMin1kb_mm39.bedpe) and generates histogram and swarm plot visualizations.


### Calculating & plotting the strength of individual interactions (LoopStrengths.ipynb) or plotting aggregate pileup analysis (LoopPileups.ipynb)
Required packages:
-	coolpuppy
-	cooltools
-	cooler

Jupyter notebooks used to calculate strengths of individual microcompartments (LoopStrengths.ipynb) or generate aggregate pileup analysis figures (LoopPileups.ipynb). Each one takes .mcool files of contacts from RCMC, a list of interactions to calculate for (.bedpe format), expected files generated by cooltools for each .mcool, and the captured region, and calculates background corrected observed/expected interaction strengths, either for each interaction individually (output as a .bedpe with additional columns for strengths for each .mcool) or as a pileup (output as a .svg of the aggregate interaction).

### Calculating compartmental identity and plotting saddleplots (CompartmentalizationAnalyses.ipynb)
Required packages:
-	coolpuppy
-	cooltools
-	cooler

Jupyter notebook used to calculate the eigenvectors indicating compartmental identity and plot them as tracks, saddleplots, and strength plots. It takes .mcool files, calculates expected files generated by cooltools for each .mcool, and calculates the eigenvectors within each captured region. These eigenvectors can be tabulated for future reference or used to generate eigenvalue tracks and saddleplots.

### Reproducibility analysis (ReproducibilityAnalysis_HiCRep.ipynb)
Required packages:
-	hicrep
-	cooltools
-	cooler
-	matplotlib

Jupyter notebook used to calculate reproducibility scores for provided mcools. The notebook uses hicrep.py to quantify the similarity between different datasets and plots a heatmap showing the reproducibility scores between all pairs of datasets and hierarchical clustering of datasets by similarity.

## How to cite
This work is shared under an MIT license. If you make use of analysis scripts or data from this work, please cite as follows:
**Goel, V.Y., et al. Dynamics of microcompartment formation at the mitosis-to-G1 transition. *bioRxiv* (2024).**

Please also refer to our polymer simulation code, also available on GitHub at https://github.com/mirnylab/microcompartments.
