# ABC-Max: Genome-wide maps of enhancer regulation connect risk variants to disease genes

This repository contains a reproducible workflow for connecting non-coding variants to their target genes through enhancer-gene maps. The workflow has been implemented as a snakemake pipeline and can be run for any set of E-G predictions and variants or credible sets. By default, the pipeline uses the genome-wide Activity-by-Contact (ABC) E-G maps across 131 cell types from Nasser et al. 2020 (https://www.biorxiv.org/content/10.1101/2020.09.01.278093v1.full). This pipeline can be run to reproduce the main results presented in Nasser et al. 2020.

All input files and parameters are specified in the configuration file (ABC-Max.config.json). Multiple sets of predictions and variants can be processed by passing lists to config["predictions"] and config["traits"]. The items on these lists must match JSON dictionaries containing all the required file and directory paths and parameters.

## Required inputs

Enhancer-Gene predictions (predFile) with the following columns
* Required columns
	* chr
	* start (start position of the enhancer element)
	* end (end position of the enhancer element)
	* TargetGene (gene symbol)
	* TargetGeneTSS (TSS of the target gene)
* Optional columns
	* CellType
	* Significance score, e.g. ABC-sore

Variant list (varList)
* Required columns
	* chr
	* position
	* variant (rsID)
* Optional columns
	* A column with scores, e.g. *p*-values

A set of background variants, e.g. all 1000genomes variants, in BED format without column names
* Required columns
	* chr
	* start
	* end
	* variants (rsID)

TODO: check that only "position" is needed, not start and end, like in the variant list.

Chromosome sizes. A tab-delimited file with two columns, chromosome name and size, without column names. A hg19 file compatible with the ABC predictions is provided.

Gene annotations files. Files compatible with the ABC predictions are provided.
* RefSeq gene BED file to pull RefSeq IDs to determine coding/noncoding
* Collapsed RefSeq gene BED file used for E-G predictions

## Optional inputs

Credible set list (csList)
* Required columns

## Dependencies

Dependencies for running the snakemake workflow are distributed as a Conda environment (ABC-Max_env.yaml). 
* snakemake

Additionally, the following R packages are required:
* optparse
* R utils
* data.table
* dplyr
* tidyr
* ggplot2
* gplots
* RColorBrewer

## Running ABC-Max

If starting with the raw-outputs of the ABC-pipeline, the first step in the ABC-Max snakemake (rule preprocessABC) removes promoter elements, filters the predictions to retain those with an ABC score >=0.15, and shrinks the elements from 500bp to 200bp. These steps are executed is the raw ABC prediction (rawPredFile) are provided, but the processed predictions (predFile) are not.

Example command:


###  1. Compute background overlap (rule: computeBackgroundOverlap)
Intersecting a background list of variants with predicted enhancers to compute background rate at which common variants overlap enhancers overall. This is done for all variants and non-coding variants.

### 2. Create variant BED files (rule: createVarFiles)
Creating .bed and .bedgraph files based on the variant list. If a significance score and threshold are provided, creating a list file for the significant variants.

### 3. Compute variant overlap (rule: overlapVariants)
Intersecting (significant) variants with predicted enhancers.

### 4. Annotate variants (rule: annotateVariants)
Calling an R script to run ABC-Max.

Example command: 
```
Rscript AnnotateCredibleSets.R --variants /oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.variant.list.txt --credibleSets /oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.cs.txt --outbase /oak/stanford/groups/akundaje/kmualim/test_code//ABC/IBD/ --trait IBD --codeDir /oak/stanford/groups/akundaje/kmualim/github/ABC-Max-pipeline/Utilities/ --predictionFile /oak/stanford/groups/akundaje/kmualim/GWAS_1/ABC/IBD/IBD.ABC.tsv.gz --bgOverlap /oak/stanford/groups/akundaje/kmualim/GWAS_1/ABC/ABC.OverlapAllSNPs.tsv.gz --cellType TRUE --TargetGene TRUE --TargetGeneTSS TRUE
```
