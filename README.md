# Benchmarking Pipeline for Enhancer Gene Links 

This repository seeks to serve as one out of three parts of a Comparison Pipeline spearheaded by the ENCODE Consortium. In this module, we seek to compare enhancer gene links with common disease variants attained from two sources: 
1. Fine-mapping inflammatory bowel disease loci to single-variant resolution (Huang et al 2017)
2. Fine-mapping in UK Biobank (Finucane et al): https://www.finucanelab.org/data

## Getting Started
1. Create the conda environment via : `conda env create -f workflow/envs/abc-max.yml`
2. Activate the conda environment : `conda activate abc-max`

## Downloading necessary input files 
* <http://mitra.stanford.edu/kundaje/kmualim/ABC_links/ABC-Max-Resources/ABC-Max-Resources.tar.gz>
* ABC-Max-Resources.tar.gz file contains all the necessary input files to reproduce Fig 1a and Fig 1b in https://www.nature.com/articles/s41586-021-03446-x as well as additional files to run ABC-Max on the fine-mapping done in IBD (Huang et al.) and UK BioBank (Finucane et al)

## Updating the config file 
1. Modify the config file to your file paths, prediction table file and trait table file. 
2. Description for main configuration input json fields are located here: https://docs.google.com/spreadsheets/d/1vf7JvtPHVBrDldiUypQKXl2QLOfKOcATjxjZcV5uZQY/edit?usp=sharing
An example of this config file is located here: ABC-Max-pipeline/config/ABC-Max.example.json
An example of this prediction table file is located here: ABC-Max-pipeline/config/ABC-Max.config-preds.tsv
An example of this traits table file is located here: ABC-Max-pipeline/config/ABC-Max.config-traits.tsv

## Running the pipeline 
To run the pipeline, see `ABC-Max-pipeline/workflow/run_snakemake.sh`
- This is an example to run the pipeline for the Stanford Slurm Cluster 
To run it on your own computer, `snakemake --snakefile Snakefile --configfile ../config/ABC-Max.example.json -j 1 --keep-target-files --rerun-incomplete`

## Example output files
Example output files are located `ABC-Max-pipeline/tests`
