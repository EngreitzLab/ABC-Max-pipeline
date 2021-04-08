# Benchmarking Pipeline for Enhancer Gene Links 

This repository seeks to serve as one out of three parts of a Comparison Pipeline spearheaded by the ENCODE Consortium. In this module, we seek to compare enhancer gene links with common disease variants attained from two papers: 
1. Fine-mapping inflammatory bowel disease loci to single-variant resolution (Huang et al 2017)
2. Functionally informed fine-mapping and polygenic localization of complex trait heritability (Weissbrod et al 2020)

The analysis done in this paper is also seen in Genome-wide maps of enhancer regulation connect risk variants to disease genes (Nasser et al 2020) that seeks to use an Enhancer-Gene Linking method, Activity-By-Contact (ABC) in mapping fine-mapped GWAS to their genes. 

## Getting Started
1. Activate conda environment via : `conda create env -f abc-max.yml`


Example command:
```
bash log.sh $CODEDIR $OUTDIR
```
