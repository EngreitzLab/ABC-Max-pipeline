# coding: utf-8
from os.path import join
import pandas as pd

# This snakefile contains the rules for conducting the ABC-Max analysis as in
# [citation].


pred_config = config["predictionsTable"]
trait_config = config["traitTable"]

preds_config_file = pd.read_table(pred_config).set_index("entry", drop=False)
trait_config_file = pd.read_table(trait_config).set_index("entry", drop=False)

all_predictions = [str(pred) for pred in list(preds_config_file.index)]
all_traits = [str(trait) for trait in list(trait_config_file.index)]

# include other snakemake files 
include: "rules/annotateVariantInputs.smk"
include: "rules/generateTSSInput.smk"
include: "rules/generateVariantOverlapInputs.smk"
include: "rules/plottingFunc.smk"

rule all:
	input:
		expand(os.path.join(config["outDir"], "{pred}/geneTSS.500bp.bed"), pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"), pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.tsv"), pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.tsv"), pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.noPromoter.tsv.gz"), pred=all_predictions),
                expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.noPromoter.tsv"), pred=all_predictions),
                expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.noPromoter.tsv"), pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/all.bg.SNPs.noPromoter.bed.gz"), pred=all_predictions),
		expand("{outdir}{pred}/{trait}/{trait}.bed", outdir=config["outDir"], trait=all_traits, pred=all_predictions),
		expand("{outdir}{pred}/{trait}/{trait}.bedgraph", outdir=config["outDir"], trait=all_traits, pred=all_predictions),
		expand("{outdir}{pred}/{trait}/{trait}.{pred}.tsv.gz", outdir=config["outDir"], trait=all_traits, pred=all_predictions),
		expand("{outdir}{pred}/{trait}/{trait}.{pred}.noPromoter.tsv.gz", outdir=config["outDir"], trait=all_traits, pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/bgVariants.count.tsv"), pred=all_predictions),
                expand(os.path.join(config["outDir"], "{pred}/bgVariants.count.noPromoter.tsv"), pred=all_predictions),
                expand(os.path.join(config["outDir"], "{pred}/bgOverlap.count.tsv"), pred=all_predictions),
                expand(os.path.join(config["outDir"], "{pred}/bgOverlap.count.noPromoter.tsv"), pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt"), trait=all_traits, pred=all_predictions),
		expand("{outdir}{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv", outdir=config["outDir"], trait=all_traits, pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"), trait=all_traits, pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.pdf"), trait=all_traits, pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.pdf"), trait=all_traits, pred=all_predictions),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.pdf"), trait=all_traits, pred=all_predictions),
#		expand("{outdir}{pred}/{trait}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"], trait=all_traits, pred=all_predictions),
#		expand("{outdir}{pred}/{trait}/GenePrecisionRecall.pdf", outdir=config["outDir"], trait=all_traits, pred=all_predictions),
		expand(os.path.join(config["outDir"], "GWAS.{trait}.cdf.pdf"), trait=all_traits), 
		expand(os.path.join(config["outDir"], "GWAS.{trait}.density.pdf"), trait=all_traits),
		expand(os.path.join(config["outDir"], "{trait}/{trait}_across_all_predictions.pdf"), trait=all_traits, pred=all_predictions)


