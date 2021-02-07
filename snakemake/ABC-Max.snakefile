# coding: utf-8
from os.path import join
import pandas as pd

# This snakefile contains the rules for conducting the ABC-Max analysis as in
# [citation].

#configfile: "ABC-Max.config.json"  ## Specify this on the command line

pred_config = "ABC-Max.config-preds.tsv"
trait_config = "ABC-Max.config-traits.tsv"

preds_config_file = pd.read_table(pred_config).set_index("entry", drop=False)
trait_config_file = pd.read_table(trait_config).set_index("entry", drop=False)
trait_config_file.head()

# Gathering all the outputs for all sets of predictions and variants
outputSet = set()

rule all:
	input:
#		outputSet,
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.tsv"), pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.bed", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.bedgraph", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.{pred}.tsv.gz", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"])
		expand(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt"), trait=config["traits"], pred=config["predictions"])
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"), trait=config["traits"], pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{trait}/{trait}_across_all_predictions.pdf"), trait=config["traits"], pred=config["predictions"])


rule computeBackgroundOverlap:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"],
		allVariants = config["bgVariants"],
		chrSizes = config["chrSizes"],
		CDS = config["CDS"]
	params:
		cellType = lambda wildcard: preds_config_file.loc[wildcard.pred, "cellType"]
	output:
		overallOverlap = os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"),
		overallOverlapCounts = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.tsv"),
		noncodingOverlap = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.tsv"),
		outDir = directory(expand("{outdir}{{pred}}", outdir=config["outDir"]))
	log: os.path.join(config["logDir"], "{pred}.bgoverlap.log")
	message: "Overlapping background variants with predictions: {wildcards.pred}"
	run:
		shell(
			"""
			# TODO: find an alternative to deal with pipefail
			set +o pipefail;
			# Intersecting a background list of variants with predicted enhancers 
			# to compute background rate at which common variants overlap enhancers
			# overall

			# TODO: if cell type column is not provided, use all predictions, or require Celltype col?
		        # make output dir
			if [ ! -d {output.outDir} ]
			then
				mkdir {output.outDir}
			fi	
			# Compute fraction of variants overlapping predictions in each cell type
			# Finding the relevant columns
			if {params.cellType}
			then
				zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType | sed 1d | sort -k 1,1 -k 2,2n | uniq | bedtools sort -i stdin -faidx {input.chrSizes} | \
				bedtools intersect -sorted -g {input.chrSizes} -a {input.allVariants} -b stdin -wa -wb | gzip > {output.overallOverlap};
			else
				zcat {input.predFile} | csvtk cut -t -f chr,start,end | sed 1d | sort -k 1,1 -k 2,2n | uniq | bedtools sort -i stdin -faidx {input.chrSizes} | \
				bedtools intersect -sorted -g {input.chrSizes} -a {input.allVariants} -b stdin -wa -wb | gzip > {output.overallOverlap};
			fi
			
			# Getting the cell type column and counting
			 zcat {output.overallOverlap} | cut -f 7 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.overallOverlapCounts};

			# Compute fraction of noncoding variants overlapping predictions in any cell type
			 zcat {output.overallOverlap} | bedtools intersect -v -a stdin -b {input.CDS} | cut -f 1-3,7 | sort | uniq | cut -f 4 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.noncodingOverlap}
			
			""")

rule createVarFiles:
	input:
		varList = lambda wildcard: trait_config_file.loc[wildcard.trait, "varList"]
	output:
		varBed = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.bed", outdir=config["outDir"]),
		varBedgraph = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.bedgraph", outdir=config["outDir"]),
		sigvarList = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.sig.varList.tsv", outdir=config["outDir"])
	log: os.path.join(config["logDir"], "{trait}.{pred}.createbed.log")
	params:
		varFilterCol = lambda wildcard: trait_config_file.loc[wildcard.trait, "varFilterCol"],
		varFilterThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait, "varFilterThreshold"],
		outDir = directory(expand("{outdir}{{pred}}/{{trait}}/", outdir=config["outDir"])),
		chrSizes = config["chrSizes"]
	message: "Creating variant BED files"
	run:
		if {params.varFilterCol} is not None:
			shell(
				"""
				# make output dir 
				if [ ! -d {params.outDir} ]
                       		then
                                	mkdir {params.outDir}
                        	fi
				# Subsetting the variant list based on significance
				# Finding the score colum
				#scoreCol=$(awk -v RS='\\t' '/{params.varFilterCol}/{{print NR; exit}}' {input.varList});

				# Filtering to retain variants exceeding the threshold
				#awk '{{ if ($($scoreCol) >= {params.varFilterThreshold}) {{ print }} }}' {input.varList}  > {output.sigvarList};
				cat {input.varList} | csvtk -t filter -f "{params.varFilterCol}>={params.varFilterThreshold}" > {output.sigvarList};
	
				# Creating the bed file
				# Finding and cutting chr, position, and variant columns
				# TODO: do not require start and stop, only position?
				cat {output.sigvarList} | csvtk cut -t -f chr,position,variant | sed '1d' | awk -F "\\t" "\$1 = \$1 FS \$2-1 FS \$2 FS \$3 FS" | cut -f1-4 | sed -e 's/8.1e+07/81000000/g'> {output.varBed};

				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cat {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				""")
		else:
			shell(
				"""
				#fi
				# Creating the bed file for all variants
				# Finding and cutting chr, pos, and var columns
				cat {input.varList} | csvtk cut -t -f chr,position,variant | sed '1d' | awk -F "\\t" "\$1 = \$1 FS \$2-1 FS \$2 FS \$3 FS" | cut -f1-4 > {output.varBed};

				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cat {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				""")
	

rule overlapVariants:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"],
		varList = lambda wildcard: trait_config_file.loc[wildcard.trait, "varList"],
		varBed = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.bed", outdir=config["outDir"]),
		varBedgraph = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.bedgraph", outdir=config["outDir"])
	output:
		overlap = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.tsv.gz", outdir=config["outDir"])
	log: os.path.join(config["logDir"], "{trait}.{pred}.overlap.log")
	params:
		chrSizes = config["chrSizes"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers"
	run:
		shell(
			"""
			# TODO: find an alternative to deal with pipefail
			set +o pipefail;

			# Creating an empty file with the final columns
			zcat {input.predFile} | head -1 | awk '{{ print $0 "\\tvariant.chr\\tvariant.start\\tvariant.end\\tQueryRegionName" }}' | gzip > {output.overlap};
	
			# Intersecting variants with predictions
			zcat {input.predFile} | sed 1d | bedtools intersect -sorted -g {params.chrSizes} -b {input.varBedgraph} -a stdin -wb | gzip >> {output.overlap}
			""")


rule annotateVariants:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"],
		varList = lambda wildcard: trait_config_file.loc[wildcard.trait, "varList"],
		csList = lambda wildcard: trait_config_file.loc[wildcard.trait, "csList"],
		predOverlapFile = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.tsv.gz", outdir=config["outDir"]),
		bgVars = config["bgVariants"],
		bgOverlap = expand("{outdir}{{pred}}/{{pred}}.OverlapAllSNPs.tsv.gz", outdir=config["outDir"]) 
	output:
		touch(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt")),
		enrichFile = expand("{outdir}{{pred}}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config["outDir"])
	log: os.path.join(config["logDir"], "{trait}.{pred}.annotate.log")
	params:
                cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred,"celltypeAnnotation"],
                codeDir = config["codeDir"],
                projectDir = config["projectDir"],
                outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
                scoreCol = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterCol"],
                scoreType = lambda wildcard: trait_config_file.loc[wildcard.trait,"varScoreType"],
                scoreThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterThreshold"],
                ctrlThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait,"varCtrlThreshold"],
                gex = config["gex"],
                promoterActivityRef = config["promoterActivityRef"],
                cellTypeCov = config["cellTypeCov"],
                specificityBackground = config["specificityBackground"],
                housekeepingList = config["housekeepingList"],
                predColMap = config["predColMap"],
                geneLists = lambda wildcard:preds_config_file.loc[wildcard.pred,"genes"],
                genesUniq = lambda wildcard: preds_config_file.loc[wildcard.pred,"genesUniq"],
                cellType = lambda wildcard: preds_config_file.loc[wildcard.pred,"cellType"],
                TargetGene = lambda wildcard: preds_config_file.loc[wildcard.pred,"TargetGene"],
                isTargetGene = lambda wildcard: preds_config_file.loc[wildcard.pred,"TargetGeneTSS"]
	message: "Annotating {wildcards.trait} variants with {wildcards.pred} predictions"
	run:
		
		#if using ABC predictions, plotting some additional features
                if any(s.startswith('ABC') for s in list({wildcards.pred})):
                        shell(
                                """
                                Rscript {params.projectDir}/AnnotateCredibleSets.R \
                                --variants {input.varList} \
                                --predictionFile {input.predOverlapFile} \
                                --backgroundVariants {input.bgVars} \
                                --bgOverlap {input.bgOverlap} \
                                --outbase {params.outDir} \
                                --trait {wildcards.trait} \
                                --credibleSets {input.csList} \
                                --codeDir {params.codeDir} \
                                --cellTypeTable {params.cellTypeTable} \
                                --gex {params.gex} \
                                --promoterActivityRef {params.promoterActivityRef} \
                                --cellTypeCov {params.cellTypeCov} \
                                --specificityBackground {params.specificityBackground} \
                                --housekeepingList {params.housekeepingList} \
                                --predColMap {params.predColMap} \
                                --genes {params.geneLists} \
                                --genesUniq {params.genesUniq} \
                                --cellType {params.cellType} \
                                --TargetGene {params.TargetGene} \
                                --TargetGeneTSS {params.isTargetGene}
                                """)
                else:
                        shell(
                                """
                                Rscript {params.projectDir}/AnnotateCredibleSets.R \
                                --variants {input.varList} \
                                --predictionFile {input.predOverlapFile} \
                                --backgroundVariants {input.bgVars} \
                                --bgOverlap {input.bgOverlap} \
                                --isABC FALSE \
                                --outbase {params.outDir} \
                                --trait {wildcards.trait} \
                                --credibleSets {input.csList} \
                                --codeDir {params.codeDir} \
                                --variantScoreCol {params.scoreCol} \
                                --scoreType {params.scoreType} \
                                --variantScoreThreshold {params.scoreThreshold} \
                                --variantCtrlScoreThreshold {params.ctrlThreshold} \
                                --cellTypeTable {params.cellTypeTable} \
                                --gex {params.gex} \
                                --promoterActivityRef {params.promoterActivityRef} \
                                --cellTypeCov {params.cellTypeCov} \
                                --specificityBackground {params.specificityBackground} \
                                --genes {params.geneLists} \
                                --genesUniq {params.genesUniq} \
                                --cellType {params.cellType} \
                                --TargetGene {params.TargetGene} \
                                --TargetGeneTSS {params.isTargetGene}
                                --housekeepingList {params.housekeepingList} \
                                --predColMap {params.predColMap} 
				""")	
# Added in functionality for data filtered for promoters BUT this is only available for predictions (like ABC) that are have promoters included 
rule runTraitEnrichment:
	input: 
		cellTypeEnrichments = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
		outfile = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.noPromoter.{trait}.tsv"),
	 	isCellType = lambda wildcard: preds_config_file.loc[wildcard.pred,"cellType"], 
		hasPromoterColumn = lambda wildcard: preds_config_file.loc[wildcard.pred,"hasPromoter"]	
	message: "Running enrichment plots"
	run:
		if {params.isCellType}=={"TRUE"} and {params.hasPromoterColumn}=={"TRUE"}:
			shell(
				"""
				Rscript {params.projectDir}PlotCellTypeEnrichment.R \
				--outdir {params.outDir} \
				--cellTypes {params.cellTypeTable} \
				--cellTypeEnrichments {input.cellTypeEnrichments} \
				--codeDir {params.projectDir} \
				--trait {wildcards.trait} 
				
                               	Rscript {params.projectDir}PlotCellTypeEnrichment.R \
                               	--outdir {params.outDir} \
                               	--cellTypes {params.cellTypeTable} \
                               	--cellTypeEnrichments {params.cellTypeEnrichments_noPromoter} \
                               	--codeDir {params.projectDir} \
                               	--trait {wildcards.trait} \
				--noPromoter TRUE 
                               	""")
		elif {params.isCellType}=={"TRUE"}:
			shell(
				"""
				Rscript {params.projectDir}PlotCellTypeEnrichment.R \
				--outdir {params.outDir} \
                                --cellTypes {params.cellTypeTable} \
                                --cellTypeEnrichments {input.cellTypeEnrichments} \
                                --codeDir {params.projectDir} \
                                --trait {wildcards.trait} \
				""")
		else:
			shell(
				"""
				touch {output.outfile}
				""")

# TODO: Should we also create aggregate plots for both data that contains promoters and data with filtered out promoters 
rule plotAggregate:
	output:
		outfile = os.path.join(config["outDir"], "{trait}/{trait}_across_all_predictions.pdf")
	params:
		predictorOfChoice = config["predictorOfChoice"],
		predictors = config["predictions"],
		projectDir = config["projectDir"],
		outDir = config["outDir"]
		
	message: "Aggregating enrichment plots across predictors"
	run:
		shell(
			"""
			python {params.projectDir}plot_aggregate.py \
			--traits {wildcards.trait} \
			--predictor_of_choice {params.predictorOfChoice} \
			--data_outdir {params.outDir} \
			--outdir {params.outDir} \
			--predictors {params.predictors}
			""")
