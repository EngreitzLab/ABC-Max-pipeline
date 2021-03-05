# coding: utf-8
from os.path import join
import pandas as pd

# This snakefile contains the rules for conducting the ABC-Max analysis as in
# [citation].

#configfile: "ABC-Max.config.json"  ## Specify this on the command line

pred_config = config["predictionsTable"]
trait_config = config["traitTable"]

preds_config_file = pd.read_table(pred_config).set_index("entry", drop=False)
trait_config_file = pd.read_table(trait_config).set_index("entry", drop=False)

# Gathering all the outputs for all sets of predictions and variants
#outputSet = set()

rule all:
	input:
#		outputSet,
		expand(os.path.join(config["outDir"], "{pred}/geneTSS.500bp.bed"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.noPromoter.tsv.gz"), pred=config["predictions"]),
                expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.noPromoter.tsv"), pred=config["predictions"]),
                expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.noPromoter.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/all.bg.SNPs.noPromoter.bed.gz"), pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.bed", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.bedgraph", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.{pred}.tsv.gz", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.{pred}.noPromoter.tsv.gz", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/bgVariants.count.tsv"), pred=config["predictions"]),
                expand(os.path.join(config["outDir"], "{pred}/bgVariants.count.noPromoter.tsv"), pred=config["predictions"]),
                expand(os.path.join(config["outDir"], "{pred}/bgOverlap.count.tsv"), pred=config["predictions"]),
                expand(os.path.join(config["outDir"], "{pred}/bgOverlap.count.noPromoter.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt"), trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"), trait=config["traits"], pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.pdf"), trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/GenePrecisionRecall.pdf", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{trait}/{trait}_across_all_predictions.pdf"), trait=config["traits"], pred=config["predictions"])


rule getGeneTSS:
        input:
                geneList = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred,"genes"]
        output:
                geneTSS = expand("{outdir}{{pred}}/geneTSS.500bp.bed", outdir=config["outDir"])
        params:
                chrSizes = config["chrSizes"]
        run:

                shell(
                        """
                        cat {input.geneList} | perl -lane 'if  ($F[5] == "+" ) {{print $F[0]."\t".$F[1]."\t".$F[1]."\t".$F[3]."\t".$F[4]."\t".$F[5]}} else {{print $F[0]."\t".$F[2]."\t".$F[2]."\t".$F[3]."\t".$F[4]."\t".$F[5]}}' > {output.geneTSS}.tmp
                        bedtools slop -b 250 -i {output.geneTSS}.tmp -g {params.chrSizes} > {output.geneTSS}
                        """)


rule computeBackgroundOverlap:
	input:
		predFile = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred, "predFile"],
		allVariants = config["bgVariants"],
		chrSizes = config["chrSizes"],
		CDS = config["CDS"]
	params:
		cellType = lambda wildcard: str(preds_config_file.loc[wildcard.pred, "cellType"]), 
		outDir = expand("{outdir}{{pred}}", outdir=config["outDir"])
	output:
		overallOverlap = os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"),
		overallOverlapCounts = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.tsv"),
		noncodingOverlap = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.tsv")
	log: os.path.join(config["logDir"], "{pred}.bgoverlap.log")
	message: "Overlapping background variants with predictions: {wildcards.pred}"
	run:
		shell(
			"""
			# TODO: find an alternative to deal with pipefail
			set +o pipefail;
			# Intersecting a background list of variants with predicted enhancers 
			# to compute background rate at which common variants overlap enhancers
			# make output dir
			if [ ! -d {params.outDir} ]
			then
				mkdir {params.outDir}
			fi
			# Compute fraction of variants overlapping predictions in each cell type
			# Finding the relevant columns
			if ({params.cellType}=='True')
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

rule computeBackgroundOverlap_noPromoters:
	input:
		overallOverlap = os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"),
		bgVars = config["bgVariants"],
		CDS = config["CDS"]
	params:
		geneTSS = expand("{outdir}{{pred}}/geneTSS.500bp.bed", outdir=config["outDir"])	
	output:
		overallOverlap_noPromoter = os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.noPromoter.tsv.gz"),
		overallOverlapCounts_noPromoter = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.noPromoter.tsv"),
		noncodingOverlap_noPromoter = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.noPromoter.tsv"),
		bgVars_noPromoter = os.path.join(config["outDir"], "{pred}/all.bg.SNPs.noPromoter.bed.gz")
	run:
		shell(
			"""
			set +o pipefail;

			zcat {input.overallOverlap} | head -1 | gzip -c > {output.overallOverlap_noPromoter}
			zcat {input.overallOverlap} | sed 1d > {input.overallOverlap}.tmp 
			bedtools intersect -v -a {input.overallOverlap}.tmp -b {params.geneTSS} | gzip -c >> {output.overallOverlap_noPromoter} 
			

			# Getting the cell type column and counting
			zcat {output.overallOverlap_noPromoter} | cut -f 7 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.overallOverlapCounts_noPromoter};
			
                        # Compute fraction of noncoding variants overlapping predictions in any cell type
			zcat {output.overallOverlap_noPromoter} | bedtools intersect -v -a stdin -b {input.CDS} | cut -f 1-3,7 | sort | uniq | cut -f 4 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.noncodingOverlap_noPromoter}
			
			rm {input.overallOverlap}.tmp

			# Remove promoter variants from bgVars 
			zcat {input.bgVars} | bedtools intersect -v -a stdin -b {params.geneTSS} | gzip > {output.bgVars_noPromoter}
			""")

rule createVarFiles:
	input:
		varList = lambda wildcard: config["traitDir"]+trait_config_file.loc[wildcard.trait, "varList"]
	output:
		varBed = os.path.join(config["outDir"],"{pred}/{trait}/{trait}.bed"),
		varBedgraph = os.path.join(config["outDir"],"{pred}/{trait}/{trait}.bedgraph"),
		sigvarList = os.path.join(config["outDir"],"{pred}/{trait}/{trait}.sig.varList.tsv"),
	log: os.path.join(config["logDir"], "{trait}.{pred}.createbed.log")
	params:
		varFilterCol = lambda wildcard: trait_config_file.loc[wildcard.trait, "varFilterCol"],
		varFilterThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait, "varFilterThreshold"],
		chrSizes = config["chrSizes"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/")
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
		predFile = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred, "predFile"],
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
			zcat {input.predFile} | sed 1d | bedtools intersect -g {params.chrSizes} -b {input.varBedgraph} -a stdin -wb | gzip >> {output.overlap}
			""")

rule overlapVariants_noPromoter:
	input:
		overlap = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.tsv.gz", outdir=config["outDir"])
	params:
		geneTSS = expand("{outdir}{{pred}}/geneTSS.500bp.bed", outdir=config["outDir"]), 
		chrSizes = config["chrSizes"]
	output:
		overlap_noPromoter = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.noPromoter.tsv.gz", outdir=config["outDir"])
	run:
		shell(
			"""
			set +o pipefail;

			zcat {input.overlap} | head -1 | gzip > {output.overlap_noPromoter}
			zcat {input.overlap} | sed 1d | bedtools intersect -g {params.chrSizes} -b {params.geneTSS} -a stdin | gzip >> {output.overlap_noPromoter}
		 					
			""")

rule generateAnnotateVariantInputs:
	input:
		bgVars = config["bgVariants"],
                bgVars_noPromoter = expand("{outdir}{{pred}}/all.bg.SNPs.noPromoter.bed.gz", outdir=config["outDir"]),
                bgOverlap = expand("{outdir}{{pred}}/{{pred}}.OverlapAllSNPs.tsv.gz", outdir=config["outDir"]),
                bgOverlap_noPromoter = expand("{outdir}{{pred}}/{{pred}}.OverlapAllSNPs.noPromoter.tsv.gz", outdir=config["outDir"])
	params:
		hasCellType = lambda wildcard: str(preds_config_file.loc[wildcard.pred,"hasCellType"])
	output: 
		bgVars_count = expand("{outdir}{{pred}}/bgVariants.count.tsv", outdir=config["outDir"]),
		bgVars_noPromoter_count = expand("{outdir}{{pred}}/bgVariants.count.noPromoter.tsv", outdir=config["outDir"]),
		bgOverlap_count = expand("{outdir}{{pred}}/bgOverlap.count.tsv", outdir=config["outDir"]), 
		bgOverlap_noPromoter_count = expand("{outdir}{{pred}}/bgOverlap.count.noPromoter.tsv", outdir=config["outDir"])
	run:
		shell(
			"""
			set +o pipefail;
			zcat {input.bgVars} | cut -f4 | sort -u | wc -l > {output.bgVars_count}
			zcat {input.bgVars_noPromoter} | cut -f4 | sort -u | wc -l > {output.bgVars_noPromoter_count}
			zcat {input.bgOverlap} | cut -f4,8 | sort -u | awk '{{count[$2]++}}END{{for(j in count) print j"\t"count[j]}}' | sort -u | cut -f2 > {output.bgOverlap_count}
			zcat {input.bgOverlap_noPromoter} | cut -f4,8 | sort -u | awk '{{count[$2]++}}END{{for(j in count) print j"\t"count[j]}}' | sort -u | cut -f2> {output.bgOverlap_noPromoter_count}
			""")
rule annotateVariants:
	input:
		predFile = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred, "predFile"],
		varList = lambda wildcard: config["traitDir"]+trait_config_file.loc[wildcard.trait, "varList"],
		csList = lambda wildcard: config["traitDir"]+trait_config_file.loc[wildcard.trait, "csList"],
		predOverlapFile = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.tsv.gz", outdir=config["outDir"]),
		bgVars = expand("{outdir}{{pred}}/bgVariants.count.tsv", outdir=config["outDir"]),
		bgVars_noPromoter = expand("{outdir}{{pred}}/bgVariants.count.noPromoter.tsv", outdir=config["outDir"]), 
		bgOverlap = expand("{outdir}{{pred}}/bgOverlap.count.tsv", outdir=config["outDir"]), 
		bgOverlap_noPromoter = expand("{outdir}{{pred}}/bgOverlap.count.noPromoter.tsv", outdir=config["outDir"])
	output:
		touch(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt")),
		enrichFile = expand("{outdir}{{pred}}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config["outDir"]),
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"])
	log: os.path.join(config["logDir"], "{trait}.{pred}.annotate.log")
	params:
		cellTypeTable = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred,"celltypeAnnotation"],
		codeDir = config["codeDir"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		predScoreCol = lambda wildcard: preds_config_file.loc[wildcard.pred,"predScoreCol"],
		minPredScore = lambda wildcard: preds_config_file.loc[wildcard.pred,"minPredScore"],
		minPredScorePromoter = lambda wildcard: preds_config_file.loc[wildcard.pred,"minPredScorePromoter"],
		varScoreCol = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterCol"],
		varScoreType = lambda wildcard: trait_config_file.loc[wildcard.trait,"varScoreType"],
		varScoreThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterThreshold"],
		genes = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred,"genes"],
		genesUniq = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred,"genesUniq"],
		geneTSS = expand("{outdir}{{pred}}/geneTSS.500bp.bed", outdir=config["outDir"]),
		hasCellType = lambda wildcard: str(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		hasTargetGeneTSS = lambda wildcard: str(preds_config_file.loc[wildcard.pred,"hasTargetGeneTSS"]),
		chr_sizes = config["chrSizes"]
	message: "Annotating {wildcards.trait} variants with {wildcards.pred} predictions"
	run:
		shell(
                """
                Rscript {params.projectDir}/Utilities/AnnotateCredibleSets.R \
                --variants {input.varList} \
                --credibleSets {input.csList} \
                --predictionFile {input.predOverlapFile} \
                --methodName {wildcards.pred} \
                --outbase {params.outDir} \
                --outEnrichment {output.enrichFile} \
                --outGenePredTable {output.genePredTable} \
                --predScoreCol {params.predScoreCol} \
                --minPredScore {params.minPredScore} \
                --minPredScorePromoters {params.minPredScorePromoter} \
                --backgroundVariants {input.bgVars} \
		--backgroundVariants_noPromoter {input.bgVars_noPromoter} \
                --bgOverlap {input.bgOverlap} \
		--bgOverlap_noPromoter {input.bgOverlap_noPromoter} \
                --trait {wildcards.trait} \
                --codeDir {params.codeDir} \
                --variantScoreCol {params.varScoreCol} \
                --variantScoreThreshold {params.varScoreThreshold} \
                --cellTypeTable {params.cellTypeTable} \
                --genes {params.genes} \
                --genesUniq {params.genesUniq} \
		--geneTSS {params.geneTSS} \
                --hasCellType {params.hasCellType} \
                --hasTargetGeneTSS {params.hasTargetGeneTSS} \
		--chr_sizes {params.chr_sizes}
				""")	

# Added in functionality for data filtered for promoters BUT this is only available for predictions (like ABC) that are have promoters included 
# TODO:  Regarding above: Add rule above to intersect variants (or predictions?) with promoter file and annotate accordingly. Then,
##             Edit AnnotateCredibleSets.R function so that it also calculates enrichments after removing all promoter variants
rule plotTraitEnrichment:
	input: 
		cellTypeEnrichments = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.eps")
	params:
		cellTypeTable = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
	 	isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]), 
#	priority: 1
	message: "Running enrichment plots"
	run:
		shell(
			"""
			Rscript {params.projectDir}Utilities/PlotCellTypeEnrichment.R \
			--outdir {params.outDir} \
			--outPdf {output.outpdf} \
			--outEps {output.outeps} \
			--cellTypes {params.cellTypeTable} \
			--cellTypeEnrichments {input.cellTypeEnrichments} \
			--codeDir {params.projectDir} \
			--trait {wildcards.trait} 
			""")

rule plotTraitEnrichment_noPromoter:
	input:
		cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.pdf"),
                outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.eps")
	params:
		cellTypeTable = lambda wildcard: config["predDir"]+preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
                projectDir = config["projectDir"],
                outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
                isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		entry = "enrichment.NoPromoters"
	message: "Running enrichment plots for predictions without promoters"
	run:
		shell(
			"""
			Rscript {params.projectDir}Utilities/PlotCellTypeEnrichment.R \
			--outdir {params.outDir} \
			--outPdf {output.outpdf} \
			--outEps {output.outeps} \
			--cellTypes {params.cellTypeTable} \
			--cellTypeEnrichments {input.cellTypeEnrichments_noPromoter} \
			--codeDir {params.projectDir} \
			--trait {wildcards.trait} \
			--entry {params.entry}
			""")

rule plotGenePrecisionRecall:
	input:
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"]),
		knownGenes = lambda wildcard: trait_config_file.loc[wildcard.trait, "knownGenes"]
	output:
		prPdf = os.path.join(config["outDir"], "{pred}/{trait}/GenePrecisionRecall.pdf")
	params:
		codeDir = config["codeDir"],
		projectDir = config["projectDir"]
	message: "Running precision-recall plot for {wildcards.trait} and {wildcards.pred} predictions"
	run:
		shell(
			"""
			Rscript {params.projectDir}Utilities/PlotGenePrecisionRecall.R \
			--outPdf {output.prPdf} \
			--genePredTable {input.genePredTable} \
			--knownGenes {input.knownGenes} \
			--codeDir {params.projectDir}
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
			python {params.projectDir}Utilities/plot_aggregate.py \
			--traits {wildcards.trait} \
			--predictor_of_choice {params.predictorOfChoice} \
			--data_outdir {params.outDir} \
			--outdir {params.outDir} \
			--predictors {params.predictors}
			""")
