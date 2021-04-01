
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
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
	 	isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]), 
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
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
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

rule plotFractionOverlap:
	input:
		cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.eps")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
                projectDir = config["projectDir"],
                outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
                isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
                entry = "enrichment"
	message: "Running fraction Overlap plots for predictions"
	run:
		shell(
			"""
                        Rscript {params.projectDir}Utilities/PlotFractionOverlap.R \
                        --outdir {params.outDir} \
                        --outPdf {output.outpdf} \
                        --outEps {output.outeps} \
                        --cellTypes {params.cellTypeTable} \
                        --cellTypeEnrichments {input.cellTypeEnrichments_noPromoter} \
                        --codeDir {params.projectDir} \
                        --trait {wildcards.trait} \
                        --entry {params.entry}
                        """)

rule plotFractionOverlap_noPromoter:
	input:
		cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.eps")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		entry = "enrichment.NoPromoters"
	message: "Running fraction overlap plots for predictions without promoters"
	run:
		shell(
			"""
			Rscript {params.projectDir}Utilities/PlotFractionOverlap.R \
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
		knownGenes = lambda wildcard: str(config["predDir"]+(trait_config_file.loc[wildcard.trait, "knownGenes"]))
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


rule plotAggregate_cdf: 
	input:
		enrichmentFiles = expand("{outdir}{pred}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config['outDir'], pred=all_predictions)
	output:
		outfile = os.path.join(config["outDir"], "GWAS.{trait}.cdf.pdf"),
		outDensity = os.path.join(config["outDir"], "GWAS.{trait}.density.pdf")
	params:
		projectDir = config["projectDir"],
		predictors = all_predictions,
		outDir = config["outDir"] 
	message: "Plotting aggregate enrichment CDF and Density across predictions"
	run:
		shell(
			"""
			Rscript {params.projectDir}Utilities/PlotEnrichmentAggregate.R \
			--names "{params.predictors}" \
			--tables "{input.enrichmentFiles}" \
			--outPdf {output.outfile} \
			--outDensity {output.outDensity} \
			--outDir {params.outDir}
			""")
	
# TODO: Should we also create aggregate plots for both data that contains promoters and data with filtered out promoters 
rule plotAggregate:
	output:
		outfile = os.path.join(config["outDir"], "{trait}/{trait}_across_all_predictions.pdf")
	params:
		predictorOfChoice = config["predictorOfChoice"],
		predictors = all_predictions,
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
