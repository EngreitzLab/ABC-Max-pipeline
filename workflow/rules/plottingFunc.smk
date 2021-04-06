
# Added in functionality for data filtered for promoters BUT this is only available for predictions (like ABC) that are have promoters included 
# TODO:  Regarding above: Add rule above to intersect variants (or predictions?) with promoter file and annotate accordingly. Then,
##             Edit AnnotateCredibleSets.R function so that it also calculates enrichments after removing all promoter variants
rule plotTraitEnrichment:
	input: 
		cellTypeEnrichments = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
#		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.eps"),
		outpdf = report(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"), caption="report/CellTypeEnrichment.rst", category="Trait Enrichment Plots", subcategory="{trait}/{pred}")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		codeDir = config["codeDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
	 	isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]), 
	message: "Running enrichment plots"
	run:
		shell(
			"""
			Rscript {params.codeDir}/PlotCellTypeEnrichment.R \
			--outdir {params.outDir} \
			--outPdf {output.outpdf} \
			--outEps {output.outeps} \
			--cellTypes {params.cellTypeTable} \
			--cellTypeEnrichments {input.cellTypeEnrichments} \
			--codeDir {params.codeDir} \
			--trait {wildcards.trait} 
			""")

rule plotTraitEnrichment_noPromoter:
	input:
		cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
#		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.pdf"),
                outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.eps"),
		outpdf = report(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.pdf"), caption="report/CellTypeEnrichment.noPromoter.rst", category="Trait Enrichment Plots", subcategory="{trait}/{pred}")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
                codeDir = config["codeDir"],
                outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
                isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		entry = "enrichment.NoPromoters"
	message: "Running enrichment plots for predictions without promoters"
	run:
		shell(
			"""
			Rscript {params.codeDir}/PlotCellTypeEnrichment.R \
			--outdir {params.outDir} \
			--outPdf {output.outpdf} \
			--outEps {output.outeps} \
			--cellTypes {params.cellTypeTable} \
			--cellTypeEnrichments {input.cellTypeEnrichments_noPromoter} \
			--codeDir {params.codeDir} \
			--trait {wildcards.trait} \
			--entry {params.entry}
			""")

rule plotFractionOverlap:
	input:
		cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
#		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.eps"),
		outpdf = report(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.pdf"), caption="report/CellTypeOverlap.rst", category="Fraction Enhancer Overlap", subcategory="{trait}/{pred}")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
                codeDir = config["codeDir"],
                outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
                isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
                entry = "enrichment"
	message: "Running fraction Overlap plots for predictions"
	run:
		shell(
			"""
                        Rscript {params.codeDir}/PlotFractionOverlap.R \
                        --outdir {params.outDir} \
                        --outPdf {output.outpdf} \
                        --outEps {output.outeps} \
                        --cellTypes {params.cellTypeTable} \
                        --cellTypeEnrichments {input.cellTypeEnrichments_noPromoter} \
                        --codeDir {params.codeDir} \
                        --trait {wildcards.trait} \
                        --entry {params.entry}
                        """)

rule plotFractionOverlap_noPromoter:
	input:
		cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
#		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.eps"),
		outpdf = report(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.pdf"), caption="report/CellTypeOverlap.noPromoter.rst", category="Fraction Enhancer Overlap", subcategory="{trait}/{pred}")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		codeDir = config["codeDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		entry = "enrichment.NoPromoters"
	message: "Running fraction overlap plots for predictions without promoters"
	run:
		shell(
			"""
			Rscript {params.codeDir}/PlotFractionOverlap.R \
			--outdir {params.outDir} \
			--outPdf {output.outpdf} \
                        --outEps {output.outeps} \
                        --cellTypes {params.cellTypeTable} \
                        --cellTypeEnrichments {input.cellTypeEnrichments_noPromoter} \
                        --codeDir {params.codeDir} \
                        --trait {wildcards.trait} \
                        --entry {params.entry}
			""")

rule plotIndividualGenePrecisionRecall:
	input:
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv",outdir=config["outDir"]),
		knownGenes = lambda wildcard: str(config["predDir"]+(trait_config_file.loc[wildcard.trait, "knownGenes"]))
	output:
#		prPdf = os.path.join(config["outDir"], "{pred}/{trait}/GenePrecisionRecall.pdf"),
		prPdf = report(os.path.join(config["outDir"], "{pred}/{trait}/GenePrecisionRecall.pdf"), caption="report/GenePrecisionRecall.rst", category="Precision and Recall", subcategory="{trait}/{pred}")
	params:
		codeDir = config["codeDir"],
		projectDir = config["projectDir"]
	message: "Running precision-recall plot for {wildcards.trait} and {wildcards.pred} predictions"
        run:
                shell(
			"""
			Rscript {params.codeDir}/PlotGenePrecisionRecall.R \
                        --outPdf {output.prPdf} \
                        --genePredTable {input.genePredTable} \
                        --knownGenes {input.knownGenes} \
                        --codeDir {params.codeDir}
                        """)

rule plotGenePrecisionRecall:
	input:
		genePredTable = expand("{outdir}{pred}/{{trait}}/GenePredictions.allCredibleSets.tsv", pred=all_predictions, outdir=config["outDir"]),
		knownGenes = lambda wildcard: str(config["predDir"]+(trait_config_file.loc[wildcard.trait, "knownGenes"]))
	output:
#		prPdf = os.path.join(config["outDir"], "GWAS.{trait}.GenePrecisionRecall.pdf"),
		prPdf = report(os.path.join(config["outDir"], "GWAS.{trait}.GenePrecisionRecall.pdf"), caption="report/GenePrecisionRecall.rst", category="Precision and Recall", subcategory="{trait}")
	params:
		codeDir = config["codeDir"],
		projectDir = config["projectDir"]
	message: "Running precision-recall plot for {wildcards.trait} across all predictions"
	run:
		shell(
			"""
			Rscript {params.codeDir}/PlotGenePrecisionRecall.R \
			--outPdf {output.prPdf} \
			--genePredTable "{input.genePredTable}" \
			--knownGenes {input.knownGenes} \
			--codeDir {params.codeDir}
			""")


rule plotAggregate_cdf: 
	input:
		enrichmentFiles = expand("{outdir}{pred}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config['outDir'], pred=all_predictions)
	output:
#		outfile = os.path.join(config["outDir"], "GWAS.{trait}.cdf.pdf"),
		outDensity = report(os.path.join(config["outDir"], "GWAS.{trait}.density.pdf"), caption="report/GWAS.density.rst", category="Enrichment Density Plots", subcategory="{trait}"),
		outfile = report(os.path.join(config["outDir"], "GWAS.{trait}.cdf.pdf"), caption="report/GWAS.cdf.rst", category="Enrichment Density Plots", subcategory="{trait}")
	params:
		codeDir = config["codeDir"],
		predictors = all_predictions,
		outDir = config["outDir"] 
	message: "Plotting aggregate enrichment CDF and Density across predictions"
	run:
		shell(
			"""
			Rscript {params.codeDir}/PlotEnrichmentAggregate.R \
			--names "{params.predictors}" \
			--tables "{input.enrichmentFiles}" \
			--outPdf {output.outfile} \
			--outDensity {output.outDensity} \
			--outDir {params.outDir}
			""")
	
# TODO: Should we also create aggregate plots for both data that contains promoters and data with filtered out promoters 
rule plotAggregate:
	output:
#		outfile = os.path.join(config["outDir"], "{trait}/{trait}_across_all_predictions.pdf"),
		outfile = report(os.path.join(config["outDir"], "{trait}/{trait}_across_all_predictions.pdf"), caption="report/ComparativeAggregatePlots.rst", category="Comparison Plots")
	params:
		predictorOfChoice = config["predictorOfChoice"],
		predictors = all_predictions,
		codeDir = config["codeDir"],
		outDir = config["outDir"]
	message: "Aggregating enrichment plots across predictors"
	run:
		shell(
			"""
			python {params.codeDir}/plot_aggregate.py \
			--traits {wildcards.trait} \
			--predictor_of_choice {params.predictorOfChoice} \
			--data_outdir {params.outDir} \
			--outdir {params.outDir} \
			--predictors {params.predictors}
			""")
