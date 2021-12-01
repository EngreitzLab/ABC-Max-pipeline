
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
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.pdf"),
                outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.noPromoter.eps")
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
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.eps")
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
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeOverlap.{trait}.noPromoter.eps")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		codeDir = config["codeDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		entry = "enrichment.NoPromoters"
	message: "Running fraction overlap plots (excluding promoter regions) for predictions without promoters"
	run:
		shell(
			"""
			set +o pipefail;

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

rule plotFractionOverlapPosteriorProb:
	input:
		allflat = expand("{outdir}{{pred}}/{trait}/data/all.flat.tsv", outdir=config["outDir"], trait=all_traits)
	
	output:
		tsv = os.path.join(config["outDir"], "{pred}/EnrichmentVsPosteriorProb.tsv"),
		noncoding = os.path.join(config["outDir"], "{pred}/EnrichmentVsPosteriorProb.nonCoding.tsv"),
		pdf = os.path.join(config["outDir"], "{pred}/EnrichmentVsPosteriorProb.pdf")
	
	params:
                traitTable = config["traitTable"],
		dataDir = config["traitDir"],
		outDir = os.path.join(config["outDir"], "{pred}/"),
                codeDir = config["codeDir"]
	
	message: "Running fraction overlap plots (including promoter regions) with increasing posterior probability threshold"
	run:
		shell(
			"""
			set +o pipefail;
			
			Rscript {params.codeDir}/PlotVariantFractionPosProb.R \
			--allFlat "{input.allflat}" \
			--datadir {params.dataDir} \
			--traitTable {params.traitTable} \
			--outdir {params.outDir} \
			--codeDir {params.codeDir} \
			--outfile {output.tsv} \
			--pdf {output.pdf} \
			--outfile_noncoding {output.noncoding}	
			""")
		

rule plotIndividualGenePrecisionRecall:
	input:
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv",outdir=config["outDir"]),
		knownGenes = lambda wildcard: str(config["geneListDir"]+(trait_config_file.loc[wildcard.trait, "knownGenes"]))
	output:
		prPdf = os.path.join(config["outDir"], "{pred}/{trait}/GenePrecisionRecall.pdf")
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
		knownGenes = lambda wildcard: str(config["geneListDir"]+(trait_config_file.loc[wildcard.trait, "knownGenes"]))
	output:
		prPdf = os.path.join(config["outDir"], "GWAS.{trait}.GenePrecisionRecall.pdf")
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
		outfile = os.path.join(config["outDir"], "GWAS.{trait}.cdf.pdf"),
		outDensity = os.path.join(config["outDir"], "GWAS.{trait}.density.pdf")
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
	
rule plottingReport: 
	input:
		cellTypeEnrichments = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv"),
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.Dedup.tsv",outdir=config["outDir"]),
		allgenePredTable = expand("{outdir}{pred}/{{trait}}/GenePredictions.allCredibleSets.tsv",outdir=config["outDir"], pred=all_predictions),
		knownGenes = lambda wildcard: str(config["geneListDir"]+(trait_config_file.loc[wildcard.trait, "knownGenes"])),
		enrichmentFiles = expand("{outdir}{pred}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config['outDir'], pred=all_predictions),
		allFlat = expand("{outdir}{{pred}}/{trait}/data/all.flat.tsv", outdir=config["outDir"], trait=all_traits)
	params:
		codeDir = config["codeDir"],
		predictors = all_predictions,
		outDir = config["outDir"],
		predDir = config["predDir"],
		traitDir = config["traitDir"],
		knownGeneMaxDistance = "1000000" ,
		isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		indivOutDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		traitTable = config["traitTable"],
		outEps = "placeholder.eps"
	output: 
		html = os.path.join(config["outDir"], "{pred}/{trait}/GWAS_enrichment_report.html")
	message: "Compiling R Markdown report in html format"
	script: os.path.join(config["codeDir"], "plottingFuncReport.Rmd")


rule getAggregateReport_input:
	input:
		enrichmentFiles = expand("{outdir}{{pred}}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv", outdir=config['outDir'], trait=all_traits)
	params:
		codeDir = config["codeDir"],
		pred = "{pred}", 
		outDir = os.path.join(config["outDir"], "{pred}/")
	output:
		outfile = os.path.join(config["outDir"], "{pred}/{pred}_aggregateTraitEnrichment.tsv"),
		numCellTypes = os.path.join(config["outDir"],"{pred}/{pred}_numEnrichedCellTypesPerTrait.tsv")
	run:
		shell(
			"""
			Rscript {params.codeDir}/plot_all_traits.R \
			--enrichmentTables "{input.enrichmentFiles}" \
			--outDir {params.outDir} \
			--pred "{params.pred}" \
			--outfile {output.outfile} \
			--numCellTypes {output.numCellTypes}
			""")	

rule getAggregatePR_input:
	input:
		allgenePredTable = expand("{outdir}{{pred}}/{trait}/GenePredictions.allCredibleSets.Dedup.tsv",outdir=config["outDir"], trait=all_traits),
		knownGenes = expand("{outdir}GeneLists.{trait}.txt",outdir=config["geneListDir"], trait=all_traits),
	params:
		codeDir = config["codeDir"],
		traitTable = config["traitTable"]
	output:
		outfile = expand("{outdir}{{pred}}/{{pred}}_aggregateGenePredictions.allCredibleSets.Dedup.tsv", outdir=config["outDir"]), 
		aggPRTable = expand("{outdir}{{pred}}/{{pred}}_aggregatePrecisionRecallTable.tsv", outdir=config["outDir"]) 
	run:
		shell(
			"""
			Rscript {params.codeDir}/plot_aggregate_PR.R \
			--codeDir {params.codeDir} \
			--allgenePredTable "{input.allgenePredTable}" \
			--outfile {output.outfile} \
			--PRTable {output.aggPRTable} \
			--traitTable {params.traitTable} \
			--knownGenes "{input.knownGenes}"
			""")

rule plottingAggregateTraitReport:
	input:
		allgenePredTable = expand("{outdir}{{pred}}/{{pred}}_aggregateGenePredictions.allCredibleSets.Dedup.tsv", outdir=config["outDir"]),
		PRtable = expand("{outdir}{{pred}}/{{pred}}_aggregatePrecisionRecallTable.tsv", outdir=config["outDir"]),
		enrichmentFiles = expand("{outdir}{pred}/{pred}_aggregateTraitEnrichment.tsv", outdir=config["outDir"], pred=all_predictions),
		enrichedCellTypes = expand("{outdir}{{pred}}/{{pred}}_numEnrichedCellTypesPerTrait.tsv", outdir=config["outDir"]),
		enrichVsPp_noncoding = expand("{outdir}{{pred}}/EnrichmentVsPosteriorProb.nonCoding.tsv", outdir=config["outDir"])
	params:
		codeDir = config["codeDir"],
		predictors = all_predictions,
		outDir = config["outDir"],
		indivOutDir = os.path.join(config["outDir"], "{pred}/"),
		predDir = config["predDir"],
		traitDir = config["traitDir"],
		knownGeneMaxDistance = "1000000" ,
		isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasCellType"]),
	output:
		html = os.path.join(config["outDir"], "{pred}/GWAS_aggregateTrait_report.html")
	message: "Compiling R Markdown report in html format"
	script: os.path.join(config["codeDir"], "plottingAggregateTraits.Rmd")


rule generateGeneListAggregate:
	input:
		knownGenes = expand("{outdir}/{geneLists}", outdir=config["geneListDir"], geneLists=trait_config_file.loc[:, "knownGenes"])
	output:
		expand("{outdir}/GeneLists.aggregate.txt", outdir=config["resources"])
	message: "Generating Aggregate Gene List File"
	run:
		shell(
			"""
			echo "GeneLists" > {output.aggregateGeneList} \
			cat {input.knownGenes} | grep -v "GeneList*" >> {output.aggregateGeneList}			
			""")


rule plottingAggregateReport:
	input:
		allgenePredTableFiles = expand("{outdir}{pred}/{pred}_aggregateGenePredictions.allCredibleSets.Dedup.tsv", pred=all_predictions, outdir=config["outDir"]),
		enrichmentFiles = expand("{outdir}{pred}/{pred}_aggregateTraitEnrichment.tsv", outdir=config["outDir"], pred=all_predictions),
		enrichedCellTypes = expand("{outdir}{pred}/{pred}_numEnrichedCellTypesPerTrait.tsv", outdir=config["outDir"], pred=all_predictions),
		enrichVsPp_noncoding = expand("{outdir}{pred}/EnrichmentVsPosteriorProb.nonCoding.tsv", outdir=config["outDir"], pred=all_predictions),
		geneLists = expand("{outdir}/GeneLists.aggregate.txt", outdir=config["resources"]) 
	params:
		codeDir = config["codeDir"],
                predictors = all_predictions,
                outDir = config["outDir"],
                predDir = config["predDir"],
                traitDir = config["traitDir"],
                knownGeneMaxDistance = "1000000" 
	output:
		html = os.path.join(config["outDir"], "GWAS_aggregatePredictions_report.html")
	message: "Compiling R Markdown report in html format"
	script: os.path.join(config["codeDir"], "plottingAggregatePredictions.Rmd")


rule plotPropertyReport:
	input:
		enhancerPredictions = expand("{outdir}{pred}_enhancerRegions_signal_coverage.tsv", outdir=config["resources"], pred=all_predictions), 
		candidateRegions = expand("{outdir}{pred}_candidateRegions_signal_coverage.tsv", outdir=config["resources"], pred=all_predictions), 
		mergeEnhancerRegions = expand("{outdir}{pred}/{pred}.mergedEnhancerRegions.tsv.gz", outdir=config["outDir"], pred=all_predictions),
		numGenes = expand("{outdir}{pred}/{pred}.metrics.numGenes.tsv", outdir=config["outDir"], pred=all_predictions), 
		numBiosamplesCounts = expand("{outdir}{pred}/{pred}.metrics.numBiosamplesCounts.tsv", outdir=config["outDir"], pred=all_predictions),
		totalUniquebp = expand("{outdir}{pred}/{pred}.metrics.uniquebp.tsv", outdir=config["outDir"], pred=all_predictions),
		numEGCounts = expand("{outdir}{pred}/{pred}.metrics.numEGCounts.tsv", outdir=config["outDir"], pred=all_predictions),
		EGDistFiles = expand("{outdir}/{pred}_distanceToTSS.tsv", outdir=config["resources"], pred=all_predictions)
	params:
		num_examples = "1000",
		allpredictions = all_predictions
	output:
		html = os.path.join(config["outDir"], "property_aggregate_report.html")
	message: "Getting enhancer properties" 
	script: os.path.join(config["codeDir"], "plottingProperties.Rmd")
	
