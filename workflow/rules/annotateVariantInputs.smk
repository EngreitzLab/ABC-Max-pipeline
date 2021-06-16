rule annotateVariants:
	input:
		varList = lambda wildcard: config["traitDir"]+str(trait_config_file.loc[wildcard.trait, "varList"]),
		csList = lambda wildcard: config["traitDir"]+str(trait_config_file.loc[wildcard.trait, "csList"]),
		predOverlapFile = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.tsv.gz", outdir=config["outDir"]),
		bgVars = expand("{outdir}{{pred}}/bgVariants.count.tsv", outdir=config["outDir"]),
		bgVars_noPromoter = expand("{outdir}{{pred}}/bgVariants.count.noPromoter.tsv", outdir=config["outDir"]), 
		bgOverlap = expand("{outdir}{{pred}}/bgOverlap.count.tsv", outdir=config["outDir"]), 
		bgOverlap_noPromoter = expand("{outdir}{{pred}}/bgOverlap.count.noPromoter.tsv", outdir=config["outDir"]),
		geneTSS = expand("{outdir}{{pred}}/geneTSS.500bp.bed", outdir=config["outDir"])
	output:
		touch(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt")),
		touch(expand("{outdir}{{pred}}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config["outDir"])),
		touch(expand("{outdir}{{pred}}/{{trait}}/data/all.flat.tsv", outdir=config["outDir"])),
		touch(expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"])),
		touch(expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.Dedup.tsv", outdir=config["outDir"])),		
		enrichFile = expand("{outdir}{{pred}}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config["outDir"]),
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"]),
		genePredTableDedup = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.Dedup.tsv", outdir=config["outDir"]),
		allflat = expand("{outdir}{{pred}}/{{trait}}/data/all.flat.tsv", outdir=config["outDir"])
	params:
		codeDir = config["codeDir"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		predScoreCol = lambda wildcard: preds_config_file.loc[wildcard.pred,"predScoreCol"],
		minPredScore = lambda wildcard: preds_config_file.loc[wildcard.pred,"minPredScore"],
		minPredScorePromoter = lambda wildcard: preds_config_file.loc[wildcard.pred,"minPredScorePromoter"],
		varScoreCol = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterCol"],
		varScoreType = lambda wildcard: trait_config_file.loc[wildcard.trait,"varScoreType"],
		varScoreThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterThreshold"],
		genes = lambda wildcard: config["dataDir"]+str(preds_config_file.loc[wildcard.pred,"genes"]),
		genesUniq = lambda wildcard: config["dataDir"]+str(preds_config_file.loc[wildcard.pred,"genesUniq"]),
		hasCellType = lambda wildcard: str(preds_config_file.loc[wildcard.pred,"hasCellType"]),
		hasTargetGeneTSS = lambda wildcard: str(preds_config_file.loc[wildcard.pred,"hasTargetGeneTSS"]),
		chr_sizes = config["chrSizes"],
		isEnhancerBed = lambda wildcard: str(preds_config_file.loc[wildcard.pred, 'EnhancerBedFile'])
	log: os.path.join(config["logDir"], "{trait}.{pred}.annotate.log")
	message: "Annotating {wildcards.trait} variants with {wildcards.pred} predictions"
	run:
		shell(
                """
		if [[ $(zcat {input.predOverlapFile} | wc -l) -ge 2 ]]
		then
                	Rscript {params.codeDir}/AnnotateCredibleSets.R \
                	--variants {input.varList} \
                	--credibleSets {input.csList} \
                	--predictionFile {input.predOverlapFile} \
                	--methodName {wildcards.pred} \
                	--outbase {params.outDir} \
                	--outEnrichment {output.enrichFile} \
                	--outGenePredTable {output.genePredTable} \
			--outGenePredTableDedup {output.genePredTableDedup} \
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
                	--genes {params.genes} \
                	--genesUniq {params.genesUniq} \
			--geneTSS {input.geneTSS} \
                	--hasCellType {params.hasCellType} \
                	--hasTargetGeneTSS {params.hasTargetGeneTSS} \
			--chr_sizes {params.chr_sizes} \
			--isEnhancerBed {params.isEnhancerBed}
		fi
				""")	

