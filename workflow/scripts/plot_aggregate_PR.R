suppressPackageStartupMessages({
	library(plyr)
	library(dplyr)
	library(tidyr)
	library(stringr)
	library(gtools)
	library(optparse)
	library(DT)
})
# import functions from R script
getPrecisionBaseline <- function(gp) mean(1 / unique(gp[,c("CredibleSet","CredibleSet.nNearbyGenes")])$TotalNearbyGenes)

doOneKnownGeneList <- function(gene.list.name, gp, predictors, maxKnownGenes=1, knownGeneMaxDistance=1000000) {
	## Current logic: Filter the gene prediction table to those credible sets with exactly one known gene nearby
	gp.plot <- gp %>% filter(DistanceToTSS <= knownGeneMaxDistance) %>%
		mutate(knownGene=TargetGene %in% subsetGeneList[[gene.list.name]]) %>%
		group_by(CredibleSet) %>% mutate(nKnownGenes=sum(knownGene)) %>% ungroup() %>%
		filter(nKnownGenes > 0 & nKnownGenes <= maxKnownGenes & CredibleSet.NoncodingWithSigVariant) %>% as.data.frame()
	return(gp.plot)
}

main <- function(){

	option_list <- list(
			    make_option(c("--codeDir"), type="character", default=NA, help="code directory"),
			    make_option(c("--allgenePredTable"), type="character", default=NA, help="Individual GenePredictions file to aggregate"),
			    make_option(c("--outfile"), type="character", default=NA, help= "Aggregated GenePredictions file across all traits"),
			    make_option(c("--PRTable"), type="character", default=NA, help="Precision Recall Table for input into PR Plotting; Compares predictions to knownGenes"),
			    make_option(c("--traitTable"), type="character", default=NA, help="Trait parameters, typically this is the traitTable config file"),
			    make_option(c("--knownGenes"), type="character", default=NA, help="List of known Genes across all traits")
			    )
	opt <- parse_args(OptionParser(option_list=option_list))
	source(paste0(opt$codeDir, "JuicerUtilities.R"))
	source(paste0(opt$codeDir, "CredibleSetTools.R"))

	names = strsplit(opt$allgenePredTable, " ") %>% unlist()
	outfile = (opt$outfile)
	aggPRTable = (opt$PRTable)
	gp.all <- read.delim(names[1], header=T, check.names=F, stringsAsFactors=F, comment.char='#')

	trait_params = read.delim(opt$traitTable, sep='\t', header=TRUE, stringsAsFactors = FALSE)
	traitsToInclude = selectUKBTraits(trait_params)
	knownGenes = strsplit(opt$knownGenes, " ") %>% unlist()
	# only keep files that exists in traitsToInclude
	subsetGeneList = read.delim(knownGenes[1], header=T, check.names=F, stringsAsFactors=F, comment.char='#')
	subsetGeneList$GeneList <- subsetGeneList[, grepl("GeneList" , colnames(subsetGeneList))]
	
	if (length(knownGenes) > 1){
		for (i in 2:length(knownGenes)){
			tryCatch(
			{
				 tempGene = read.delim(file=knownGenes[i], header=T, check.names=F, stringsAsFactors=F, comment.char='#')
				 tempGene$GeneList <- tempGene[, grepl("GeneList" , colnames(tempGene))]
				 subsetGeneList = rbind(subset(subsetGeneList, select=GeneList), subset(tempGene, select=GeneList))
			},
			error=function(cond) {
				message(paste("Gene File does not exist::", i))
			},
			finally={
				message(paste("Gene File Processed:",i))
			})
		}
	}

	if (length(names) > 1){
		for (i in 2:length(names)){
			tryCatch(
			{
				temp = read.delim(file=names[i], header=T, check.names=F, stringsAsFactors=F, comment.char='#')
				gp.all = rbind(gp.all, temp)
			},
			error=function(cond) {
				message(paste("File does not exist::", i))
			},
			finally={
				message(paste("Processed:",i))
			})
		}
	}
	predictors <- colnames(gp.all)[greplany(c("GeneScore.","GenePrediction.","GenePredictionMax."), colnames(gp.all))]
	subsetGeneList <- subsetGeneList %>% distinct(GeneList)
	for (gl in colnames(subsetGeneList)){
		gp.plot <- doOneKnownGeneList(gl, gp.all, predictors, trait_params, subsetGeneList)
		nRecall <- sum(gp.plot$knownGene)
		aggPR <- getPRTable(gp.plot, predictors)
	}
	aggPR$nRecall <- nRecall
	write.table(aggPR, file=aggPRTable, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
	write.table(gp.all, file=outfile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

main()
