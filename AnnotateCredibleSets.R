##################################################################
## Jesse Engreitz
## May 14, 2017
## Annotate GWAS credible sets with enhancer-gene predictions
## Added to codebase Dec 9, 2018

##  This executable script takes in a set of variant predictions (variants overlapped with ABC),
##  plus information about a set (or subset) of these variants and their credible sets,
##  and outputs a series of useful data tables that can be processed or analyzed in various
##  ways by other scripts.

# TODO: refactor so that finemapping info is optional

suppressPackageStartupMessages(library("optparse"))

# Required columns for predictions and overlaps:
# match ABC.Score/Score.Fraction

# Parse options from commans line
option.list <- list(
  make_option("--variants", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.variant.list.txt", help="File containing variants to consider"),
  make_option("--credibleSets", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.cs.txt", help="File containing credible set annotations"),
  make_option("--predictionFile", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/abc.tsv.gz", help="Prediction overlap file (.txt.gz) output from intersecting variants with predicted enhancers"),
  make_option("--isABC", type="logical", default=FALSE, help="Using ABC predictions? If TRUE, some additional metrics are plotted."),
  make_option("--outbase", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_out/IBD/", help="Output file basename"),
  # TODO: edit downstream steps so that these are only used for ABC predictions
  make_option("--cutoff", type="numeric", default=0.015, help="Cutoff on ABC score for distal elements"),
  make_option("--cutoffTss", type="numeric", default=0.1, help="Cutoff on ABC score for tss/promoter elements"),
  make_option("--trait", type="character", default="IBD", help="Name of the trait or disease"),
  make_option("--variantScoreCol", type="character", default="PosteriorProb", help="Score determining variant significance. If set to NULL, all variants in the input file will be included"),
  make_option("--scoreType", type="character", default="PP", help="Type of variantScoreCol, used in plot labels."),
  make_option("--variantScoreThreshold", type="numeric", default=0.1, help="Score cutoff for desired variants to analyze, e.g. PP>=0.1"),
  make_option("--variantCtrlScoreThreshold", type="numeric", default=0.01, help="Score cutoff for for a control set of variants, e.g. PP<0.01"),
  make_option("--backgroundVariants", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/all.bg.SNPs.bed.gz", help="A set of background variants to use in the enrichment analysis. Bed format with chr, start, end, rsID"),
  make_option("--genes", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/RefSeqCurated.170308.bed", help="RefSeq gene BED file; this is to pull RefSeq IDs to determine coding/noncoding"),
  make_option("--genesUniq", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed", help="Collapsed RefSeq gene BED file used for E-G predictions"),
  make_option("--cellTypeTable", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/ABC-GWAS/data/CellTypes.Annotated.ABCPaper.txt", help="Table with annotations of cell types, with columns 'CellType', 'Categorical.*', 'Binary.*' for plotting enrichments"),
  make_option("--relevantCellTypes", type="character", default="Binary.IBDRelevant", help="Column in cell type table containing mask for cell types that are 'relevant' to the trait"),
  make_option("--tissueCategory", type="character", default="Categorical.IBDTissueAnnotations", help="Column in the cell type table containing tissue type categories"),
  # Binary?
  make_option("--predColForStats", type="character", default="ConnectionStrengthRank.Binary.IBDRelevant", help="Prediction to use for making stats"),
  make_option("--gex", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/GeneTSSActivityQuantile.tsv", help="Filename with table of gene expression (or promoter activity) quantiles"),
  make_option("--gexQuantileCutoff", type="numeric", default=0.4, help="Gene expression quantile cutoff to count gene as expressed"),
  make_option("--promoterActivityRef", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/GeneList.txt", help="File from 1 cell type to use to extract distribution of promoter activity across genes (used for TSS-activity weighted predictions)"),
  make_option("--cellTypeCov", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/covariance.all.tsv"),
  make_option("--specificityBackground", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/allSpecificityScores.SparseMatrix.rds"),
  make_option("--removePromoterVariants", type="logical", default=FALSE, help="Remove credible sets with promoter variants from the filter.cs list"),
  make_option("--removeCodingVariants", type="logical", default=TRUE, help="Remove credible sets with coding variants from the filter.cs list"),
  make_option("--removeSpliceSiteVariants", type="logical", default=TRUE, help="Remove credible sets with splice site variants from the filter.cs list"),
  make_option("--housekeepingList", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/Human.HKGeneList.txt", help="List of housekeeping / ubiquitously expressed genes; will ignore ABC connections to these genes"),
  make_option("--predColMap", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Test_data/ABCColMap.txt", help="Needed for ABC predictions from new codebase 191221; pass 'NULL' to ignore"),
  make_option("--codeDir", type="character", default="/home/groups/engreitz/Projects/ABC-Max-pipeline/Utilities/", help="Directory to code base")
)

opt <- parse_args(OptionParser(option_list=option.list))

#if (DEBUG <- FALSE) {
#  opt$codeDir <- "/home/groups/engreitz/Projects/ABC-Max-pipeline/LanderLab-EP-Prediction/"
#  opt$variants <- "data/variant.list.IBD.FM.txt"
#  opt$predictionFile <- "analysis/IBD/abc.tsv.gz"
#  opt$outbase <- "analysis/IBD/"
#  opt$credibleSets <- "data/all.cs.IBD.FM.txt"
#  opt$cellTypeTable <- "data/CellTypes.Annotated.txt"
#  opt$geneLists <- "data/IBD.GeneLists.181229.txt"
#}

# Loading libraries and utilities
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "JuicerUtilities.R")))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "CredibleSetTools.R")))

# Function to save progress
saveProgress <- function() save.image(file=paste0(opt$outbase,"AnnotateCredibleSets.RData"))
saveProgress()

##############################################################################
## Load common data

# All genes
genes <- readBed(opt$genes)
genes$symbol <- unlist(lapply(strsplit(as.character(as.matrix(genes$name)), ";"), "[", 1))
genes <- addTSSToBED(genes)

# Genes used in the predictions
genes.uniq <- readBed(opt$genesUniq)
genes.uniq <- addTSSToBED(genes.uniq)

# Housekeeping genes to ignore
hk.list <- read.delim(opt$housekeepingList, header=F, stringsAsFactors=F)[,1]

#  Vector of all diseases/traits in the variant list
# TODO: assuming one trait per variant list, DONE
variant.list <- read.delim(opt$variants, check.names=F)
#diseases <- unique(variant.list$Disease)

# Finding a vector of relevant cell types
# TODO: Requirements for this file?
cell.type.annot <- read.delim(opt$cellTypeTable, check.names=F)
cell.type.list <- cell.type.annot$CellType

# All credible sets
all.cs <- read.delim(opt$credibleSets, check.names=F, stringsAsFactors=F)
all.cs$CredibleSet <- factor(all.cs$CredibleSet)
# Optionally removing all coding, splice site, and promoter variants. These sets
# are used in the enrichment analysis.
filter.cs <- subset(all.cs, (!AnyCoding | !opt$removeCodingVariants) &
                      (!AnySpliceSite | !opt$removeSpliceSiteVariants) &
                      (!AnyPromoter | !opt$removePromoterVariants))  ## Consider whether to factor this out of this script

variant.list.filter <- subset(variant.list, CredibleSet %in% filter.cs$CredibleSet)

# Finding CSs and variants that exceed the provided threhold. If no variant score
# and threshold are provided, using all variants in downstream steps
# eval(parse(text = opt$variantScoreCol))
sigScore.cs = NULL
variant.list.sigScore = NULL
if (!(is.null(opt$variantScoreCol)) & !(is.null(opt$variantScoreThreshold))) {
  sigScore.cs <- subset(filter.cs, CredibleSet %in% subset(variant.list.filter, eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold)$CredibleSet)
  variant.list.sigScore <- subset(variant.list.filter, CredibleSet %in% sigScore.cs$CredibleSet) ## this is the list of variants in credible sets with at least 1 variant with posterior prob >10%
}

# Overlapping variants with predictions
# TODO: edit loadVariantOverlap() so that predColMap is no longer needed
# TODO: generalize so that this can be run for any set of predictions with
# particular column names
if (opt$isABC){
  predColMap <- if (!is.null(opt$predColMap)) read.delim(opt$predColMap, stringsAsFactors=F) else NULL
  overlap <- loadVariantOverlap(opt$predictionFile, genes.uniq, genes, variant.names=variant.list$variant, colMap=predColMap)
  overlap <- filterVariantOverlap(overlap, opt$cutoff, opt$cutoffTss, hk.list)
  overlap <- addPromoterWeightedPrediction(overlap, opt$promoterActivityRef)
} else {
  overlap <- loadVariantOverlap(opt$predictionFile, genes.uniq, genes, variant.names=variant.list$variant)
}

# Annotating overlaps
all.flat <- annotateVariantOverlaps(overlap, variant.list, all.cs)
all.flat <- subset(all.flat, CellType %in% cell.type.list)  ## IMPORTANT CHANGE
filter.flat <- subset(all.flat, CredibleSet %in% filter.cs$CredibleSet)

sigScore.flat = NULL
if (!(is.null(opt$variantScoreCol)) & !(is.null(opt$variantScoreThreshold))) {
  sigScore.flat <- subset(filter.flat, CredibleSet %in% sigScore.cs$CredibleSet)
}

variant.by.cells <- getVariantByCellsTable(filter.flat)
variant.by.genes <- getVariantByGenesTable(filter.flat)
genes.by.cells <- getGenesByCellsTable(filter.flat)

# Finding all diseases/traits in the credible set table
# Assume that the cs has only one trait
#traits <- unique(all.cs$Disease)
trait <- opt$trait

# Not using the gene lists?
#gene.lists <- loadGeneLists(opt$geneLists)

#if (!is.null(opt$e2gAltMethods)) {
#  alt.overlap <- loadVariantOverlap(opt$e2gAltMethods, genes.uniq, genes, variant.names=variant.list$variant, overwriteTSS=TRUE)
#  alt.overlap <- annotateVariantOverlaps(alt.overlap, variant.list, all.cs)
#  gene.lists <- c(gene.lists, getGeneListsFromE2GOverlap(alt.overlap, opt$posteriorProb))
#} else {
#  alt.overlap <- NULL
#}

# Redundant?
#if (!is.null(opt$geneScores)) {
#  gene.scores <- read.delim(opt$geneScores, check.names=F, stringsAsFactors=F)
#  gene.lists <- c(gene.lists, getGeneListsFromGeneScores(alt.overlap, opt$posteriorProb))
#}

# Creating an output directory and writing the result to files
dir.create(paste0(opt$outbase,"data/"))
write.tab(all.flat, file=paste0(opt$outbase, "data/all.flat.tsv"))
write.tab(filter.flat, file=paste0(opt$outbase, "data/filter.flat.tsv"))
if (!(is.null(sigScore.flat))) {
  write.tab(sigScore.flat, file=paste0(opt$outbase, "data/sigScore.flat.tsv"))
}

# Reading in the gene TSS activity quantile information
gex <- read.delim(opt$gex, check.names=F)

#pp.label <- paste0("pp",opt$posteriorProb*100)

############################################################################
## Basic stats about credible sets
# TODO: if no score/threshold is provided, use all variants.
# TODO: fix posteriorProb here, DONE
stats <- list()
getBasicStats <- function(all.cs, filter.cs, sigScore.cs, variant.list.sigScore, sigScore.flat) {
  res <- list()
  res$`Number of credible sets (all.cs)` <- nrow(all.cs)
  res$`Number of all.cs sets that overlap any enhancer` <- with(all.flat, length(unique(CredibleSet)))
  res$`Number of credible sets with no coding or splice variants (filter.cs)` <- nrow(filter.cs)
  res$`Number of filter.cs sets that overlap any enhancer` <- with(filter.flat, length(unique(CredibleSet)))
  #if (!is.null(opt$ctrlProb)) res$`Num Variants with PP below ctrl threshold in sigScore.cs` <- with(variant.list.sigScore, sum(PosteriorProb < opt$ctrlProb))
  
  if (!is.null(opt$variantScoreCol)) {
    res$`Number of credible sets with no coding or splice variants, and at least one variant with score above threshold (sigScore.cs)` <- nrow(sigScore.cs)
    res$`Number of sigScore.cs sets that have any variant (regardless of significance score) that overlaps any enhancer` <- with(subset(all.flat, CredibleSet %in% sigScore.cs$CredibleSet), length(unique(CredibleSet)))
    res$`Number of variants in sigScore.cs sets` <- nrow(variant.list.sigScore)
    res$`Num Variants with score above threshold in sigScore.cs` <- with(variant.list.sigScore, sum(eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold))
    res$`Num variants above score threshold in sigScore.cs that overlap any enhancer` <- with(sigScore.flat, length(unique(QueryRegionName[eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold])))
    res$`Num sigScore credible sets with variants above score threshold that overlap any enhancer` <- with(sigScore.flat, length(unique(CredibleSet[eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold])))
    for (bin in cell.bins) {
      res[[paste0("Num variants above the score threshold in sigScore.cs that overlap an enhancer in any ", bin, " cell type")]] <- 
        with(sigScore.flat, length(unique(QueryRegionName[eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold & CellType %in% with(cell.type.annot, CellType[get(bin)])])))
      res[[paste0("Num credible sets with variants above the score threshold that overlap an enhancer in any ", bin, " cell type")]] <- 
        with(sigScore.flat, length(unique(CredibleSet[eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold & CellType %in% with(cell.type.annot, CellType[get(bin)])])))
    }
  }
  
  
  return(res)
}

stats[['all']] <- getBasicStats(all.cs, filter.cs, sigScore.cs, variant.list.sigScore, sigScore.flat)


saveProgress()


##############################################################################
## Calculate enrichment per cell type by comparing the significant/provided
## variants with a set of background variants

# A set of background variants
bgVars <- read.delim(gzfile(opt$backgroundVariants), check.names=F, header=F)

# Use the set of background variants instead of ctrlPP
# Remove the column from output that uses ctrlPP
#if (!is.null(opt$ctrlProb)) {
edir <- paste0(opt$outbase, "enrichment/")
dir.create(edir)

pdf(file=paste0(opt$outbase, "Enrichment.CellType.vsScore.pdf"), width=5, height=5)
# Note: assuming one trait
#for (trait in unique(variant.list.filter$Disease)) {
#  tryCatch({
#curr.vl <- subset(variant.list.filter, Disease==trait)
#variant.by.cells <- getVariantByCellsTable(filter.flat)
# TODO: use bgVars instead of ctrlPP, DONE
# TODO: change "variant" column to rsID?
# TODO: if the prediction file includes a column "Promoter" that has both TRUE 
# and FALSE values, calculate enrichment with and without promoters
if (opt$promoters){
  # With promoters
  enrich <- computeCellTypeEnrichment(variant.by.cells,
                                      variant.list.filter,
                                      cell.type.annot,
                                      score.col=opt$variantScoreCol,
                                      min.score=opt$variantScoreThreshold,
                                      bg.vars=bgVars)
  # Without promoters
  enrich.nopromoter <- computeCellTypeEnrichment(getVariantByCellsTable(subset(filter.flat, !Promoter), opt$variantScoreCol),
                                                 subset(curr.vl, !Promoter),
                                                 cell.type.annot,
                                                 score.col=opt$variantScoreCol,
                                                 min.score=opt$variantScoreThreshold, bg.vars=bgVars)
  enrich <- merge(enrich, enrich.nopromoter %>% select(CellType,vsGenome.enrichment,vsGenome.log10pBinom,vsGenome.Significant), by="CellType", suffixes=c("",".NoPromoters"))
} else {
  enrich <- computeCellTypeEnrichment(variant.by.cells,
                                      variant.list.filter,
                                      cell.type.annot,
                                      score.col=opt$variantScoreCol,
                                      min.score=opt$variantScoreThreshold,
                                      bg.vars=bgVars)
}

write.tab(enrich, file=paste0(edir, "/Enrichment.CellType.vsScore.", trait,".tsv"))
cell.type.annot <- addEnrichmentSignificantCellTypes(cell.type.annot, enrich, trait)

# TODO: Assuming one trait
for (cat in cell.categories) {
  print(cat)
  ## Aggregate for each of the cell type categories
  plotCellTypeEnrichmentBarplot(enrich, cat, main=paste(trait, cat), sort.by.group=TRUE)
  plotCellTypeEnrichmentBarplot(enrich, cat, main=paste(trait, cat), sort.by.group=FALSE)
  
  # TODO: use score.col and bgVars
  enrich.grouped <- computeCellTypeEnrichment(variant.by.cells, variant.list.filter, 
                                              cell.type.annot, cell.group.by=cat, 
                                              score.col=opt$variantScoreCol,
                                              min.score=opt$variantScoreThreshold,
                                              bg.vars=bgVars) # NOT SUPPORTED YET, ldsc=ldsc[[trait]])
  write.tab(enrich.grouped, file=paste0(edir, "/Enrichment.CellType.vsScore.", trait, ".", cat, ".tsv"))    
  # TODO: change plot axis label
  plotCellTypeEnrichmentBarplot(enrich.grouped, "CellType", main=paste(trait, cat))
  # }
  #}, error = function(e) { print(e); print(paste0("Failed enrichment barplots for ", trait)) })
}
dev.off()


# Redundant?
#if (!is.null(ldsc)) cell.type.annot$`Binary.AnyDisease_FMOverlap_Enriched` <- apply(cell.type.annot[,paste0("Binary.",diseases,"_FMOverlap_Enriched"),drop=F], 1, any)

# What to do with these cell.bins?
cell.bins <- colnames(cell.type.annot)[grepl("Binary.", colnames(cell.type.annot))]
write.tab(cell.type.annot, file=paste0(opt$outbase,"CellTypes.Annotated.txt"))

saveProgress()

# TODO: if an option is provided on command line, create these barplots
pdf(file=paste0(opt$outbase, "OverlapGroupedByPosteriorProb.pdf"), height=6, width=8)
posterior.prob.breaks <- c(0,0.001,0.01,0.1,1)
for (cat in c(cell.categories, cell.bins)) {
  tryCatch({
    freq <- plotOverlapByPosteriorProb(variant.list.filter, filter.flat, posterior.prob.breaks, cell.type.annot, cat)
  }, error = function(e) print(paste0("Failed plotOverlapByPosteriorProb for category ", cat, "; ", e)))
}
dev.off()


##############################################################################
## A merged metric based on posterior probabilities of all variants ... this 
##  doesn't work as well
#tmp <- merge(variant.by.cells, variant.list[,c("variant","PosteriorProb")], by.x="QueryRegionName", by.y="variant")
#tmp <- subset(tmp, QueryRegionName %in% subset(variant.list, Disease == "IBD")$variant)
#res <- do.call(rbind, by(tmp, tmp$CellType, function(x) return(data.frame(CellType=x$CellType[1], sum=sum(x$PosteriorProb), weighted=sum(x$max.ABC * x$PosteriorProb)))))


##############################################################################
## Analysis per credible set
##  How many credible sets can we "explain" using various thresholds?
plotCredibleSetHeatmaps <- function(flat, cs, v.list, score.col=NULL, score.cutoff=NULL) {
  pdf(file=paste0(opt$outbase, "CredibleSetCellTypeHeatmap.pdf"), width=16, height=16, onefile=TRUE)
  # If a score and threshold are provided, subset
  to.plot <- flat
  if (!is.null(score.cutoff)){
    to.plot <- subset(flat, eval(parse(text = score.col)) >= score.cutoff)
  }
  tab <- plotCredibleSetCellTypeHeatmap(to.plot, main="all")
  sorted.cells <- colnames(tab)
  # Edit: assuming one disease
  #for (trait in traits) {
  #to.plot <- subset(flat, Disease==trait & eval(parse(text = score.col)) >= score.cutoff)
  if (nrow(to.plot) > 1) {
    tab <- plotCredibleSetCellTypeHeatmap(to.plot, main=trait)
  }
  dev.off()
  
  pdf(file=paste0(opt$outbase, "CredibleSetHeatmaps.pdf"), width=8, height=20, onefile=T)
  dir.create(paste0(opt$outbase, "CredibleSetHeatmaps"))
  for (cs.name in cs$CredibleSet) {
    tryCatch({
      plotGeneByCellTypeHeatmap(cs.name, subset(flat, PosteriorProb >= score.cutoff), v.list, genes, sorted.cells, main=cs.name, write.matrix=paste0(opt$outbase, "CredibleSetHeatmaps/",cs.name,".tsv"))
    }, error = function(e) paste0("Failed on ", cs.name))
  }
  dev.off()
  return(sorted.cells)
}

# If a variant score threhold was provided, running for significant variants.
# Else, using all variants
if (!is.null(sigScore.cs)){
  sorted.cells <- plotCredibleSetHeatmaps(sigScore.flat, sigScore.cs, variant.list.sigScore, score.col=opt$variantScoreCol, score.cutoff=opt$variantScoreThreshold)
} else {
  sorted.cells <- plotCredibleSetHeatmaps(all.flat, all.cs, variant.list)
}
#sorted.cells <- plotCredibleSetHeatmaps(filter.flat, filter.cs, variant.list.filter, 0.01, "pp01")


##############################################################################
## Output gene prioritization table
saveProgress()

pred.col.stats <- opt$predColForStats 

#pred.cols <- list(
#  DistanceRank=1,
#  ConnectionStrengthRank.Binary.AnyDisease_FMOverlap_Enriched=1,
#  `GeneList.Prediction.eRNA-Andersson2014`=TRUE,
#  `GeneList.Prediction.PCHiC-Javierre2016`=TRUE,
#  `GeneList.Prediction.ENCODE2012`=TRUE,
#  `GeneList.Prediction.Chun2017-eQTL`=TRUE,
#  `GeneList.Prediction.COGS0.5`=TRUE,
#  `GeneList.Prediction.COGS0.9`=TRUE,
#  `GeneList.Prediction.multiXcan.IBDCombined.AllCellTypes`=TRUE,
#  `GeneList.Prediction.OpenTargets.AllCellTypes`=TRUE,
#  `GeneList.Prediction.Chen2016.eQTL`=TRUE,
#  `GeneList.Prediction.RoadmapMergedLiu2017.Max`=TRUE,
#  `GeneList.Prediction.RoadmapMergedLiu2017`=TRUE,
#  `GeneList.Prediction.EnhancerAtlasGao2020`=TRUE,
#  `GeneList.Prediction.EnhancerAtlasGao2020.Max`=TRUE,
#  `GeneList.Prediction.Cao2017-JEME`=TRUE,
#  `GeneList.Prediction.Cao2017-JEME.Max`=TRUE,
#  `GeneList.Prediction.Sheffield2013`=TRUE,
#  `GeneList.Prediction.Sheffield2013.Max`=TRUE,
#  `GeneList.Prediction.Whalen2016-TargetFinder`=TRUE,
#  `GeneList.Prediction.Whalen2016-TargetFinder.Max`=TRUE,
#  `GeneList.Prediction.Granja2019`=TRUE,
#  `GeneList.Prediction.Granja2019.Max`=TRUE
#  )
#

# TODO: Do not add the gene lists
# TODO: use score instead of PP
# TODO: if a score thredhold is not provided, use all variants
if(!is.null(sigScore.cs)){
  gp.sigScore <- getGenePrioritizationTable(sigScore.flat, sigScore.cs, genes, genes.uniq, sorted.cells, cell.type.annot, cell.bins, var.score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold)
  #gp.sigScore <- addE2GMethodsToGP(gp.sigScore, alt.overlap, opt$variantScoreThreshold))
  write.tab(gp.sigScore, file=paste0(opt$outbase,"GenePredictions.tsv"))
  best.genes.sigScore <- getBestGenesFromPrioritizitionTable(gp.sigScore, pred.col.stats)
  write.tab(best.genes.sigScore, file=paste0(opt$outbase, "GenePredictions.Best2Genes.tsv"), col.names=F)
  
  # This is for comparing between prediction sets
  #gene.list.comparison.sigScore <- compareABCPredictionsToGeneLists(gp.sigScore, cell.bins, pred.cols=pred.cols[names(pred.cols) %in% colnames(gp.sigScore)])
  #write.tab(gene.list.comparison.sigScore, file=paste0(opt$outbase, "GenePredictions.GeneListComparison.tsv"))
  
  pair.stats <- getGeneCellTypePairAnalysis(sigScore.flat, gp.sigScore, "ConnectionStrengthRank", 1, cell.type.annot, unique.cell.col="Categorical.CellTypeMerged")
  stats$all <- c(stats$all, pair.stats)
  stats$all$`E-G distance for best 2 genes - min` <- quantile(subset(sigScore.flat, TargetGene %in% best.genes.sigScore & eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold)$distance, 0)
  stats$all$`E-G distance for best 2 genes - median` <- quantile(subset(sigScore.flat, TargetGene %in% best.genes.sigScore & eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold)$distance, 0.5)
  stats$all$`E-G distance for best 2 genes - max` <- quantile(subset(sigScore.flat, TargetGene %in% best.genes.sigScore & eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold)$distance, 1)
  stats$all$`E-G distance for best 2 genes - mean` <- mean(subset(sigScore.flat, TargetGene %in% best.genes.sigScore & eval(parse(text = opt$variantScoreCol)) >= opt$variantScoreThreshold)$distance)
  
  gp.all <- getGenePrioritizationTable(all.flat, all.cs, genes, genes.uniq, sorted.cells, cell.type.annot, cell.bins, var.score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold)
  #gp.all <- addE2GMethodsToGP(gp.all, alt.overlap, opt$posteriorProb)
  write.tab(gp.all, file=paste0(opt$outbase,"GenePredictions.all.tsv"))
  best.genes.all <- getBestGenesFromPrioritizitionTable(gp.all)
  write.tab(best.genes.all, file=paste0(opt$outbase, "GenePredictions.all.Best2Genes.tsv"), col.names=F)
  #gene.list.comparison <- compareABCPredictionsToGeneLists(gp.all, cell.bins, pred.cols=pred.cols[names(pred.cols) %in% colnames(gp.sigScore)])
  #write.tab(gene.list.comparison, file=paste0(opt$outbase, "GenePredictions.all.GeneListComparison.tsv"))
} else {
  gp.all <- getGenePrioritizationTable(all.flat, all.cs, genes, genes.uniq, sorted.cells, cell.type.annot, cell.bins, gene.lists, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold)
  #gp.all <- addE2GMethodsToGP(gp.all, alt.overlap, opt$posteriorProb)
  write.tab(gp.all, file=paste0(opt$outbase,"GenePredictions.all.tsv"))
  best.genes.all <- getBestGenesFromPrioritizitionTable(gp.all)
  write.tab(best.genes.all, file=paste0(opt$outbase, "GenePredictions.all.Best2Genes.tsv"), col.names=F)
  #gene.list.comparison <- compareABCPredictionsToGeneLists(gp.all, cell.bins, pred.cols=pred.cols[names(pred.cols) %in% colnames(gp.sigScore)])
  #write.tab(gene.list.comparison, file=paste0(opt$outbase, "GenePredictions.all.GeneListComparison.tsv"))
  
  pair.stats <- getGeneCellTypePairAnalysis(all.flat, gp.all, "ConnectionStrengthRank", 1, cell.type.annot, unique.cell.col="Categorical.CellTypeMerged")
  stats$all <- c(stats$all, pair.stats)
  stats$all$`E-G distance for best 2 genes - min` <- quantile(subset(all.flat, TargetGene %in% best.genes.all >= opt$variantScoreThreshold)$distance, 0)
  stats$all$`E-G distance for best 2 genes - median` <- quantile(subset(all.flat, TargetGene %in% best.genes.all >= opt$variantScoreThreshold)$distance, 0.5)
  stats$all$`E-G distance for best 2 genes - max` <- quantile(subset(all.flat, TargetGene %in% best.genes.all >= opt$variantScoreThreshold)$distance, 1)
  stats$all$`E-G distance for best 2 genes - mean` <- mean(subset(all.flat, TargetGene %in% best.genes.all >= opt$variantScoreThreshold)$distance)
  
}

## TO move
abcmax <- getABCMaxTable(gp.all, all.flat, cell.type.annot, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, cellTypeFlag=opt$predColForStats)
write.tab(abcmax, file=paste0(opt$outbase, "GenePredictions.ABCMaxSummary.tsv"))

## TO move
#abcmax <- getABCMaxTable(gp.all, all.flat, cell.type.annot, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, cellTypeFlag="AnyDisease_LDSC_Enriched")
#write.tab(abcmax, file=paste0(opt$outbase, "GenePredictions.ABCMaxSummary.LDSC.tsv"))

## Repeat using promoter activity-weighted score for ranking, rather than ABC score
#gp.sigScore.pw <- getGenePrioritizationTable(filter.flat, filter.cs, genes, genes.uniq, sorted.cells, cell.type.annot, cell.bins, gene.lists, min.pp=opt$posteriorProb, score.col="ABCWeightedByPromoter")
#write.tab(gp.sigScore.pw, file=paste0(opt$outbase,"GenePredictions.PromoterWeighted.tsv"))
#best.genes.sigScore.pw <- getBestGenesFromPrioritizitionTable(gp.sigScore.pw)
#write.tab(best.genes.sigScore.pw, file=paste0(opt$outbase, "GenePredictions.PromoterWeighted.Best2Genes.tsv"), col.names=F)
#gene.list.comparison.pw <- compareABCPredictionsToGeneLists(gp.sigScore.pw, cell.bins)
#write.tab(gene.list.comparison.pw, file=paste0(opt$outbase, "GenePredictions.PromoterWeighted.GeneListComparison.tsv"))

saveProgress()


#############################################################################
## To do: 
##  Count genes that have multiple independent variants pointing to them in the
##  same cell type
# TODO: get tissue annotation column from the argument parser, if needed
#tmp <- merge(filter.flat, cell.type.annot[,c("CellType","Categorical.IBDTissueAnnotations")])

if (!is.null(opt$tissueCategory)){
  tmp <- merge(filter.flat, cell.type.annot[,c("CellType", opt$tissueCategory)])
  top.genes.per.variant <- as.data.frame(tmp %>% group_by(QueryRegionName,CellType) %>% top_n(2, -Score.Fraction) %>% ungroup())
  top.genes.per.variant$LocusID <- unlist(lapply(strsplit(as.character(as.matrix(top.genes.per.variant$CredibleSet)),"-"),"[[",2))
  multiple.variants <- as.data.frame(unique(top.genes.per.variant[,c("CredibleSet","CellType","TargetGene")] %>% group_by(CellType,TargetGene)) %>% tally())
  multiple.variants <- merge(multiple.variants, unique(top.genes.per.variant[,c("TargetGene","LocusID")]))
  multiple.variants <- subset(multiple.variants, n>1) %>% arrange(TargetGene)
  write.tab(multiple.variants, file=paste0(opt$outbase, "GenePredictions.MultipleIndependentVariants.tsv"))
  
  # Same analysis but divided by cell category
  multiple.variants.cat <- as.data.frame(unique(top.genes.per.variant[,c("CredibleSet",opt$tissueCategory,"TargetGene")] %>% group_by(eval(parse(text = opt$tissueCategory)),TargetGene)) %>% tally())
  multiple.variants.cat <- merge(multiple.variants.cat, unique(top.genes.per.variant[,c("TargetGene","LocusID")]))
  multiple.variants.cat <- subset(multiple.variants.cat, n>1) %>% arrange(TargetGene)
  colnames(multiple.variants.cat) <- c("TargetGene", "TissueCategory", "n", "LocusID")
  write.tab(multiple.variants.cat, file=paste0(opt$outbase, "GenePredictions.MultipleIndependentVariants.ByTissueCategory.tsv"))
} else {
  tmp <- merge(filter.flat, cell.type.annot[,c("CellType")])
  top.genes.per.variant <- as.data.frame(tmp %>% group_by(QueryRegionName,CellType) %>% top_n(2, -Score.Fraction) %>% ungroup())
  top.genes.per.variant$LocusID <- unlist(lapply(strsplit(as.character(as.matrix(top.genes.per.variant$CredibleSet)),"-"),"[[",2))
  multiple.variants <- as.data.frame(unique(top.genes.per.variant[,c("CredibleSet","CellType","TargetGene")] %>% group_by(CellType,TargetGene)) %>% tally())
  multiple.variants <- merge(multiple.variants, unique(top.genes.per.variant[,c("TargetGene","LocusID")]))
  multiple.variants <- subset(multiple.variants, n>1) %>% arrange(TargetGene)
  write.tab(multiple.variants, file=paste0(opt$outbase, "GenePredictions.MultipleIndependentVariants.tsv"))
}

# TODO: if score threhold is not provided, use all variants
if (!is.null(sigScore.cs)){
  stats$all <- c(stats$all, list(
    `Number of loci with multiple independent credible sets (sigScore)`=sum(table(table(unlist(lapply(strsplit(as.character(as.matrix(sigScore.cs$CredibleSet)), "-"),"[[",2))))[-1]),
    `Number of those loci where multiple sets point to same gene in same cell type`=length(unique(multiple.variants$LocusID)),
    `Number of those loci where multiple sets point to same gene in similar cell type (e.g., myeloid)`=length(unique(multiple.variants.cat$LocusID))
  ))
} else {
  stats$all <- c(stats$all, list(
    `Number of loci with multiple independent credible sets (all)`=sum(table(table(unlist(lapply(strsplit(as.character(as.matrix(all.cs$CredibleSet)), "-"),"[[",2))))[-1]),
    `Number of those loci where multiple sets point to same gene in same cell type`=length(unique(multiple.variants$LocusID)),
    `Number of those loci where multiple sets point to same gene in similar cell type (e.g., myeloid)`=length(unique(multiple.variants.cat$LocusID))
  ))
}


#intersect(subset(multiple.variants, n>1)$TargetGene, best.genes.all)
# [1] "PDGFB"   "CBX7"    "TNFSF15" "DDC"     "IL2RA"   "NKX2-3"  "IL10"
# [8] "IKZF1"   "IL15RA"  "ADO"     "EGR2"    "RNF186"  "SBNO2"   "TMCO4"
# [15] "SLC26A3"

#############################################################################
## Metrics of cell-type specificity of enhancers and genes
# TODO: if a score is not provided, use all variants. Use the score column and 
# threholds instead of pp
if (!is.null(sigScore.flat)){
  cts <- getCellTypeSpecificityStats(sigScore.flat, best.genes.sigScore, gex, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold)
  stats$all <- c(stats$all, cts$stats)
} else {
  cts <- getCellTypeSpecificityStats(all.flat, best.genes.all, gex)
  stats$all <- c(stats$all, cts$stats)
}


#############################################################################
## Histograms of variant and credible set effects
# TODO: using score and threhold instead of pp
pdf(file=paste0(opt$outbase, "VariantHistograms.pdf"), width=5, height=4.5)
var.stats <- plotVariantHistograms(filter.flat, variant.list.filter, filter.cs, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, NULL)
stats$all <- c(stats$all, var.stats)
for (cat in cell.categories) {
  var.stats <- plotVariantHistograms(filter.flat, variant.list.filter, filter.cs, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, cat)
  stats$all <- c(stats$all, var.stats)
}
dev.off()

# ======================================
# Characteristics of enhancers with risk variants
# ======================================

# Plotting characteristics for relevant cell types and all cell types
enhancer.properties <- getEnhancerProperties(subset(filter.flat, CellType %in% relevant.cell.types), variant.list.filter)
enhancer.properties.all <- getEnhancerProperties(filter.flat, variant.list.filter)

# TODO: if the predictions are ABC, plot the ABC metrics. Else, plot some metrics
# from the general prediction file format

# If a score column, significance threhold, and a control threshold are provided,
# plotting elements overlapping the significant and control variants. If only
# a significance threshold is provided, plotting using the significant variants.
# Else, plotting using all variants. By default, plots are constructed using
# the PP and ctrlPP threholds.

# If score threholds are provided, must provide opt$scoreType for plot labels
pdf(file=paste0(opt$outbase, "EnhancerProperties.RelevantCellTypes.pdf"), width=4.5, height=5)
plotEnhancerProperties(enhancer.properties, score.type=opt$scoreType, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, ctrl.score=opt$variantCtrlScoreThreshold)
dev.off()
pdf(file=paste0(opt$outbase, "EnhancerProperties.AllCellTypes.pdf"), width=4.5, height=5)
plotEnhancerProperties(enhancer.properties.all, score.type=opt$scoreType, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, ctrl.score=opt$variantCtrlScoreThreshold)
dev.off()

# ======================================
# Write output tables
# ======================================

# If a score column and a significance threhold were provided...
if (!is.null(variant.list.sigScore)){
  write.tab(variant.list.sigScore, file=paste0(opt$outbase, "sig.variant.list.tsv"))
  write.tab(sigScore.flat, file=paste0(opt$outbase, "sig.flat.tsv"))
}

write.tab(variant.list.filter, file=paste0(opt$outbase, "filter.variant.list.tsv"))
write.tab(filter.flat, file=paste0(opt$outbase, "filter.flat.tsv"))

write.tab(variant.list, file=paste0(opt$outbase, "all.variant.list.tsv"))
write.tab(all.flat, file=paste0(opt$outbase, "all.flat.tsv"))

# Writing summary statistics
writeStats <- function() {
  getStatsToWrite <- function(stat.list) {
    data.frame(Text=names(stat.list), Value=unlist(stat.list))
  }
  for (name in names(stats)) 
    write.tab(getStatsToWrite(stats[[name]]), file=paste0(opt$outbase, "stats.", name, ".tsv"))
}
writeStats()

saveProgress()

# ======================================
# Plot performance of the gene predictors
# ======================================

# TODO: run this if comparing multiple prediction sets

# If a score and a threhold were provided, using those
#gp.plot <- gp.all
#if (!is.null(gp.sigScore)) gp.plot <- gp.sigScore
#
#pdf(file=paste0(opt$outbase, "GenePredictions.Enrichment.pdf"), width=5, height=5)
#plotGeneRankEnrichment(gp.plot, cell.type.annot, cell.bins, xlim=c(1,10))
#dev.off()
#
#pdf(file=paste0(opt$outbase, "GenePredictions.Enrichment.Max20.pdf"), width=5, height=5)
#plotGeneRankEnrichment(subset(gp.plot, DistanceRank <= 20), cell.type.annot, cell.bins, xlim=c(1,10))
#dev.off()
#
#pdf(file=paste0(opt$outbase, "GenePredictions.Enrichment.200Kb.pdf"), width=5, height=5)
#plotGeneRankEnrichment(subset(gp.plot, PromoterDistanceToBestSNP <= 200000), cell.type.annot, cell.bins, xlim=c(1,10))
#dev.off()



