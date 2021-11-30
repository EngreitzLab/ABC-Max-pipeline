##################################################################
## Jesse Engreitz
## May 14, 2017
## Annotate GWAS credible sets with enhancer-gene predictions
## Added to codebase Dec 9, 2018

##  This executable script takes in a set of variant predictions (variants overlapped with ABC),
##  plus information about a set (or subset) of these variants and their credible sets,
##  and outputs a series of useful data tables that can be processed or analyzed in various
##  ways by other scripts.

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

# TODO: match ABC.Score/ABC.Score throughout
# TODO: refactor so that finemapping info is optional

# Parse options from commans line
option.list <- list(
  make_option("--variants", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.variant.list.txt", help="File containing variants to consider"),
  make_option("--credibleSets", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.cs.txt", help="File containing credible set annotations"),
  make_option("--predictionFile", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/abc.tsv.gz", help="Prediction overlap file (.txt.gz) output from intersecting variants with predicted enhancers"),
  make_option("--methodName", type="character", default="ABC", help="Name of E-G prediction method"),
  make_option("--predScoreCol", type="character", default="ABC.Score", help="Name of the column with the prediction score, e.g. ABC-score."),
  make_option("--outbase", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_out/IBD/", help="Output file basename"),
  make_option("--outEnrichment", type="character", default="./Enrichment.CellType.vsScore.tsv", help="Output cell type enrichment table"),
  make_option("--outGenePredTable", type="character", default="./GenePredictions.allCredibleSets.tsv", help="Output gene prediction table filename"),
  make_option("--outGenePredTableDedup", type="character", default="./GenePredictions.allCredibleSets.tsv", help="Output gene prediction table filename"),
  make_option("--hasCellType", type="logical", default=TRUE, help="Do predictions have an associated cellType column?"),
  make_option("--hasTargetGeneTSS", type="logical", default=TRUE, help="Do predictions have an associated targetGeneTSS column?"),
  # TODO: edit downstream steps so that these are only used for ABC predictions
  make_option("--minPredScore", type="numeric", default=NA, help="Cutoff on prediction score for distal elements"),
  make_option("--minPredScorePromoters", type="numeric", default=NA, help="Cutoff on prediction score for tss/promoter elements"),
  make_option("--trait", type="character", help="Name of the trait or disease"),
  make_option("--variantScoreCol", type="character", default="PosteriorProb", help="Score determining variant significance. If set to NULL, all variants in the input file will be included"),
  # make_option("--scoreType", type="character", default="PP", help="Type of variantScoreCol, used in plot labels."),
  make_option("--variantScoreThreshold", type="numeric", default=0.1, help="Score cutoff for desired variants to analyze, e.g. PP>=0.1"),
  # make_option("--variantCtrlScoreThreshold", type="numeric", default=0.01, help="Score cutoff for for a control set of variants, e.g. PP<0.01"),
  # TODO: Move the background variant calculation to another script, because it needs to be run once for every prediction method, not once for every prediction method x trait combination
  make_option("--backgroundVariants", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/all.bg.SNPs.bed.gz", help="A set of background variants to use in the enrichment analysis. Bed format with chr, start, end, rsID"),
  make_option("--backgroundVariants_noPromoter", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/all.bg.SNPs.noPromoter.bed.gz", help="A set of background variants to use in the enrichment analysis. Bed format with chr, start, end, rsID"),
  make_option("--bgOverlap", type="character", default="/oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_out/ABC.OverlapAllSNPs.tsv.gz", help="Background variant overlap with predictions. Bed format with chr, start, end, rsID, enh-chr, enh-start, enh-end, CellType"),
  make_option("--bgOverlap_noPromoter", type="character", default="/oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_out/ABC.OverlapAllSNPs.noPromoter.tsv.gz", help="Background variant overlap with predictions. Bed format with chr, start, end, rsID, enh-chr, enh-start, enh-end, CellType"),
  make_option("--genes", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/RefSeqCurated.170308.bed", help="RefSeq gene BED file; this is to pull RefSeq IDs to determine coding/noncoding"),
  make_option("--genesUniq", type="character", default="/oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed", help="Collapsed RefSeq gene BED file used for E-G predictions"),
  make_option("--geneTSS", type="character", default="/oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed", help="TSS of genes file calculated"),
  make_option("--chr_sizes", type="character", help="Chromosome sizes"),
  make_option("--isEnhancerBed", type="logical", help="Calculating enrichment of an enhancer bedfile?"),

  make_option("--genePredMaxDistance", type="numeric", default=1000000, help="Gene prediction table: Include genes within this distance"),
  make_option("--biosampleEnrichThreshold", type="numeric", default=0.001, help="Gene prediction and enrichment tables: Bonferroni-adjusted p-value to call biosample as significantly enriched for overlapping variants (set to 1 to include all biosamples in gene prediction table)"),

  #make_option("--cellTypeTable", type="character", help="Table with annotations of cell types, with columns 'CellType', 'Categorical.*', 'Binary.*' for plotting enrichments"),
  #make_option("--relevantCellTypes", type="character", default="Binary.IBDRelevant", help="Column in cell type table containing mask for cell types that are 'relevant' to the trait"),
  #make_option("--tissueCategory", type="character", default="Categorical.IBDTissueAnnotations", help="Column in the cell type table containing tissue type categories"),
  #make_option("--predColForStats", type="character", default="ConnectionStrengthRank.Binary.IBDRelevant", help="Prediction to use for making stats"),
  #make_option("--gex", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/GeneTSSActivityQuantile.tsv", help="Filename with table of gene expression (or promoter activity) quantiles"),
  #make_option("--gexQuantileCutoff", type="numeric", default=0.4, help="Gene expression quantile cutoff to count gene as expressed"),
  #make_option("--promoterActivityRef", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/GeneList.txt", help="File from 1 cell type to use to extract distribution of promoter activity across genes (used for TSS-activity weighted predictions)"),
  #make_option("--cellTypeCov", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/covariance.all.tsv"),
  #make_option("--specificityBackground", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/allSpecificityScores.SparseMatrix.rds"),
  make_option("--removePromoterVariants", type="logical", default=FALSE, help="Remove credible sets with promoter variants from the filter.cs list"),
  make_option("--removeCodingVariants", type="logical", default=TRUE, help="Remove credible sets with coding variants from the filter.cs list"),
  make_option("--removeSpliceSiteVariants", type="logical", default=TRUE, help="Remove credible sets with splice site variants from the filter.cs list"),
  make_option("--codeDir", type="character", help="Directory to code base")
)

opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)
setwd(opt$outbase)

# Loading libraries and utilities
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "JuicerUtilities.R")))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "CredibleSetTools.R")))

# Function to save progress
saveProgress <- function() save.image(file=paste0(opt$outbase,"AnnotateCredibleSets.RData"))
saveProgress()

##############################################################################
## Load common data

# convert all boolean to caps strings
## TODO: Change these to booleans instead of string, DONE
opt$hasCellType <- as.logical(opt$hasCellType)
opt$hasTargetGeneTSS <- as.logical(opt$hasTargetGeneTSS)
genes.tss.file = opt$geneTSS
# All genes
## TODO: Move all of this gene processing logic somewhere else, so that the gene files input into this script are fully ready ; DONE
genes <- readBed(opt$genes)
genes <- addTSSToBED(genes)
# Genes used in the predictions
genes.uniq <- readBed(opt$genesUniq)
genes.uniq <- addTSSToBED(genes.uniq)
## Only include protein-coding genes
genes.uniq <- subset(genes.uniq, name %in% genes$name)
# Housekeeping genes to ignore
# hk.list <- read.delim(opt$housekeepingList, header=F, stringsAsFactors=F)[,1]

#  Vector of all diseases/traits in the variant list
# TODO: assuming one trait per variant list, DONE
variant.list <- read.delim(opt$variants, check.names=F)
#diseases <- unique(variant.list$Disease)

# Finding a vector of relevant cell types
# TODO: Requirements for this file? CellType | Categorical.IBDTissueAnnotations2
# TODO: only use this celltype table is using the ABC predictions. Else?
#if (opt$cellTypeTable != "nan") { 
#  cell.type.annot <- read.delim(opt$cellTypeTable, check.names=F)
#  cell.type.list <- cell.type.annot$CellType
#}

if (opt$hasCellType) {
  ## JME: This is not the best way to get this list of cell types
  ## KSM: Updated to grab celltypes from prediction file instead
  bgOverlap <- read.delim(opt$bgOverlap, check.names=F, header=F)
  cell.type.list <- unique(bgOverlap$V1)
} else {
  ## JME: I don't understand what this does?
  ## KSM: Initializing the variable 
  cell.type.list <- NULL
}


# All credible sets
all.cs <- read.delim(opt$credibleSets, check.names=F, stringsAsFactors=F)
all.cs$CredibleSet <- factor(all.cs$CredibleSet)
all.cs$MaxVariantScore <- sapply(all.cs$CredibleSet, function(cs) max(subset(variant.list, CredibleSet == cs)[,opt$variantScoreCol]))

all.cs$AnyCoding <- as.logical(all.cs$AnyCoding)
all.cs$AnySpliceSite <- as.logical(all.cs$AnySpliceSite)
all.cs$AnyPromoter <- as.logical(all.cs$AnyPromoter)

# Optionally removing all coding, splice site, and promoter variants. These sets
# are used in the enrichment analysis.
filter.cs <- subset(all.cs, (!AnyCoding | !opt$removeCodingVariants) &
                      (!AnySpliceSite | !opt$removeSpliceSiteVariants) &
                      (!AnyPromoter | !opt$removePromoterVariants))  ## Consider whether to factor this out of this script

variant.list.filter <- subset(variant.list, CredibleSet %in% filter.cs$CredibleSet)

# Finding CSs and variants that exceed the provided threhold. If no variant score
# and threshold are provided, using all variants in downstream steps
sigScore.cs = NULL
variant.list.sigScore = NULL
if (!(is.null(opt$variantScoreCol)) & !(is.null(opt$variantScoreThreshold))) {
  sigScore.cs <- subset(filter.cs, CredibleSet %in% subset(variant.list.filter, get(opt$variantScoreCol) >= opt$variantScoreThreshold)$CredibleSet)
  variant.list.sigScore <- subset(variant.list.filter, CredibleSet %in% sigScore.cs$CredibleSet) ## this is the list of variants in credible sets with at least 1 variant with posterior prob >10%
}

# Overlapping variants with predictions
# the required column names
overlap <- loadVariantOverlap(opt$predictionFile, genes.uniq, genes, variant.names=variant.list$variant, isTargetGeneTSS=opt$hasTargetGeneTSS)
write.table(overlap, file=paste0(opt$outbase, "data/overlap.full.tsv"))
overlap <- filterVariantOverlap(overlap, opt$predScoreCol, opt$minPredScore, opt$minPredScorePromoters, opt$hasTargetGeneTSS)
write.table(overlap, file=paste0(opt$outbase, "data/overlap.tsv"))

# Annotating overlaps
all.flat <- annotateVariantOverlaps(overlap, variant.list, all.cs)
if (opt$hasCellType){ 
  all.flat <- subset(all.flat, CellType %in% cell.type.list)  ## IMPORTANT CHANGE
}
filter.flat <- subset(all.flat, CredibleSet %in% filter.cs$CredibleSet)

if (!(is.null(opt$variantScoreCol))){
  filter.flat <- subset(filter.flat, opt$variantScoreCol>=opt$variantScoreThreshold)
}

# If a threshold is provided, getting annotations for significant variants
sigScore.flat = NULL
if (!(is.null(opt$variantScoreCol)) & !(is.null(opt$variantScoreThreshold))) {
  sigScore.flat <- subset(filter.flat, CredibleSet %in% sigScore.cs$CredibleSet)
}
# Finding all diseases/traits in the credible set table
# Assume that the cs has only one trait
#traits <- unique(all.cs$Disease)
trait <- opt$trait

# specify file names 
all.flat.file=paste0(opt$outbase, "data/all.flat.tsv")
filter.flat.file=paste0(opt$outbase, "data/filter.flat.tsv")
variant.list.filter.file=paste0(opt$outbase, "data/variant.filter.tsv")

# Creating an output directory and writing the result to files
dir.create(paste0(opt$outbase,"data/"))

filter.flat <- filter.flat[, c((1:ncol(filter.flat))[-1], 1)]
write.tab(all.flat, file=all.flat.file)
write.tab(filter.flat, file=filter.flat.file)
if (!(is.null(sigScore.flat))) {
  write.tab(sigScore.flat, file=paste0(opt$outbase, "data/sigScore.flat.tsv"))
}
variant.list.filter$pos_end <- variant.list.filter$position + 1
colnames <- names(variant.list.filter)
variant.list.filter <- variant.list.filter[, c(1, 2, length(colnames), 3:length(colnames)-1)]
write.tab(variant.list.filter, file=variant.list.filter.file)
# Reading in the gene TSS activity quantile information
# gex <- read.delim(opt$gex, check.names=F)


##############################################################################
## Calculate enrichment per cell type by comparing the significant/provided
## variants with a set of background variants

# Use the set of background variants instead of ctrlPP
# Remove the column from output that uses ctrlPP
#if (!is.null(opt$ctrlProb)) {
edir <- paste0(opt$outbase, "enrichment/")
dir.create(edir)

#pdf(file=paste0(opt$outbase, "Enrichment.CellType.vsScore.pdf"), width=5, height=5)

# Note: assuming one trait
#for (trait in unique(variant.list.filter$Disease)) {
#  tryCatch({
#curr.vl <- subset(variant.list.filter, Disease==trait)
variant.by.cells <- getVariantByCellsTable(filter.flat, isCellType=opt$hasCellType, isEnhancerBed=opt$isEnhancerBed)
#write.tab(variant.by.cells, file="variant.by.cells.tsv")
## Get Corresponding background variant, background Overlap and prediction File without regions that
# intersect with promoters
# A set of background variants
bgVars <- read.delim(opt$backgroundVariants, check.names=F, header=F)  
bgVars_noPromoter <- read.delim(opt$backgroundVariants_noPromoter, check.names=F, header=F)
bgOverlap <-read.delim(opt$bgOverlap, check.names=F, header=F)
bgOverlap_noPromoter <- read.delim(opt$bgOverlap_noPromoter, check.names=F, header=F)

# get filter.flat with no Promoters
noPromoters = paste0(opt$outbase, "data/filter.flat.noPromoters.tsv")
getNoPromoterPredictions(filter.flat.file, noPromoters, genes.tss.file)
variant.file.noPromoter <- read.table(noPromoters, header=T,fill = TRUE)
variant.by.cells.noPromoter <- getVariantByCellsTable(variant.file.noPromoter, isCellType=opt$hasCellType, isEnhancerBed=opt$isEnhancerBed)

# variant.list.filter with no Promoters
variant.list.noPromoter.file=paste0(opt$outbase, "data/variant.filter.noPromoters.tsv")
getNoPromoterPredictions(variant.list.filter.file, variant.list.noPromoter.file, genes.tss.file)
variant.list.noPromoters <- read.table(variant.list.noPromoter.file, header=T, fill = TRUE)

# With promoters
enrich <- computeCellTypeEnrichment(variant.by.cells,
                                    variant.list.filter,
                                    cell.type.list,
                                    trait,
                                    score.col=opt$variantScoreCol,
                                    min.score=opt$variantScoreThreshold,
                                    bg.vars=bgVars$V1,
                                    bg.overlap=bgOverlap$V2,
                                    isCellType=opt$hasCellType,
                                    enrichment.threshold=opt$biosampleEnrichThreshold)


# Without promoters
enrich.nopromoter <- computeCellTypeEnrichment(variant.by.cells.noPromoter,
                                               variant.list.noPromoters,
                                               cell.type.list,
                                               trait,
                                               score.col=opt$variantScoreCol,
                                               min.score=opt$variantScoreThreshold, 
                                               bg.vars=bgVars_noPromoter$V1, 
                                               bg.overlap=bgOverlap_noPromoter$V2, 
                                               isCellType=opt$hasCellType,
                                               enrichment.threshold=opt$biosampleEnrichThreshold)

enrich <- merge(enrich, enrich.nopromoter, by="CellType", suffixes=c("",".NoPromoters"))

print("Saving")
write.tab(enrich, file=opt$outEnrichment)

saveProgress()


####################################################################################
## Write gene predictions table
# TODO: need to modify to run on predictions that don't have a Score column

if (!opt$isEnhancerBed){
	enriched.cell.types <- subset(enrich, Significant)$CellType
	gp.all <- getGenePrioritizationTable(
	  all.flat, 
	  all.cs, 
	  genes, 
	  genes.uniq, 
	  enriched.cell.types, 
	  score.col=opt$predScoreCol, 
	  score.min=opt$minPredScore,
	  var.score.col=opt$variantScoreCol, 
	  var.score.min=opt$variantScoreThreshold,
	  max.distance=opt$genePredMaxDistance,
	  method.name=opt$methodName)

	writeGenePrioritizationTable(gp.all, file=paste0(opt$outGenePredTable))
	
	# dedup variantGenePairs
	all.flat.dedup <- dedupVariantGenePairs(all.flat)
	gp.dedup <- getGenePrioritizationTable(
	  all.flat.dedup,
          all.cs,
          genes,
          genes.uniq,
          enriched.cell.types,
          score.col=opt$predScoreCol,
          score.min=opt$minPredScore,
          var.score.col=opt$variantScoreCol,
          var.score.min=opt$variantScoreThreshold,
          max.distance=opt$genePredMaxDistance,
          method.name=opt$methodName)	
	writeGenePrioritizationTable(gp.dedup, file=paste0(opt$outGenePredTableDedup))
}
saveProgress()


