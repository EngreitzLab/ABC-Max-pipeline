## Jesse Engreitz
## December 9, 2018
## Functions to support GWAS - EP prediction analysis

library(dplyr)
library(tidyr)

## for dplyr < 0.7
pull <- function(x,y) {x[,if(is.name(substitute(y))) deparse(substitute(y)) else y, drop = FALSE][[1]]}

####################################################################################
## Functions to load in variant predictions from ABC

# TODO: edit so that colMap is no longer needed
loadVariantOverlap <- function(overlap.file, genes.uniq, genes, variant.names=NULL, colMap=NULL, overwriteTSS=FALSE) {
  ## Loads overlap file, and filters to variants in variant.names
  x <- read.delim(gzfile(overlap.file), check.names=F)
  
  if (!is.null(colMap)) {
    colMap <- unfactor(colMap)
    ## Added for remapping ABC columns from new ABC prediction format to column names expected by this R codebase JME 191221
    for (i in 1:nrow(colMap)) {
      oldCol <- colMap$FileColumns[i]
      newCol <- colMap$CodeColumns[i]
      if (oldCol == "") {
        x[newCol] <- NA
      } else if (oldCol %in% colnames(x)) {
        colnames(x)[colnames(x) == oldCol] <- newCol
      } else {
        stop(paste0("Mapping column names in loadVariantOverlap: Column name ", oldCol, " not found."))
      }
    }
  }
  
  if (!is.null(variant.names)) {
    tmp <- data.frame(variant=factor(as.character(as.matrix(variant.names)), levels=levels(x$QueryRegionName)))
    tmp <- subset(tmp, !is.na(variant))
    x <- merge(x, tmp, by.x="QueryRegionName", by.y="variant")
  }
  
  ## Merge and add various other annotations
  if (overwriteTSS) {
    if ("TargetGeneTSS" %in% colnames(x)) x <- x %>% select(-TargetGeneTSS)
    x <- merge(x, with(genes.uniq, data.frame(TargetGene=name, TargetGeneTSS=tss)), by="TargetGene")
  }
  #x$isOwnTSS <- with(x, TargetGeneTSS >= start & TargetGeneTSS <= end)
  
  codingSymbols <- subset(genes, grepl(";NM_", name))$symbol
  x$TargetGeneIsCoding <- as.character(as.matrix(x$TargetGene)) %in% codingSymbols
  
  return(x)
}


getCodingGenes <- function(genes=NULL, genes.uniq=NULL) {
  if (is.null(genes)) {
    genes <- readBed("/seq/lincRNA/data/hg19/RefSeqCurated.170308.bed")
    genes$symbol <- unlist(lapply(strsplit(as.character(as.matrix(genes$name)), ";"), "[", 1))
  }
  if (is.null(genes.uniq)) genes.uniq <- readBed("/seq/lincRNA/data/hg19/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed")
  codingSymbols <- as.character(as.matrix(subset(genes, grepl(";NM_", name))$symbol))
  result <- subset(genes.uniq, name %in% codingSymbols)
  result$tss <- result$start; result$tss[result$strand == "-"] <- result$end[result$strand == "-"]
  return(result)
}


filterVariantOverlap <- function(overlap, cutoff, tss.cutoff, hk.list) {
  ## Implements a different cutoff for distal enhancers versus distal promoters
  overlap <- overlap %>% filter( (((class != "tss" & class != "promoter") | isOwnTSS) & ABC.Score >= cutoff) | 
                                   ((class == "tss" | class == "promoter") & ABC.Score >= tss.cutoff)) %>% filter( !(TargetGene %in% hk.list) )
  return(overlap)
}


addPromoterWeightedPrediction <- function(overlap, promoter.activity.ref.file, accessibility.col='DHS.RPKM.TSS1Kb', remove.outliers=20) {
  ## Use the quantile metrics for promoter activity in a given cell type to index into a reference set of promoter activities (e.g. from K562 cells)
  ref <- read.delim(promoter.activity.ref.file)
  promoter.activity <- sort(sqrt((0.0001+ref$`H3K27ac.RPKM.TSS1Kb`)*(0.0001+ref[,accessibility.col])), decreasing=T)
  
  ## Chop off top 20 promoter outliers (cutting off bottom won't change anything, because we're already cutting out about bottom ~50% of genes from the predictions)
  promoter.activity <- promoter.activity[(remove.outliers+1):length(promoter.activity)]
  
  ## Rescale so max is 100
  promoter.activity <- promoter.activity * 100 / max(promoter.activity)
  
  ## Index into activity
  overlap$TargetGeneTSSActivity <- rev(promoter.activity)[round(length(promoter.activity)*overlap$TargetGenePromoterActivityQuantile)]
  
  ## This metric: range 0-100;  100 = enhancer is predicted to explain an absolute amount of transcription equal to most highly transcribed gene in genome
  overlap$ABCWeightedByPromoter <- with(overlap, ABC.Score * TargetGeneTSSActivity)
  
  return(overlap)
}


annotateVariantOverlaps <- function(overlap, variant.list, all.cs, var.cols=c("CredibleSet","Disease","PosteriorProb","Coding","SpliceSite","Promoter","LocusID")) {
  cols.to.remove <- which(colnames(variant.list) %in% c("chr","position","start","end"))  
  all.flat <- overlap
  all.flat <- merge(all.flat, variant.list[,-cols.to.remove], by.x="QueryRegionName", by.y="variant")
  return(all.flat)
}

loadGeneLists <- function(files) {
  do.call(c, lapply(strsplit(files,",")[[1]], loadOneGeneList))
}

loadOneGeneList <- function(file) {
  ## Returns list of gene sets (vectors of strings)
  if (is.null(file)) return(list())
  tab <- read.delim(file, comment.char="#", check.names=F, stringsAsFactors=F)
  lists <- list()
  for (col in colnames(tab)[grepl("GeneList.", colnames(tab))]) {
    lists[[col]] <- na.omit(tab[,col])
    lists[[col]] <- unique(lists[[col]][lists[[col]] != ""])
  } 
  return(lists)
}


#######################################################################################
## Functions for making and loading data from permutation tables written by MakeVariantCountTables.R

# Changed to use score col instead of PP
# Changed ABC.Score to ABC.score, is this correct?
# TODO: use column name based on selection
getVariantByCellsTable <- function(overlap, score.col="PosteriorProb") {
  write.tab(overlap, file="overlap.tsv")
  variant.by.cells <- overlap %>% group_by(QueryRegionName) %>% summarise( n.genes=n(), TargetGenes=paste(TargetGene,collapse=','), PosteriorProb=max(PosteriorProb) ) %>% as.data.frame()
  write.tab(variant.by.cells, file="variant.by.cells.tsv")
  return(variant.by.cells)
}

# Changed ABC.Score to ABC.score
# TODO: use column name based on selection
getVariantByGenesTable <- function(overlap) {
  variant.by.genes <- overlap %>% group_by(QueryRegionName,TargetGene) %>% summarise( n.cells=n(), CellTypes=paste0(CellType,collapse=',')) %>% as.data.frame()
  return(variant.by.genes)
}

getGenesByCellsTable <- function(overlap) {
  genes.by.cells <- overlap %>% group_by(TargetGene, CellType) %>% summarise( n.vars=n() ) %>% as.data.frame()
  return(genes.by.cells)
}


loadPermutationOverlapStats <- function(dir, nperm, filename="%.variants.by.cells.tsv.gz", FUN=countCellTypeOverlaps, ...) {
  result <- list()
  for (i in 1:nperm) {
    dat <- read.delim(gzfile(paste0(dir,"/",gsub("%",i, filename))))
    result[[i]] <- do.call(FUN, c(list(dat=dat), list(...)))
    result[[i]]$permutation <- i
  }
  return(result)
}


countCellTypeOverlaps <- function(dat, cell.type.list, variant.names=NULL, cell.groups=NULL, weightByPIP=FALSE) {
  cell.type.list <- na.omit(cell.type.list)
  
  if (!is.null(variant.names)) {
    dat <- dat %>% filter(QueryRegionName %in% variant.names)
  }
  
  if (nrow(dat) == 0) {
    result <- data.frame(CellType=cell.type.list, n=0)
  } else {
    
    if (!is.null(cell.groups)) {
      ## Convert cell types to cell type groups
      df <- data.frame(CellType=cell.type.list, CellGroup=cell.groups)
      dat <- merge(dat, df, by="CellType")
      dat$CellType <- dat$CellGroup
      cell.type.list <- unique(cell.groups)
    }
    
    dat$CellType <- refactor(dat$CellType, cell.type.list)
    write.tab(dat, file="dat.tsv")
    if (!weightByPIP) {
      dat <- unique(dat[,c("QueryRegionName","CellType")])
      result <- dat %>% group_by(CellType) %>% tally()
    } else {
      dat <- unique(dat[,c("QueryRegionName","CellType","PosteriorProb")])
      result <- dat %>% group_by(CellType) %>% summarise(n=sum(PosteriorProb))
    }
    result <- merge(data.frame(CellType=cell.type.list), result, all.x=T)
    result$n[is.na(result$n)] <- 0
  }
  return(result)
}


computeCellTypeEnrichmentByPermutation <- function(variant.names, label, cell.type.annot, outdir, permdir, nperm=1000, color.col='category1') {
  ## variant.names (factor)
  ## label (name for plotting)
  ## cell.type.annot (data.frame with 'cell_type' column)
  
  permute.celltypes <- loadPermutationOverlapStats(permdir, nperm, FUN=countCellTypeOverlaps, cell.type.list=cell.type.annot$cell_type, variant.names=variant.names)
  permute.celltypes <- permute.celltypes %>% rbind_all() %>% spread(key="permutation",value="n")
  
  pp10.celltypes <- countCellTypeOverlaps(variants.by.cells, cell.type.list=cell.type.annot$CellType, variant.names=variant.names)
  #  pp10.celltypes <- merge(cell.type.annot, pp10.celltypes, all.x=TRUE, by.x="cell_type", by.y="CellType")
  pp10.celltypes$permute.mean <- apply(permute.celltypes[,-1], 1, mean, na.rm=T)
  pp10.celltypes$permute.sd <- apply(permute.celltypes[,-1], 1, sd, na.rm=T)
  pp10.celltypes$permute.max <- apply(permute.celltypes[,-1], 1, max, na.rm=T)
  pp10.celltypes$q.permute <- with(pp10.celltypes, sapply(1:length(n), function(i) 1-ecdf(as.vector(as.matrix(permute.celltypes[i,-1])))(n[i])))
  pp10.celltypes$enrichment <- with(pp10.celltypes, n / permute.mean)
  pp10.celltypes$log10.p <- -with(pp10.celltypes, mapply(FUN=ppois, n, permute.mean, log.p=TRUE, lower.tail=F))/log(10)
  pp10.celltypes$log10.p[is.infinite(pp10.celltypes$log10.p)] <- NA
  pp10.celltypes <- pp10.celltypes %>% arrange(desc(log10.p))
  write.tab(pp10.celltypes, file=paste0(outdir, "/Enrichment.CellType.", label,".tsv"))
  return(to.plot)
}


computeCellTypeEnrichment <- function(variants.by.cells, variant.list, cell.type.annot, trait, cell.group.by=NULL, score.col="PosteriorProb", min.score=0.1, bg.vars=NULL, bg.overlap=NULL, noPromoter=FALSE) {
  ## Computes the enrichment of variants with high vs low posterior probabilities in each cell type
  ## variants.by.cells    data.frame output by getVariantByCellsTable
  ## variant.list         data.frame with all variant info. Will subset the variants.by.cells df by this list of variants
  ## cell.type.annot      data.frame containing cell types and columns with categorical annotations
  ## cell.group.by        column name from cell.type.annot to group the enrichment calculation
  ## hi.pp                Posterior probability threshold for 'high pp' variants
  ## lo.pp                Posterior probability threshold for 'low pp' variants
  
  # If a score threshold is probivided, selecting significant variants. Else, using all variants.
  # TODO: require promoter column?
  if (!(is.null(min.score))){
    hi.vars <- subset(variant.list, get(score.col) >= min.score)$variant
  } else {
    hi.vars <- variant.list$variant
  }
  
  #lo.vars <- subset(variant.list, (PosteriorProb < lo.pp) & (!noPromoter | !Promoter))$variant
  # Using bg variants instead of the lo.pp threshold
  #bg <- bg.vars
  #bg.vars <- bg.vars$V4
  
  stopifnot(nrow(hi.vars) > 0)
  stopifnot(nrow(bg.vars) > 0)
  
  message(cell.group.by)
  group.list <- if (!is.null(cell.group.by)) cell.type.annot[[cell.group.by]] else NULL
  #hi.count <- countCellTypeOverlaps(variants.by.cells, cell.type.list=cell.type.list, variant.names=hi.vars, cell.groups=group.list)
  #bg.count <- countCellTypeOverlaps(variants.by.cells, cell.type.list=cell.type.list, variant.names=bg.vars, cell.groups=group.list)
  if (!is.null(hi.vars)) {
       hi.count	 <- variants.by.cells %>% filter(QueryRegionName %in% hi.vars)
    }
  hi.count <- data.frame(n=length(hi.count$QueryRegionName))
  bg_V8 <- bg.overlap %>% distinct(V4)
  #weighted.count <- countCellTypeOverlaps(variants.by.cells, cell.type.list=cell.type.list, variant.names=variant.list$variant, cell.groups=group.list, weightByPIP=TRUE)
  write.tab(hi.count, file="TraitCount_NoCellType.tsv")
  write.tab(bg_V8, file="BackgroundCount_NoCellType.tsv")
  hi.count$n.ctrl <- length(unique(bg_V8$V4))
  hi.count$total <- length(hi.vars)
  hi.count$total.ctrl <- length(bg.vars$V4)
#  hi.count$total.non_unique <- length(bg.vars$V4)
  hi.count$Disease <- trait
  # What to use as prop.snps?
  hi.count$prop.snps <- with(hi.count, n/total)
  hi.count$enrichment <- with(hi.count, n/total / ( (n.ctrl)/total.ctrl ))
  hi.count$log10.p <- -with(hi.count, mapply(FUN=phyper, n, total, total.ctrl, n+n.ctrl, log.p=TRUE, lower.tail=F))/log(10)
  hi.count$log10.p[is.infinite(hi.count$log10.p)] <- NA
  hi.count$p <- with(hi.count, mapply(FUN=phyper, n, total, total.ctrl, n+n.ctrl, log.p=FALSE, lower.tail=F))
  hi.count$Significant <- p.adjust((hi.count$p), method="bonferroni") < 0.001
  hi.count$enrichment[with(hi.count, n==0 & n.ctrl==0)] <- NA
  #hi.count$n.weighted <- weighted.count$n
  #hi.count$total.weighted <- sum(variant.list[score.col])
#  hi.count$vsGenome.enrichment <- with(hi.count, n/total / LDSC.Prop._SNPs)
  hi.count$vsGenome.log10pBinom <- with(hi.count, pbinom(n, total, prop.snps, lower.tail=FALSE, log.p=TRUE)) / log(10)
  hi.count$vsGenome.Significant <- p.adjust(10^(hi.count$vsGenome.log10pBinom), method="bonferroni") < 0.001
  
   
  #if (!is.null(ldsc)) {
  #  ## TO DO: add support for cell.group.by and LDSC ... by adding Prop._SNPs
  #  stopifnot(is.null(cell.group.by))
  #  tomerge <- subset(ldsc, !is.na(CellType))
  #  colnames(tomerge)[colnames(tomerge) != "CellType"] <- paste0("LDSC.",colnames(tomerge)[colnames(tomerge) != "CellType"])
  #  hi.count <- merge(hi.count, tomerge, by="CellType", all.x=TRUE)
  #  hi.count$vsGenome.enrichment <- with(hi.count, n/total / LDSC.Prop._SNPs)
  #  hi.count$vsGenome.log10pBinom <- with(hi.count, pbinom(n, total, LDSC.Prop._SNPs, lower.tail=FALSE, log.p=TRUE)) / log(10)
  #  hi.count$vsGenome.Significant <- p.adjust(10^(hi.count$vsGenome.log10pBinom), method="bonferroni") < 0.001
  #  hi.count$vsGenome.PIPWeighted.enrichment <- with(hi.count, n.weighted/total.weighted / LDSC.Prop._SNPs)
  #}
  
  #hi.count <- merge(hi.count, cell.type.annot, all.x=TRUE, by="CellType")
  #hi.count <- hi.count %>% arrange(desc(log10.p))
  
  return(hi.count)
}


##################################################################################
## Code to compute cell type enrichments accounting for cell-type specificity of ABC
##  enhancer activity, and for similarities between cell types in the datasets
## Uses Mahalanobis distance, and is paired with scripts/CalculateCellTypeSimilarity.R

testMeanByPermutation <- function(obs, background, totalObs=length(obs), totalBackground=length(background), nPerm=100000, pseudocount=0.00001) {
  ## Permutation test for specificity scores (wilcox test does not perform well with many zero ranks;  t-test not powered because of sparse data)
  ## Returns permutation p-value for mean(obs) > mean(background)
  ## obs:  Observed values (e.g. specificity scores)
  ## background:  Background values (e.g. specificity scores for all 1000G variants... or all nonzero values, if padBackground is specified)
  ## totalBackground:  Supplement background values with zeros up to this total length
  ## nPerm:  Number of permutations to run
  
  ## in practice, sometimes the provided obs and background values also contain zeros: if all of the background is zero, then fail:
  if (all(background == 0)) return(NA)
  
  test.stat <- function(xNonZero,yNonZero,totalX=length(xNonZero),totalY=length(yNonZero)) (sum(c(xNonZero,pseudocount))/totalX) / (sum(c(yNonZero,pseudocount))/totalY)
  
  n.obs <- length(obs)
  z <- c(obs, background)
  n.allNonzero <- length(z)
  
  if (totalBackground == length(background) & totalObs == length(obs)) {
    permFun <- function(i) {
      indices <- sample.int(n.allNonzero, n.obs)
      test.stat(z[indices], z[-indices])
    }
    perm.dist <- sapply(1:nPerm, permFun)
    result <- mean(perm.dist > test.stat(obs, background))
  } else {
    ## This code takes the mean of a mostly zeroes vector, without ever allocating the vector to contain zeros
    totalAll <- totalBackground + totalObs
    
    permFun <- function(i) {
      indices.x <- sample.int(totalAll, totalObs)
      #indices.y <- (1:totalAll)[-indices.x]
      mask.x <- indices.x <= n.allNonzero
      vals.x <- z[indices.x[mask.x]]
      #mask.y <- indices.y <= n.allNonzero
      vals.y <- z[-indices.x[mask.x]]
      test.stat(vals.x,vals.y,totalObs,totalBackground)
    }
    perm.dist <- sapply(1:nPerm, permFun)
    result <- mean(perm.dist > test.stat(obs, background, totalObs, totalBackground))
  }
  return(result)
}

getSpecificityScores <- function(scores, cov.mat) {
  ## Manually computes mahalanobois distance, in a more efficient manner for this purpose than the existing R function
  rowz <- (as.matrix(scores) %*% cov.mat) * as.matrix(scores)
  specificityScores <- rowz / apply(rowz, 1, sum)
  specificityScores[is.nan(specificityScores)] <- 0
  specificityScores <- sqrt(specificityScores)
  stopifnot(all(range(specificityScores, na.rm=T) == c(0,1)))
  return(specificityScores)
}

getCellTypeSpecificityScores <- function(flat, cellTypeCov, binarize=FALSE) {
  stopifnot(is.factor(flat$CellType))
  sharedCellTypes <- intersect(levels(flat$CellType), colnames(cellTypeCov))
  flat <- subset(flat, CellType %in% sharedCellTypes)
  flat$CellType <- factor(as.character(as.matrix(flat$CellType)), levels=sharedCellTypes)
  flat$QueryRegionName <- factor(as.character(as.matrix(flat$QueryRegionName)), levels=unique(as.character(as.matrix(flat$QueryRegionName))))
  cellTypeCov <- cellTypeCov[sharedCellTypes,sharedCellTypes]
  
  mat <- flat %>% select(CellType, QueryRegionName, `activity_base`) %>% unique() %>% spread(CellType, `activity_base`, fill=0, drop=FALSE)
  rownames(mat) <- mat$QueryRegionName
  if (binarize) { mat <- mat %>% select(-QueryRegionName) %>% as.matrix(); mat[mat > 0] <- 1 }
  else mat <- mat %>% select(-QueryRegionName) %>% as.matrix() %>% asinh()
  specificityScores <- getSpecificityScores(mat, cellTypeCov)
  return(specificityScores)
}


getCellTypeSpecificityFreq <- function(specificityScores, variant.list) {
  cellTypeSums <- apply(specificityScores, 2, sum, na.rm=T)
  cellTypeFreq <- cellTypeSums / nrow(variant.list)
  result <- data.frame(CellType=colnames(specificityScores), SpecificitySum=cellTypeSums, SpecificityFreq=cellTypeFreq)
  return(result)
}


getCellTypeSpecificityEnrichments <- function(flat, variant.list, cellTypeCov, specificityBackground, nBackground=nrow(specificityBackground), hi.pp=0.1, lo.pp=0.01, nPerm=100, binarize=FALSE) {
  #### Compute hiPP vs loPP enrichments and p-values
  sharedCellTypes <- intersect(levels(flat$CellType), colnames(cellTypeCov))
  specificityScores.hi <- getCellTypeSpecificityScores(subset(flat, PosteriorProb >= hi.pp), cellTypeCov, binarize=binarize)
  specificityScores.lo <- getCellTypeSpecificityScores(subset(flat, PosteriorProb <= lo.pp), cellTypeCov, binarize=binarize)
  freq.hi <- getCellTypeSpecificityFreq(specificityScores.hi, subset(variant.list, PosteriorProb >= hi.pp))
  freq.lo <- getCellTypeSpecificityFreq(specificityScores.lo, subset(variant.list, PosteriorProb <= lo.pp))
  
  result <- merge(freq.hi, freq.lo, by="CellType", suffixes=c(".hi",".lo"))
  result$spec.hiVsLo.enrichment <- with(result, SpecificityFreq.hi/SpecificityFreq.lo)
  
  padZeros <- function(v, n) c(v, rep(0,n-length(v)))
  n.hiPp <- with(variant.list, sum(PosteriorProb >= hi.pp))
  n.loPp <- with(variant.list, sum(PosteriorProb <= lo.pp))
  
  pvals.hiVsLo <- data.frame(
    CellType=sharedCellTypes, 
    spec.hiVsLo.p.Wilcox=sapply(sharedCellTypes, function(ct) {
      hivals <- padZeros(specificityScores.hi[,ct], n.hiPp)
      lovals <- padZeros(specificityScores.lo[,ct], n.loPp)
      #  y <- c(hivals, lovals)
      #  A <- factor(rep(c('hi','lo'),c(n.hiPp,n.loPp)))
      wilcox.test(hivals, lovals, alternative="greater", na.rm=T, paired=TRUE)$p.value
      #  #wilcox_test(y~A, alternative="greater", distribution="exact") %>% pvalue() %>% as.numeric()
    }),
    spec.hiVsLo.pPerm=sapply(sharedCellTypes, function(ct)
      testMeanByPermutation(specificityScores.hi[,ct], specificityScores.lo[,ct], totalObs=n.hiPp, totalBackground=n.loPp, nPerm=nPerm))
  )
  #pvals.hiVsLo$spec.hiVsLo.fdr <- p.adjust(pvals.hiVsLo$spec.hiVsLo.pPerm)
  
  result <- merge(result, pvals.hiVsLo)
  
  #### Compute hiPP vs genome background enrichments and p-values
  res.vsGenome <- do.call(rbind, lapply(sharedCellTypes, function(ct) {
    #background.ct <- subset(specificityBackground, CellType == ct)$SpecificityScore
    background.ct <- specificityBackground[,ct]
    #background.ct <- background.ct[background.ct > 0]   ## This makes the testMeanByPermutation step go way faster
    #ct.string <- as.character(as.matrix(ct))
    currScores <- specificityScores.hi[,ct] #.string]
    hivals <- padZeros(currScores, n.hiPp)
    bkgroundsum <- (sum(background.ct)/nBackground)
    data.frame(CellType=ct,
               spec.GenomeBackground=bkgroundsum,
               spec.vsGenome.enrichment=(sum(currScores)/n.hiPp) / bkgroundsum,
               #spec.vsGenome.p.T=t.test(hivals, background.ct, alternative="greater")$p.value,
               spec.vsGenome.p.Wilcox=wilcox.test(hivals, background.ct, paired=TRUE, alternative="greater")$p.value)
    #vsGenome.pPerm=testMeanByPermutation(currScores, background.ct, totalObs=n.hiPp, totalBackground=nBackground, nPerm=nPerm))
  }))
  result <- merge(result, res.vsGenome)
  
  return(result)
}


###############################################################

plotCellTypeEnrichmentBarplots <- function(dat, color.cols, main="", ...) {
  for (color.col in color.cols) {
    plotCellTypeEnrichmentBarplot(dat, color.col, main=paste(main, color.col), ...)
  }
}

# TODO: fix log10.p
plotCellTypeEnrichmentBarplot <- function(dat, color.col, sort.by.group=TRUE, main="") {
  to.plot <- subset(dat, !is.na(get(color.col)) & !is.na(log10.p))
  
  suppressPackageStartupMessages(library(RColorBrewer))
  to.plot$color <- rainbow(length(levels(to.plot[[color.col]])))[as.numeric(to.plot[[color.col]])]
  
  if (sort.by.group) {
    to.plot <- to.plot %>% group_by(!!as.name(color.col)) %>% arrange(log10.p) %>% ungroup() %>% as.data.frame()
    to.plot$labels <- ""
    for (cat in unique(to.plot[[color.col]])) to.plot$labels[ round(median(which(to.plot[[color.col]] == cat))) ] <- cat
  } else {
    to.plot <- to.plot %>% arrange(log10.p) %>% as.data.frame()
    labels <- NULL
  }
  
  b <- with(to.plot, barplot(log10.p, col=color, border=NA, las=2, names.arg=labels, ylab='-log10 p_geom', main=main))
  abline(h=-log10(0.05/nrow(to.plot)), lty=2, col='gray')
  
  b <- with(to.plot, barplot(-log10(p.adjust(10^-log10.p)), col=color, border=NA, las=2, names.arg=labels, ylab='-log10 p_geom BH-corrected', main=main))
  abline(h=-log10(0.05), lty=2, col='gray')
  
  enrich.plot <- to.plot
  if (sort.by.group) enrich.plot <- to.plot %>% group_by(!!as.name(color.col)) %>% arrange(enrichment) %>% ungroup() %>% as.data.frame()
  else enrich.plot <- to.plot %>% arrange(enrichment) %>% as.data.frame()
  b <- with(enrich.plot, barplot(enrichment, col=color, border=NA, las=2, names.arg=labels, ylab='Enrichment (high vs low PP)', main=main))
  abline(h=1, lty=2, col='gray')
  
  if (all(c("vsGenome.log10pBinom","vsGenome.enrichment") %in% colnames(to.plot))) {
    if (sort.by.group) to.plot <- to.plot %>% group_by(!!as.name(color.col)) %>% arrange(-vsGenome.log10pBinom) %>% ungroup() %>% as.data.frame()
    else to.plot <- to.plot %>% arrange(-vsGenome.log10pBinom) %>% as.data.frame()
    b <- with(to.plot, barplot(-log10(p.adjust(10^vsGenome.log10pBinom)), col=color, border=NA, las=2, names.arg=labels, ylab='-log10 vsGenome p_binom BH-corrected', ...))
    abline(h=-log10(0.05), lty=2, col='gray')
    
    enrich.plot <- to.plot
    if (sort.by.group) enrich.plot <- to.plot %>% group_by(!!as.name(color.col)) %>% arrange(vsGenome.enrichment) %>% ungroup() %>% as.data.frame()
    else enrich.plot <- to.plot %>% arrange(vsGenome.enrichment) %>% as.data.frame()
    b <- with(enrich.plot, barplot(vsGenome.enrichment, col=color, border=NA, las=2, names.arg=labels, ylab='Enrichment (high PP vs genome)', main=main))
    abline(h=1, lty=2, col='gray')
  }
  
  if (all(c("LDSC.Enrichment_p","LDSC.Enrichment") %in% colnames(to.plot))) {
    if (sort.by.group) to.plot <- to.plot %>% group_by(!!as.name(color.col)) %>% arrange(-LDSC.Enrichment_p) %>% ungroup() %>% as.data.frame()
    else to.plot <- to.plot %>% arrange(-LDSC.Enrichment_p) %>% as.data.frame()
    b <- with(to.plot, barplot(-log10(LDSC.Enrichment_p), col=color, border=NA, las=2, names.arg=labels, ylab='-log10 LDSC p', ...))
    abline(h=-log10(0.05), lty=2, col='gray')
    
    enrich.plot <- to.plot
    if (sort.by.group) enrich.plot <- to.plot %>% group_by(!!as.name(color.col)) %>% arrange(LDSC.Enrichment) %>% ungroup() %>% as.data.frame()
    else enrich.plot <- to.plot %>% arrange(LDSC.Enrichment) %>% as.data.frame()
    b <- with(enrich.plot, barplot(LDSC.Enrichment, col=color, border=NA, las=2, names.arg=labels, ylab='Enrichment (%Heritability / %SNPs)', ...))
    abline(h=1, lty=2, col='gray')
  }
  return(to.plot)
}


plotOverlapByPosteriorProb <- function(variant.list, flat, posterior.prob.breaks, cell.type.annot=NULL, cat.col=NULL) {
  variant.list$PP.breaks <- cut(variant.list$PosteriorProb, posterior.prob.breaks)
  
  getTable <- function(vl, flat) {
    variant.list$ABC <- variant.list$variant %in% unique(flat$QueryRegionName)
    freq <- as.data.frame(table(subset(variant.list, ABC)$PP.breaks) / table(variant.list$PP.breaks))
    colnames(freq) <- c("PosteriorProb", "FractionOverlapping")
    freq$Enrichment <- freq$FractionOverlapping / freq$FractionOverlapping[1]
    freq
  }
  
  if (!is.null(cell.type.annot)) {
    ## Get the counts for each group of cell types specified by cat.col
    flat <- merge(flat, cell.type.annot[,c("CellType",cat.col)], by="CellType")
    freq <- do.call(rbind, lapply(unique(flat[,cat.col]), function(bin) {
      curr <- getTable(variant.list, subset(flat, get(cat.col) == bin))
      curr$Category <- bin
      curr
    }))
    
    if (grepl("Binary.", cat.col)) {
      ## Analysis for variants that are UNIQUE to each group and those that are SHARED
      newlabs <- do.call(rbind, by(flat[,c("QueryRegionName",cat.col)], flat$QueryRegionName, function(x) {
        labels=unique(x[,cat.col])
        if (length(labels) == 1) {
          return(data.frame(QueryRegionName=x$QueryRegionName[1], labels=paste0(labels[1],"-Only")))
        } else {
          return(data.frame(QueryRegionName=x$QueryRegionName[1], labels="Both"))
        }
      }, simplify=FALSE))
      flat <- merge(flat, newlabs, by="QueryRegionName")
      
      freq <- rbind(freq, do.call(rbind, lapply(unique(newlabs$labels), function(bin) {
        curr <- getTable(variant.list, subset(flat, labels == bin))
        curr$Category <- bin
        curr
      })))
    }
  } else {
    freq <- getTable(vl)
    freq$Category <- 'all'
  }
  
  #barplot(freq, xlab="Posterior Probability", ylab="Fraction overlapping ABC enhancer", beside=T)
  suppressPackageStartupMessages(library(ggplot2))
  g <- ggplot(freq, aes(factor(Category), FractionOverlapping, fill = PosteriorProb)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_brewer(palette = "OrRd") + theme_bw() + ggtitle(cat.col)
  print(g)
#  write.tab(freq, "OverlapGroupedByPosteriorProb.tsv")  
  g <- ggplot(freq, aes(factor(Category), Enrichment, fill = PosteriorProb)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_brewer(palette = "OrRd") + theme_bw() + ggtitle(cat.col)
  print(g)
  return(freq)
}



######################################################################
## Functions for plotting heatmaps

plotCredibleSetCellTypeHeatmap <- function(all.flat, return.raw=FALSE, histogram=TRUE, ...) {
  all.flat <- unfactor(all.flat)
  tab <- with(all.flat, table(CredibleSet, CellType))
  tab[tab > 0] <- 1
  tab <- as.matrix(tab)
  
  # add gene list to credible set names
  tmp <- unique(all.flat[,c("CredibleSet","TargetGene")])
  cs.gene.map <- with(tmp, tapply(TargetGene, CredibleSet, paste0, collapse=','))
  rownames(tab) <- paste0(rownames(tab),"  ",cs.gene.map[rownames(tab)])
  
  par(mar=c(12,0,4,8))
  tab.sorted <- plotCredibleSetCellTypeHeatmapHelper(tab, histogram=histogram, breaks=c(-1, 0.5, 1), col=c('white','red'), ...)
  
  if (return.raw) return(tab)
  else return(tab.sorted)
}

plotCredibleSetCellTypeHeatmapHelper <- function(tab, histogram=TRUE, margins=c(16,16), ...) {
  #tab <- matrix(0, nrow=length(unique(all.flat$CellType)), ncol=length(unique(all.flat$CredibleSet)))
  ## SUggest providing col (colors) and breaks (vector)
  suppressPackageStartupMessages(library(gplots))
  
  if (nrow(tab) < 2 | ncol(tab) < 2) return(tab)
  res <- heatmap.2(tab, trace='none', key=F, margins=margins, cex.main=0.8, cexCol=1, cexRow=1, ...)
  
  tab <- tab[res$rowInd, res$colInd]
  
  if (histogram) {
    par(mar=c(3,1,1,3))
    h.row <- apply(tab, 1, sum)
    h.col <- apply(tab, 2, sum)
    barplot(h.col, axes=T, ylim=c(0,max(h.col)), space=0, col='gray')
    barplot(h.row, axes=T, xlim=c(0,max(h.row)), space=0, col='gray', horiz=T)
  }
  return(tab)
}




plotGeneByCellTypeHeatmap <- function(cs.name, filter.flat, variant.list, genes, cell.list, buffer=1000000, write.matrix=NULL, ...) {
  cs.variants <- subset(variant.list, CredibleSet == cs.name)
  curr.pred <- subset(filter.flat, CredibleSet == cs.name)
  if (nrow(curr.pred) == 0) return()
  curr.pred$CellType <- as.character(as.matrix(curr.pred$CellType))
  
  ## Plot for all genes within 1 Mb of variants or out to furthest gene
  #plot.genes <- unique(subset(genes, chr==as.character(as.matrix(cs.variants$chr[1])) & 
  #  end >= min(min(cs.variants$position)-buffer,min(curr.pred$TargetGeneTSS)) & 
  #  start <= max(max(cs.variants$position)+buffer,max(curr.pred$TargetGeneTSS)))$symbol)
  plot.genes <- subset(genes, symbol %in% as.character(as.matrix(unique(curr.pred$TargetGene))))
  plot.genes <- plot.genes[!duplicated(plot.genes$symbol),]
  plot.genes <- plot.genes[order(plot.genes$start),]
  plot.genes <- plot.genes$symbol
  
  tab <- matrix(0, nrow=length(cell.list)+1, ncol=length(plot.genes)+1)
  rownames(tab) <- c(cell.list,""); colnames(tab) <- c(plot.genes,"")
  for (i in 1:nrow(curr.pred)) {
    target.gene <- as.character(as.matrix(curr.pred[i,"TargetGene"]))
    if (!(target.gene %in% plot.genes)) {
      cat(paste0("Found gene that doesn't match: ", target.gene, "\n"))
    } else {
      tryCatch({
        tab[curr.pred[i,"CellType"], target.gene] <- curr.pred[i,"ABC.Score"]*100
      }, error = function(e) { print("Caught in plotGeneByCellTypeHeatmap"); browser() })
    }
  }
  
  if (ncol(tab) < 2) {
    print(paste0("Skipping ",cs.name," because only 1 gene is predicted"))
  } else {
    plotCellTypeHeatmap(tab, ...)
    if (!is.null(write.matrix)) {
      write.table(tab, file=write.matrix, sep='\t', quote=F)
    }
  }
}


plotCellTypeHeatmap <- function(tab, ...) {
  ## Plot a heatmap with genes or variants on X-axis and cell types on Y-axis
  require(gplots)
  require(RColorBrewer)
  par(mar=c(8.1,4.1,4.1,8.1))
  palette <- colorRampPalette(c("white","red"))(n = 12)
  heatmap.2(tab, Rowv=F, Colv=F, trace='none', dendrogram='none', col=palette, key=F, 
            breaks=c(0:10,30,100), margins=c(6,12), cex.main=0.8, cexCol=1, cexRow=1, ...)
}






####################################################################################
## Gene Prioritization Table


# Using score.col and min.score instead of min.PP
# TODO: refactor
# "ABC.Score" = "ABC.Score"
getGenePrioritizationTable <- function(
  all.clean, all.cs, 
  genes, genes.uniq, 
  cell.types, cell.type.annot, cell.bins, 
  #gene.lists=list(),
  score.col="ABC.Score",
  var.score.col="PosteriorProb",
  min.score=0.1, 
  distance=1000000, 
  contact.col="hic_contact_pl_scaled_adj") {
  ## Make an integrated gene table for prioritization comparisons
  ## For each credible set, list each gene that is either within the specified distance or is
  ##  connected by ABC. Compute various prioritization metrics for each gene
  
  ## Only include protein-coding genes
  genes <- subset(genes, !grepl("NR_",name) & symbol %in% genes.uniq$name)
  
  dat <- do.call(rbind, lapply( as.character(as.matrix(all.cs$CredibleSet)), 
                                function(cs.name) {
                                  range <- as.numeric(subset(all.cs, CredibleSet == cs.name)[,c("start","end")] + c(-distance,distance))
                                  curr.chr <- as.character(as.matrix(subset(all.cs, CredibleSet == cs.name)$chr))
                                  curr.genes <- unique(subset(genes.uniq, (chr == curr.chr & tss >= range[1] & tss <= range[2]) | 
                                                                (name %in% as.character(as.matrix(unique(subset(all.clean, get(score.col) >= min.score & CredibleSet == cs.name)$TargetGene)))))$name)
                                  if (length(curr.genes) > 0) curr.genes <- curr.genes[sapply(curr.genes, function(gene) any(grepl(paste0(gene, ";NM_"), genes$name)))]
                                  if (length(curr.genes) == 0) {
                                    print(paste0("No nearby genes found for ",cs.name))
                                    return(NULL)
                                  }
                                  
                                  mat <- matrix(0, nrow=length(curr.genes), ncol=length(cell.types), dimnames=list(row=curr.genes, col=cell.types))
                                  mat.contact <- matrix(0, nrow=length(curr.genes), ncol=length(cell.types), dimnames=list(row=curr.genes, col=cell.types))
                                  for (i in with(all.clean, which(CredibleSet == cs.name & PosteriorProb  >= 0.1 & as.character(as.matrix(TargetGene)) %in% curr.genes))) {
                                    #if (! (all.clean$TargetGene[i] %in% rownames(mat))) print(all.clean$TargetGene[i])
                                    #if (! (all.clean$CellType[i] %in% colnames(mat))) print(all.clean$CellType[i])
                                    tryCatch({
                                      mat[as.character(as.matrix(all.clean$TargetGene[i])), as.character(as.matrix(all.clean$CellType[i]))] <- all.clean[i,score.col]
                                      mat.contact[as.character(as.matrix(all.clean$TargetGene[i])), as.character(as.matrix(all.clean$CellType[i]))] <- all.clean[i,contact.col]
                                    }, error = function(e) print(paste("failed on", i)))
                                  }
                                  
                                  annot <- data.frame(CredibleSet=cs.name, Disease=subset(all.cs, CredibleSet == cs.name)$Disease,
                                                      CodingSpliceOrPromoterVariants=0, AnyCoding=0, AnySpliceSite=0, AnyPromoter=0, TargetGene=curr.genes, PromoterDistanceToBestSNP=0, IsClosestGeneToBestSNP=0, ConnectionStrengthRank=NA, CellTypeCountRank=NA, ContactRank=NA)   
                                  # TO add: CodingVariant=NA, SpliceSiteVariant=NA, PromoterVariant=NA, PromoterDistanceToAnySNP=0, 
                                  ## For genes with no predictions, rank is assigned as NA.  For others, lower rank is better (best = 1)
                                  annot$PromoterDistanceToBestSNP <- sapply(curr.genes, function(gene) abs(subset(all.cs, CredibleSet == cs.name)$BestSNPPos - unique(subset(genes.uniq, name == gene)$tss)))
                                  annot <- annot %>% group_by(CredibleSet, add=FALSE) %>% mutate(DistanceRank=rank(PromoterDistanceToBestSNP)) %>% as.data.frame() 
                                  
                                  annot$GeneBodyDistanceToBestSNP <- sapply(curr.genes, function(gene) {
                                    pos <- subset(all.cs, CredibleSet == cs.name)$BestSNPPos
                                    tmpgene <- subset(genes.uniq, name == gene)[1,,drop=F]
                                    startdist <- tmpgene$start - pos
                                    enddist <- tmpgene$end - pos
                                    if (sign(startdist) != sign(enddist)) { 
		  		      return(0) }
                                    else {
                                      return(min(abs(startdist), abs(enddist)))
				    }
                                  })
                                  annot <- annot %>% group_by(CredibleSet, add=FALSE) %>% mutate(DistanceToGeneBodyRank=rank(GeneBodyDistanceToBestSNP)) %>% as.data.frame() 
                                  
                                  if (nrow(mat) > 0) {
                                    cnx <- apply(mat, 1, max)
                                    annot$ConnectionStrengthRank[cnx != 0] <- rank(-cnx, ties.method="min")[cnx != 0]
                                    cnx <- apply(mat, 1, function(row) sum(row > 0))
                                    annot$CellTypeCountRank[cnx != 0] <- rank(-cnx, ties.method="min")[cnx != 0]  ## This should really be normalized to 'meta-cell types'
                                    annot$MaxABC <- apply(mat, 1, max)
                                    
                                    annot$ContactMean <- apply(mat.contact, 1, mean)
                                    annot$ContactRank[annot$ContactMean != 0] <- rank(-annot$ContactMean, ties.method="min")[annot$ContactMean != 0]
                                    
                                    for (bin in cell.bins) {
                                      curr.cells <- intersect(colnames(mat), as.character(as.matrix(with(cell.type.annot, CellType[get(bin)]))))
                                      cnx <- apply(mat[,curr.cells,drop=F], 1, max)
                                      annot[[paste0("MaxABC.",bin)]] <- cnx
                                      annot[[paste0("ConnectionStrengthRank.",bin)]] <- NA
                                      annot[[paste0("ConnectionStrengthRank.",bin)]][cnx != 0] <- rank(-cnx, ties.method="min")[cnx != 0]
                                      annot[[paste0("CellTypeCount.",bin)]] <- apply(mat[,curr.cells,drop=F], 1, function(row) sum(row != 0))
                                    }
                                  } else {
			            print("nrow(mat) seems to be smaller than 0")
                                    return(NULL)
                                  }
                                  
                                  #for (listname in names(gene.lists)) {
                                  #  annot[[listname]] <- toupper(curr.genes) %in% toupper(gene.lists[[listname]])
                                  #}
                                  
                                  annot$CodingSpliceOrPromoterVariants <- with(subset(all.cs, CredibleSet == cs.name), any(AnyCoding, AnySpliceSite, AnyPromoter))
                                  annot$AnyCoding <- subset(all.cs, CredibleSet == cs.name)$AnyCoding
                                  annot$AnySpliceSite <- subset(all.cs, CredibleSet == cs.name)$AnySpliceSite
                                  annot$AnyPromoter <- subset(all.cs, CredibleSet == cs.name)$AnyPromoter
                                  annot$IsClosestGeneToBestSNP <- annot$TargetGene == as.character(as.matrix(subset(all.cs, CredibleSet == cs.name)$BestSNPNearestGene))
                                  cbind(annot, mat)
                                }
  ))
  dat
}


namedVectorsToTable <- function(...) {
  ## Example:  namedVectorsToTable( table(sample(10,100,replace=TRUE)), table(sample(10,100,replace=TRUE)) )
  lists <- list(...)
  for (l in names(lists)) names(lists[[l]])[is.na(names(lists[[l]]))] <- "NA"
  cols <- unique(unlist(lapply(lists, names)))
  mat <- matrix(NA, nrow=length(lists), ncol=length(cols), dimnames=list(names(lists), cols))
  for (l in names(lists)) {
    mat[l, names(lists[[l]])] <- lists[[l]]
  }
  return(mat)
}

getPredictedCellTypes <- function(gp, cellTypes) {
  apply(gp, 1, function(row) {
    scores <- as.numeric(row[cellTypes]); names(scores) <- cellTypes
    scores <- scores[scores > 0]
    scores <- sort(scores)
    paste0(names(scores), collapse=',')
  })
}

getPredictedVariants <- function(gp, all.flat, cellTypes, score.col="PosteriorProb", min.score=0.1) {
  tmp <- all.flat %>% filter(CellType %in% cellTypes & get(score.col) >= min.score) %>% group_by(CredibleSet, TargetGene) %>% summarise(Variants=paste0(unique(QueryRegionName), collapse=',')) %>% as.data.frame()
  tmp <- unfactor(tmp)
  gp$origOrder <- 1:nrow(gp)
  gp.tmp <- merge(gp, tmp, all.x=TRUE, by=c("CredibleSet","TargetGene"))
  return(gp.tmp$Variants[order(gp.tmp$origOrder)])
}

# cellTypeFlag="AnyDisease_FMOverlap_Enriched" aka "IBD_FMOverlap_Enriched"
# TODO: get correct cell type column from main script
getABCMaxTable <- function(gp, all.flat, cell.type.annot, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, cellTypeFlag="IBD_FMOverlap_Enriched") {
  gp.pred <- subset(gp, get(paste0("ConnectionStrengthRank.Binary.",cellTypeFlag)) == 1)
  enriched.cell.types <- as.character(as.matrix(subset(cell.type.annot, get(paste0("Binary.",cellTypeFlag)))$CellType))
  gp.pred$CellTypes <- getPredictedCellTypes(gp.pred, enriched.cell.types)
  gp.pred$Variants <- getPredictedVariants(gp.pred, all.flat, enriched.cell.types, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold)
  abcmax <- gp.pred[,c("Disease","TargetGene","CredibleSet","Variants","CellTypes")]
  abcmax <- abcmax[order(abcmax$TargetGene),]
  return(abcmax)
  }

getBestGenesFromPrioritizitionTable <- function(gp, pred.col.stats="ConnectionStrengthRank", top.ranks=1:2) {
  return(unique(subset(gp[c(pred.col.stats,"TargetGene")], get(pred.col.stats) %in% top.ranks)$TargetGene))
}

reshapeGenePrioritizationTable <- function(gp, all.cs, ranks=1:3) {
  gp.flat <- gp %>% gather(key="RankMethod", value="Rank", starts_with("ConnectionStrengthRank"), DistanceRank) %>% filter(Rank %in% ranks)
  gp.new <- gp.flat %>% select(CredibleSet, Disease, TargetGene, CodingSpliceOrPromoterVariants, RankMethod, Rank)
  result <- do.call(rbind, by(gp.new, paste0(gp.new$CredibleSet,gp.new$RankMethod), function(x) {
    for (rank in ranks) {
      i <- which(x$Rank == rank)
      x[,paste0("GeneRank",rank)] <- ifelse(length(i) > 0, paste0(x$TargetGene[i], collapse=','), NA)
    }
    x <- x %>% select(-TargetGene, -Rank)
    return(unique(x))
  }, simplify=FALSE)) 
  result <- result %>% arrange(RankMethod,CredibleSet)
  result <- merge(result, all.cs, by="CredibleSet")
  return(result)
}

compareABCPredictionsToGeneLists <- function(gene.pred.table, cell.bins=c(), pred.cols=list(
  IsClosestGeneToBestSNP=TRUE,
  ConnectionStrengthRank=1,
  ConnectionStrengthRank=1:2,
  ConnectionStrengthRank=1:3,
  CellTypeCountRank=1)) {
  
  ## Returns a table with one row per pair of (gene list, prediction method). Prediction methods are entries in pred.cols input list
  ##   nCtrl = # of genes in credible sets where at least one gene has a valid, non-NA prediction in the given pred.col
  ##   nList = # of genes in those credible sets that are also in the given gene list
  ##   nCtrlInTopRanks = # of genes in those credible sets that have any of the values specified in pred.col
  ##   nListInTopRanks = similar to above, for the genes in a given gene list
  ##   fractionCtrlInTopRanks
  ##   fractionListInTopRanks
  ##   enrichment
  ##   log10p.geom
  ##   PredictionMethod
  ##   allCsList = # of credible sets where at leasdt one nearby geene is in the given gene list
  ##   nCs = # of credible sets where at least one gene has a valid, non-NA prediction in the given pred.col
  ##   nCsList = # of above sets where at least one nearby gene is in the given gene list
  ##   nCsListInTopRanks = # of above credible sets where at least one ABC-linked gene in the given gene list has any of the values specified in pred.col
  ##   nLinkedExpectedInTopRanks = Expected # of credible sets connecting to known gene assuming random choices among genes linked with ABC
  ##   nRandomExpectedInTopRanks = Expected # of credible sets connecting to known gene assuming random choice among genes within 1 Mb
  
  ## Add more prediction columns
  for (bin in cell.bins) {
    pred.cols <- c(pred.cols, setNames(list(1, 1:2), rep(paste0("ConnectionStrengthRank.",bin),2)))
  }
  
  # TODO:
  gene.lists <- colnames(gene.pred.table)[grepl("GeneList.", colnames(gene.pred.table))]
  if (length(gene.lists) == 0) return(data.frame())
  if (all(apply(gene.pred.table[,gene.lists,drop=F], 2, sum, na.rm=T) == 0)) return(data.frame())
  
  all.stats <- do.call(rbind, lapply(1:length(pred.cols), function(i) {
    pc <- names(pred.cols)[i]
    cs.to.include <- unique(subset(gene.pred.table, !is.na(get(pc)))$CredibleSet)
    gp.to.include <- subset(gene.pred.table[,unique(c("CredibleSet","TargetGene",names(pred.cols),gene.lists))], CredibleSet %in% cs.to.include)
    
    ranks.ctrl <- table(gp.to.include[,pc], useNA="always")
    top.ranks <- as.character(as.matrix(pred.cols[[i]]))
    
    rank.stats <- do.call(rbind, lapply(gene.lists, function(gl) {
      ranks.gl <- table(subset(gp.to.include[,c(pc,gl)], get(gl))[,pc], useNA="always")
      nCtrl <- sum(ranks.ctrl, na.rm=T)
      nList <- sum(ranks.gl, na.rm=T)
      nRankCtrl <- sum(ranks.ctrl[top.ranks], na.rm=T)
      nRankList <- sum(ranks.gl[top.ranks], na.rm=T)
      
      gpt.slim <- gene.pred.table[,c(gl,"CredibleSet","TargetGene",pc)]  ## this speeds things up a lot; making subsets of large tables is time-consuming
      
      names.allCsList <- unique(subset(gpt.slim, get(gl))$CredibleSet)
      names.CsList <- unique(subset(gp.to.include[,c(gl,"CredibleSet")], get(gl))$CredibleSet)
      names.CsListInTopRanks <- unique(subset(gpt.slim, get(gl) & get(pc) %in% top.ranks)$CredibleSet)
      names.CsListTargetGenes <- unique(subset(gp.to.include[,c(pc,gl,"TargetGene")], get(gl))$TargetGene)
      names.TargetGenes <- unique(subset(gpt.slim, get(gl) & get(pc) %in% top.ranks)$TargetGene)
      nTargetGenes <- length(names.TargetGenes)
      
      tmp <- subset(gene.pred.table[,c(gl,pc,"CredibleSet")], CredibleSet %in% names.CsList)
      
      res <- data.frame(
        GeneList=gl,
        nCtrl=nCtrl,
        nList=nList,
        nCtrlInTopRanks=nRankCtrl,
        nListInTopRanks=nRankList,
        fractionCtrlInTopRanks=nRankCtrl/nCtrl,
        fractionListInTopRanks=nRankList/nList,
        enrichment=(nRankList/nList) / (nRankCtrl/nCtrl),
        log10p.geom=phyper(nRankList, nRankCtrl, nCtrl-nRankCtrl, nList, lower.tail=FALSE, log.p=TRUE)/log(10),
        nAnyLink=with(gene.pred.table, sum(get(gl) & !is.na(get(pc)))),
        allCsList = length(unique(subset(gpt.slim, get(gl))$CredibleSet)),
        nCs = length(cs.to.include),
        nCsList = length(names.CsList),
        nCsListInTopRanks = length(names.CsListInTopRanks),
        nCsAnyLink=length(unique(subset(gpt.slim, get(gl) & get(pc))$CredibleSet)),
        nRandomExpectedInTopRanks = sum( unlist(by(tmp, tmp$CredibleSet, function(x) min(length(top.ranks), nrow(x)) / nrow(x), simplify=FALSE)) ),
        nLinkedExpectedInTopRanks = sum( unlist(by(tmp, tmp$CredibleSet, function(x) min(length(top.ranks), sum(!is.na(x[,pc]))) / sum(!is.na(x[,pc])), simplify=FALSE)), na.rm=T),
        names.allCsList = paste0(names.allCsList, collapse=","),
        names.CsList = paste0(names.CsList, collapse=","),
        names.CsListInTopRanks = paste0(names.CsListInTopRanks, collapse=","),
        names.CsListTopGenes = paste0(names.CsListTargetGenes, collapse=","),
        names.TargetGenes = paste0(names.TargetGenes, collapse=","),
        nTargetGenes=nTargetGenes
      )
      return(res)
    }))
    rank.stats$PredictionMethod <- paste0(pc,'=',paste0(top.ranks,collapse=','))
    return(rank.stats)
  }))
  return(all.stats)
}


getGeneCellTypePairAnalysis <- function(flat, gene.pred.table, pred.col, top.ranks, cell.type.annot, unique.cell.col="Categorical.CellTypeMerged") {
  ## Compare the number of predicted gene/cell type pairs
  pc <- pred.col
  cs.to.include <- unique(subset(gene.pred.table, !is.na(get(pc)))$CredibleSet)
  
  gp <- subset(gene.pred.table, CredibleSet %in% cs.to.include & as.character(as.matrix(get(pred.col))) %in% top.ranks)
  flat <- merge(flat, gp[,c("CredibleSet","TargetGene")], by=c("CredibleSet","TargetGene"))
  flat <- merge(flat, cell.type.annot[,c("CellType",unique.cell.col)], by="CellType")
  flat <- unique(flat[,c("CredibleSet","TargetGene",unique.cell.col)])
  combos <- do.call(rbind, by(flat, flat$CredibleSet, function(f) with(f, data.frame(
    nGenes=length(unique(TargetGene)),
    csGenes=with(gene.pred.table, sum(CredibleSet == as.character(as.matrix(f$CredibleSet[1])))),
    nCellTypes=length(unique(get(unique.cell.col))),
    totalCellTypes=length(unique(cell.type.annot[,unique.cell.col])),
    nCombos=nrow(f))), simplify=FALSE))
  combos$totalCombos <- with(combos, csGenes * totalCellTypes)
  combos$fractionCombos <- with(combos, nCombos / totalCombos)
  
  stats <- with(combos, list(
    `Possible combinations of cell types and genes`=sum(totalCombos),
    `Prioritized combinations of cell types and genes`=sum(nCombos),
    `Fold-reduction in search space for gene-cell type combinations`=1 / (sum(nCombos)/sum(totalCombos))
  ))
  return(stats)
}


###################################################################
## Analysis of cell-type specificity of predictions
# Using score.col and min.score instead of pp
# If a score is not provided, using all variants
getCellTypeSpecificityStats <- function(flat, best.genes, gex, score.col=NULL, min.score=NULL) {
  
  # If no score threshold is provided, using all variants
  if (!is.null(score.col) & !is.null(min.score)){
    cts <- as.data.frame(subset(flat, get(score.col) >= min.score & TargetGene %in% best.genes) %>% group_by(QueryRegionName,TargetGene) %>% tally())
    ctsExpressedMid <- as.data.frame(subset(flat, get(score.col) >= min.score  & TargetGene %in% best.genes & TargetGenePromoterActivityQuantile >= 0.7) %>% group_by(QueryRegionName,TargetGene) %>% tally())
  } else {
    cts <- as.data.frame(subset(flat, TargetGene %in% best.genes) %>% group_by(QueryRegionName,TargetGene) %>% tally())
    ctsExpressedMid <- as.data.frame(subset(flat, TargetGene %in% best.genes & TargetGenePromoterActivityQuantile >= 0.7) %>% group_by(QueryRegionName,TargetGene) %>% tally())
  }
  cts <- merge(cts, ctsExpressedMid, by=c("QueryRegionName","TargetGene"), all.x=T)
  colnames(cts)[colnames(cts) == "n.x"] <- "nPredWhenExpressed"
  colnames(cts)[colnames(cts) == "n.y"] <- "nPredWhenExpressedMid"
  cts$nPredWhenExpressedMid[is.na(cts$nPredWhenExpressedMid)] <- 0
  cts$nGeneExpressed <- sapply(as.character(as.matrix(cts$TargetGene)), function(gene) sum(subset(gex, name == gene)[,-1] >= opt$gexQuantileCutoff))
  cts$nGeneExpressedMid <- sapply(as.character(as.matrix(cts$TargetGene)), function(gene) sum(subset(gex, name == gene)[,-1] >= 0.7))  ## 0.7 corresponds to about half of 'expressed' genes in a cell type
  
  res <- list(
    `Total cell types`=ncol(gex)-1,
    `Total variant-gene links considered for cell type specificity analysis`=nrow(cts),
    `Average number of cell types per variant-gene link`=mean(cts$nPredWhenExpressed),
    `Median number of cell types per variant-gene link`=median(cts$nPredWhenExpressed),
    `Average number of cell types per variant-gene link (mid+ expression only)`=mean(cts$nPredWhenExpressedMid),
    `Median number of cell types per variant-gene link (mid+ expression only)`=median(cts$nPredWhenExpressedMid),
    `Average number of cell types where gene is expressed`=mean(cts$nGeneExpressed),
    `Median number of cell types where gene is expressed`=median(cts$nGeneExpressed),
    `Average number of cell types where gene is expressed (mid+ expression only)`=mean(cts$nGeneExpressedMid),
    `Median number of cell types where gene is expressed (mid+ expression only)`=median(cts$nGeneExpressedMid))
  
  return(list(tab=cts, stats=res))
}



###################################################################
## Code for variant histograms

plotVariantHistograms <- function(flat, variant.list, cs, score.col, min.score, cell.group=NULL) {
  flat <- subset(flat, get(score.col) >= min.score)
  variant.list <- subset(variant.list, get(score.col) >= min.score)
  
  ## Redo the factor so that count table only includes these variants
  variant.list$variant <- factor(as.character(as.matrix(variant.list$variant)))
  flat$QueryRegionName <- refactor(flat$QueryRegionName, variant.list$variant)
  stopifnot(all(!is.na(flat$QueryRegionName)))
  
  if (is.null(cell.group)) {
    variant.by.cells <- getVariantByCellsTable(flat)
    genes.by.cells <- getGenesByCellsTable(flat)
  } else {
    flat <- merge(flat, cell.type.annot, by="CellType")
    variant.by.cells <- flat %>% group_by_(.dots=c("QueryRegionName",cell.group)) %>% summarise( n.genes=n() )
    genes.by.cells <- flat %>% group_by_(.dots=c("TargetGene", cell.group)) %>% summarise( n.vars=n() )
  }
  cells.per.variant <- table(variant.by.cells$QueryRegionName)
  
  variant.by.genes <- getVariantByGenesTable(flat)
  genes.per.variant <- table(variant.by.genes$QueryRegionName)
  
  main <- paste0("Grouping: ", cell.group)
  hist(cells.per.variant, border=NA, breaks=50, col='gray', xlab="# Cell Types/Groups per Variant", main=main)
  hist(cells.per.variant[cells.per.variant != 0], border=NA, breaks=50, col='gray', xlab="# Cell Types/Groups per Variant (no zeroes)", main=main)
  hist(genes.per.variant, border=NA, breaks=50, col='gray', xlab="# Genes per Variant", main=main)
  hist(genes.per.variant[genes.per.variant != 0], border=NA, breaks=50, col='gray', xlab="# Genes per Variant (no zeroes)", main=main)
  hist(table(as.character(as.matrix(genes.by.cells$TargetGene))), border=NA, breaks=50, col='gray', xlab="# Cell Types/Groups per Gene", main=main)
  
  n.cells <- if (is.null(cell.group)) length(levels(flat$CellType)) else length(levels(flat[,cell.group]))
  prefix <- paste0("[Min. score: ", min.score, "; Cell Grouping: ", cell.group, "] ")
  cell.suffix <- paste0(" (of ", n.cells, ")")
  stats <- list()
  stats[[paste0(prefix, "Median cell types/groups per variant, given variant is active in at least 1 cell type")]] <- paste0(median(cells.per.variant[cells.per.variant != 0]), cell.suffix)
  stats[[paste0(prefix, "Mean cell types/groups per variant, given variant is active in at least 1 cell type")]] <- paste0(mean(cells.per.variant[cells.per.variant != 0]), cell.suffix)
  stats[[paste0(prefix, "Median genes per variant, given variant is active in at least 1 cell type")]] <- median(genes.per.variant[genes.per.variant != 0])
  stats[[paste0(prefix, "Mean genes per variant, given variant is active in at least 1 cell type")]] <- mean(genes.per.variant[genes.per.variant != 0])
  stats[[paste0(prefix, "Fraction variants that regulate more than 1 gene, given variant regulates at least 1 gene")]] <- sum(genes.per.variant[genes.per.variant != 0] > 1) / sum(genes.per.variant != 0)
  return(stats)
}


###################################################################
## Properties of enhancers

getEnhancerProperties <- function(flat, variant.list.filter) {
  #flat$Activity <- with(flat, sqrt(normalized_dhs * normalized_h3k27ac))
  v.stats <- flat %>% group_by(QueryRegionName) %>% summarize(
    maxA=max(activity_base), 
    maxC=max(hic_contact_pl_scaled_adj), 
    maxABC=max(ABC.Score),
    nGenes=length(unique(TargetGene)), 
    nCellTypes=length(unique(CellType)), 
    minDistance=min(distance), 
    maxLength=max(end-start), 
    totalLinks=n(),
    PosteriorProb=PosteriorProb[1])
  v.stats <- as.data.frame(v.stats)
  return(v.stats)
}

plotEnhancerProperties <- function(v.stats, score.type=NULL, score.col=NULL, min.score=NULL, ctrl.score=NULL) {
  
  # Using threholds if they are provided
  hi <- v.stats
  if (!(is.null(score.col) & !is.null(min.score))){
    hi <- unfactor(subset(v.stats, get(score.col)  >= min.score))
  }
  
  lo <- NULL
  if (!is.null(ctrl.score)){
    lo <- unfactor(subset(v.stats, get(score.col)  < ctrl.score))
  }
  
  # TODO: a separate function for ABC and the general prediction file format
  # TODO: add a command line option to input the type of the score, e.g. "PP", DONE
  if (!is.null(lo)){
    plotCumulativeDistribution(hi$maxA, lo$maxA, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="Activity")
    plotCumulativeDistribution(hi$maxC, lo$maxC, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="Contact")
    plotCumulativeDistribution(hi$maxABC, lo$maxABC, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="ABC Score")
    plotCumulativeDistribution(hi$nGenes, lo$nGenes, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="Unique target genes")
    plotCumulativeDistribution(hi$nCellTypes, lo$nCellTypes, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="Unique cell types")
    plotCumulativeDistribution(hi$minDistance, lo$minDistance, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="Distance to closest target gene")
    plotCumulativeDistribution(hi$maxLength, lo$maxLength, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="Length of element")
    plotCumulativeDistribution(hi$totalLinks, lo$totalLinks, paste0(score.type,">=",min.score), paste0(score.type,"<",ctrl.score), main="Target genes x cell types")
  } else if (!(is.null(score.col) & !is.null(min.score))){
    plotCumulativeDistribution(hi$maxA, paste0(score.type,">=",min.score), main="Activity")
    plotCumulativeDistribution(hi$maxC, paste0(score.type,">=",min.score), main="Contact")
    plotCumulativeDistribution(hi$maxABC, paste0(score.type,">=",min.score), main="ABC Score")
    plotCumulativeDistribution(hi$nGenes, paste0(score.type,">=",min.score), main="Unique target genes")
    plotCumulativeDistribution(hi$nCellTypes, paste0(score.type,">=",min.score), main="Unique cell types")
    plotCumulativeDistribution(hi$minDistance, paste0(score.type,">=",min.score), main="Distance to closest target gene")
    plotCumulativeDistribution(hi$maxLength, paste0(score.type,">=",min.score), main="Length of element")
    plotCumulativeDistribution(hi$totalLinks, paste0(score.type,">=",min.score), main="Target genes x cell types")
  } else {
    plotCumulativeDistribution(hi$maxA, main="Activity")
    plotCumulativeDistribution(hi$maxC, main="Contact")
    plotCumulativeDistribution(hi$maxABC, main="ABC Score")
    plotCumulativeDistribution(hi$nGenes, main="Unique target genes")
    plotCumulativeDistribution(hi$nCellTypes, main="Unique cell types")
    plotCumulativeDistribution(hi$minDistance, main="Distance to closest target gene")
    plotCumulativeDistribution(hi$maxLength, main="Length of element")
    plotCumulativeDistribution(hi$totalLinks, main="Target genes x cell types")
  }
}


loadLdscResults <- function(file, CellTypes) {
  if (length(file) == 0) return(NULL)
  if (is.null(file) | file == "NULL") return(NULL)
  input <- read.delim(file, stringsAsFactors=F, check.names=F)
  
  if (all(c("Disease","Category","Enrichment_p","Prop._h2") %in% colnames(input))) {
    result <- as.list(by(input, input$Disease, function(x) {
      x$CellType <- factor(as.character(as.matrix(x$Category)), levels=CellTypes)
      x$Significant <- x$Enrichment_p < 0.05 & x$Enrichment > 0
      return(x)
    }))    
  }
  else if ("Disease" %in% colnames(input) & "Category" %in% colnames(input) & "Enrichment_p" %in% colnames(input)) {
    ## This is the format from the UKBiobank results
    ldsc <- read.delim("/seq/lincRNA/RAP/GWAS/191223_Partitioning/AvgHiC.ABC0.015.minus150/ldsc/dz/Liu2015-IBD.results.txt", stringsAsFactors=F, check.names=F)
    
    result <- as.list(by(input, input$Disease, function(x) {
      x <- merge(x, ldsc[,c("Category","Prop._SNPs")])
      x$CellType <- factor(as.character(as.matrix(x$Category)), levels=CellTypes)
      x$Significant <- x$Enrichment_p < 0.05 & x$Enrichment > 0
      return(x)
    }))
    names(result) <- unique(as.character(as.matrix(input$Disease)))
  } else if ("ldsc" %in% colnames(input) & "Disease" %in% colnames(input)) {
    result <- lapply(input$ldsc, function(file) {
      ldsc <- read.delim(file, check.names=F, stringsAsFactors=F)
      ldsc$Category <- ldsc$Category <- gsub("L2_0","",ldsc$Category)
      ldsc$CellType <- factor(as.character(as.matrix(ldsc$Category)), levels=CellTypes)
      ldsc$Significant <- ldsc$Enrichment_p < 0.05 & ldsc$Enrichment > 0
      return(ldsc)
    })
    names(result) <- as.character(as.matrix(input$Disease))
  } else {
    stop(paste0("Required columns not found in LDSC file: ",file))
  }
  return(result)
}


addLdscSignificantCellTypes <- function(cell.type.annot, ldsc, diseases) {
  if (is.null(ldsc)) return(cell.type.annot)
  diseases <- intersect(as.character(as.matrix(diseases)), names(ldsc))
  for (disease in diseases) {
    mask <- as.character(as.matrix(cell.type.annot$CellType)) %in% as.character(as.matrix(subset(ldsc[[disease]], Significant)$CellType))
    cell.type.annot[,paste0("Binary.",disease,"_LDSC_Enriched")] <- mask
  }
  cell.type.annot[,paste0("Binary.AnyDisease_LDSC_Enriched")] <- apply(cell.type.annot[,paste0("Binary.",diseases,"_LDSC_Enriched"),drop=F], 1, any)
  return(cell.type.annot)
}

addEnrichmentSignificantCellTypes <- function(cell.type.annot, enrich, trait) {
  if ("Significant" %in% colnames(enrich)) {
    mask <- as.character(as.matrix(cell.type.annot$CellType)) %in% as.character(as.matrix(subset(enrich, Significant)$CellType)) 
    cell.type.annot[,paste0("Binary.",trait,"_FMOverlap_Enriched")] <- mask
  }
  return(cell.type.annot)
}


plotGeneRankEnrichment <- function(gp, cell.type.annot, cell.bins, xlim=c(1,10)) {
  ## TO do: Integrate this function with compareABCPredictionsToGeneLists (which computes many of the same statistics)
  
  ## loop over cell type annots and cell.bins and rank columns
  suppressPackageStartupMessages(library(ggplot2))
  print("Here") 
  #rank.cols <- colnames(gp)[greplany(c("ConnectionStrengthRank.","CellTypeCount.","MaxABC."), colnames(gp))]
  rank.cols <- colnames(gp)[greplany(c("ConnectionStrengthRank","DistanceRank"), colnames(gp))]
  cat.cols <- colnames(gp)[grepl("GeneList.", colnames(gp)) & !grepl("Smillie", colnames(gp))]
  pred.cols <- colnames(gp)[grepl("GeneList.Prediction.", colnames(gp))]
  print(rank.cols)
  print(pred.cols)
  print(cat.cols) 
  #moved to getGenePrioritizationTable
  #gp <- gp %>% group_by(CredibleSet, add=FALSE) %>% mutate(DistanceRank=rank(PromoterDistanceToBestSNP)) %>% as.data.frame() 
  
  int_breaks <- function(x, n = 3) pretty(x, n)[pretty(x, n) %% 1 == 0]
  #pr <- data.frame()
  print(cat.cols) 
  for (cat.col in cat.cols) {
    ## Get credible sets with at least one "true" nearby
    #csWithCategory <- unique(subset(gp, get(cat.col))$CredibleSet)
    ## Get credible sets with exactly 1 "true" nearby
    csWithCategory <- table(subset(gp, get(cat.col))$CredibleSet)
    csWithCategory <- names(csWithCategory[csWithCategory == 1])
    
    distanceRank.allcs <- getGeneRankEnrichment(gp, "DistanceRank", cat.col, cs.list=csWithCategory, Method="Distance")
    distanceGeneBody.allcs <- getGeneRankEnrichment(gp, "DistanceToGeneBodyRank", cat.col, cs.list=csWithCategory, Method="DistanceToGeneBody")
    print(length(pred.cols)) 
    ## Plot performance of other provided predictors
    if (length(pred.cols)>0) {
      pred <- do.call(rbind, lapply(pred.cols, function(col) getGenePredictionStats(gp, col, cat.col, csWithCategory)))
      write.tab(pred, "PrecisionRecallTable.tsv")
      r <- ggplot(pred, aes(x=Recall, y=Precision, color=PredictionMethod)) + geom_point(size=3)
      r <- r + xlim(0,1) + ylim(0,1)
      r <- r + theme_classic()
      r <- r + theme(plot.title=element_text(size = 10),
                     aspect.ratio=1)
      r <- r + ggtitle(paste0("CustomPredictors\n",cat.col))
      print(r)
    }
    
    #browser()
    #pred$Category <- cat.col
    #pr <- rbind(pr, pred[,c("Category","PredictionMethod","Precision","Recall")])
    
    ## Plot performance of rank predictions (ABC ConnectionStrengthRank)
    for (rank.col in rank.cols) {
      curr.cs <- getCsToEvaluateRanking(gp, rank.col, cat.col, exactlyOne=TRUE)
      all.ranks <- getGeneRankEnrichment(gp, rank.col, cat.col, curr.cs, Method="ABC") ## fill out here
      if (is.null(all.ranks)) next
      
      #pdf(file=paste0(opt$outbase, "GenePrediction.Enrichment.pdf"), width=5, height=5)
      toplot <- subset(all.ranks, rank>=xlim[1] & rank<=xlim[2])
      p <- ggplot(toplot, aes(x=rank, y=Enrichment.Category)) + geom_bar(stat='identity')
      p <- p + theme_classic()
      p <- p + xlim(xlim)
      p <- p + ylab("Enrichment")
      p <- p + scale_x_discrete(name="Gene Rank", breaks=seq(xlim[1],xlim[2]), labels=waiver()) + scale_y_continuous(breaks = int_breaks)
      p <- p + theme(text = element_text(size = 15), plot.title=element_text(size = 10))
      p <- p + ggtitle(paste0(rank.col,"\n",cat.col))
      print(p)
      
      toplot2 <- toplot %>% gather(key="Count",value="value",Freq.Random,Freq.InCategory)
      toplot2$Count <- factor(as.character(as.matrix(toplot2$Count)), levels=c("Freq.InCategory","Freq.Random"), labels=c("ABC","Random"))
      q <- ggplot(toplot2, aes(x=rank, y=value, fill=Count)) + geom_bar(stat='identity', position=position_dodge(), width=0.75)
      q <- q + theme_classic()
      q <- q + xlim(xlim)
      q <- q + ylab("Count")
      q <- q + scale_x_discrete(name="Gene Rank", breaks=seq(xlim[1],xlim[2]), labels=waiver()) + scale_y_continuous(breaks = int_breaks)
      q <- q + theme(text = element_text(size = 15), plot.title=element_text(size = 10))
      q <- q + scale_fill_manual(values=c('ABC'='red', 'Random'='black'))
      q <- q + ggtitle(paste0(rank.col,"\n",cat.col))
      print(q)
      
      distance.ranks <- getGeneRankEnrichment(gp, "DistanceRank", cat.col, cs.list=curr.cs, "Distance")
      geneBody.ranks <- getGeneRankEnrichment(gp, "DistanceToGeneBodyRank", cat.col, cs.list=curr.cs, "DistanceToGeneBody")
      toplot3 <- rbind(toplot, distance.ranks, geneBody.ranks)
      toplot3 <- subset(toplot3, rank>=xlim[1], rank<=xlim[2])
      r <- ggplot(toplot3, aes(x=Recall, y=Precision, color=Method)) + geom_point(size=3)
      r <- r + xlim(0,1) + ylim(0,1)
      r <- r + theme_classic()
      colors <- c('red','black','gray'); names(colors) <- c("ABC", "Distance", "DistanceToGeneBody")
      r <- r + scale_color_manual(values=colors)
      r <- r + geom_hline(yintercept=min(distance.ranks$Precision), linetype='dashed', color='gray')
      r <- r + theme(text = element_text(size = 20), 
                     plot.title=element_text(size = 10),
                     aspect.ratio=1)
      r <- r + ggtitle(paste0(rank.col,"\n",cat.col,"\n(", length(curr.cs), " Credible Sets)"))
      print(r)
      
      abcRank.allcs <- getGeneRankEnrichment(gp, rank.col, cat.col, cs.list=csWithCategory, Method="ABC")
      
      toplot3 <- rbind(abcRank.allcs, distanceRank.allcs, distanceGeneBody.allcs)
      toplot3 <- subset(toplot3, rank>=xlim[1], rank<=xlim[2])
      r <- ggplot(toplot3, aes(x=Recall, y=Precision, color=Method)) + geom_point(size=3)
      r <- r + xlim(0,1) + ylim(0,1)
      r <- r + theme_classic()
      colors <- c('red','black','gray'); names(colors) <- c("ABC", "Distance", "DistanceToGeneBody")
      r <- r + geom_hline(yintercept=min(distanceRank.allcs$Precision), linetype='dashed', color='gray')
      r <- r + scale_color_manual(values=colors)
      r <- r + theme(text = element_text(size = 20), 
                     plot.title=element_text(size = 10),
                     aspect.ratio=1)
      r <- r + ggtitle(paste0(rank.col,"\n",cat.col,"\n(All credible sets: ",length(csWithCategory),")"))
      print(r)
      
      #dev.off()
    }
  }
}


getCsToEvaluateRanking <- function(gp, rank.col, cat.col, exactlyOne=FALSE) {
  if (!exactlyOne) {
    unique(gp$CredibleSet)[sapply(unique(gp$CredibleSet), 
                                  function(cs) with(subset(gp[,c("CredibleSet",cat.col,rank.col)], CredibleSet == cs), 
                                                    any(get(cat.col), na.rm=T) & any(get(rank.col), na.rm=T)))]  
  } else {
    unique(gp$CredibleSet)[sapply(unique(gp$CredibleSet), 
                                  function(cs) with(subset(gp[,c("CredibleSet",cat.col,rank.col)], CredibleSet == cs), 
                                                    sum(get(cat.col), na.rm=T) == 1 & any(get(rank.col), na.rm=T)))]      
  }
}


getGeneRankEnrichment <- function(gp, rank.col="ConnectionStrengthRank.Binary.IBDRelevant", cat.col="GeneList.JME_Combined", cs.list=NULL, Method=NA) {
  if (is.null(cs.list)) {
    curr.cs <- getCsToEvaluateRanking(gp, rank.col, cat.col, exactlyOne=TRUE)
    #    curr.cs <- unique(subset(gp, get(cat.col) & !is.na(get(rank.col)))$CredibleSet)    
  } else {
    curr.cs <- cs.list
  }
  
  if (length(curr.cs) == 0) return(NULL)
  
  cat.ranks <- data.frame(table(subset(gp, get(cat.col) & CredibleSet %in% curr.cs)[,rank.col]))
  
  any.ranks <- data.frame(table(subset(gp, CredibleSet %in% curr.cs)[,rank.col]))
  colnames(any.ranks) <- c("rank","Freq.Any")
  
  random.ranks <- gp %>% filter(CredibleSet %in% curr.cs) %>% group_by(CredibleSet) %>% transmute(csRank=1:length(CredibleSet),weightedRank=1/length(CredibleSet))
  weighted.ranks <- random.ranks %>% group_by(csRank) %>% summarise(sum(weightedRank)) %>% as.data.frame()
  colnames(weighted.ranks) <- c("rank","Freq.Random")
  #table(random.ranks$csRank)
  
  all.ranks <- merge(unfactor(any.ranks), unfactor(weighted.ranks), all=TRUE, by="rank")
  if (nrow(cat.ranks) > 0) {
    colnames(cat.ranks) <- c("rank","Freq.InCategory")
    all.ranks <- merge(all.ranks, unfactor(cat.ranks), all=TRUE, by="rank")
  } else {
    all.ranks$Freq.InCategory <- 0
  }
  all.ranks$rank <- as.numeric(all.ranks$rank)
  all.ranks <- all.ranks[order(all.ranks$rank),]
  
  all.ranks$Enrichment.Category <- with(all.ranks, Freq.InCategory / Freq.Random)
  all.ranks[is.na(all.ranks)] <- 0
  all.ranks$Precision <- with(all.ranks, cumsum(Freq.InCategory) / cumsum(Freq.Any))
  all.ranks$Recall <- with(all.ranks, cumsum(Freq.InCategory)/max(sum(Freq.InCategory),length(cs.list)))
  all.ranks$Method <- Method
  return(all.ranks)
}


getGenePredictionStats <- function(gp, pred.col, cat.col="GeneList.JME_Combined", cs.list=NULL) {
  "Evaluate how well a binary predictor identifies the genes in the cat.col"
  
  if (is.null(cs.list)) {
    curr.cs <- getCsToEvaluateRanking(gp, rank.col, cat.col, exactlyOne=TRUE)
    #    curr.cs <- unique(subset(gp, get(cat.col) & !is.na(get(rank.col)))$CredibleSet)    
  } else {
    curr.cs <- cs.list
  }
  
  if (length(curr.cs) == 0) return(NULL)
  
  touse <- gp %>% filter(CredibleSet %in% curr.cs)
  
  result <- with(touse, data.frame(
    PredictionMethod=gsub("GeneList.Prediction.","",pred.col), 
    Precision=sum(get(pred.col) & get(cat.col)) / sum(get(pred.col)),
    Recall=sum(get(pred.col) & get(cat.col)) / sum(get(cat.col))))
  
  return(result)
}


getGeneListsFromE2GOverlap <- function(alt.overlap, posteriorProb) {
  ## Simply grabs all genes linked to any variant above the given posterior probability, regardless of cell type
  ## In future, might want to create ranks out of the scores, and consider cell types if provided
  tmp <- subset(alt.overlap, PosteriorProb >= posteriorProb)
  if (nrow(tmp) == 0) return(list())
  results <- with(tmp, tapply(as.character(as.matrix(TargetGene)), Method, unique))
  for (i in 1:length(results)) if (is.null(results[[i]])) results[[i]] <- vector()
  #browser()
  #names(results) <- paste0("GeneList.Prediction.",unique(tmp$Method))
  return(results)
}


addE2GMethodsToGP <- function(gp, alt.overlap, posteriorProb) {
  tmp <- subset(alt.overlap, PosteriorProb >= posteriorProb) %>% mutate(Method=as.character(as.matrix(Method)))
  gp$origOrder <- 1:nrow(gp)
  for (currMethod in unique(as.character(as.matrix(tmp$Method)))) {
    print(currMethod)
    tomerge <- tmp %>% filter(Method==currMethod) %>% select(TargetGene, score) %>% 
      group_by(TargetGene) %>% summarise(tmpScore=max(as.numeric(score))) %>% as.data.frame()
    if (!all(is.na(tomerge$tmpScore))) {
      gp <- merge(gp, tomerge, all.x=TRUE)
      gp <- gp %>% group_by(CredibleSet) %>% mutate(tmpMax=(tmpScore==max(tmpScore, na.rm=T))) %>% as.data.frame()
      gp$tmpMax[is.na(gp$tmpScore)] <- FALSE
      gp[,paste0("GeneList.Prediction.",currMethod,".Max")] <- gp$tmpMax
      gp[,paste0(currMethod,".Score")] <- gp$tmpScore
      gp <- gp %>% select(-tmpScore, -tmpMax)
    }
    
    ## Temporary fix - these should get added in the main getGenePriori... function
    col <- paste0("GeneList.Prediction.",currMethod)
    if (!(col %in% colnames(gp))) gp[,col] <- gp$TargetGene %in% tomerge$TargetGene
  }
  gp <- gp %>% arrange(origOrder) %>% select(-origOrder)
  return(gp)
}

addGeneScoresToGP <- function(gp, gene.scores) {
  scoreCols <- colnames(gene.scores)[grepl("GeneScore",colnames(gene.scores))]
  gp$origOrder <- 1:nrow(gp)
  for (currMethod in scoreCols) {
    print(currMethod)
    tomerge <- data.frame(TargetGene=gene.scores$TargetGene, score=gene.scores[,currMethod]) %>% filter(!is.na(TargetGene)) %>% group_by(TargetGene) %>% summarise(tmpScore=max(score)) %>% as.data.frame()
    if (!all(is.na(tomerge$tmpScore))) {
      gp <- merge(gp, tomerge, all.x=TRUE)
      gp <- gp %>% group_by(CredibleSet) %>% mutate(tmpMax=(tmpScore==max(tmpScore, na.rm=T))) %>% as.data.frame()
      gp$tmpMax[is.na(gp$tmpScore)] <- FALSE
      gp[,paste0(currMethod,".Max")] <- gp$tmpMax
      gp[,paste0(currMethod,".Score")] <- gp$tmpScore
      gp <- gp %>% select(-tmpScore, -tmpMax)
    }
    
    ## Temporary fix - these should get added in the main getGenePriori... function
    col <- currMethod
    gp[,col] <- gp$TargetGene %in% subset(tomerge, !is.na(tmpScore))$TargetGene
  }
  gp <- gp %>% arrange(origOrder) %>% select(-origOrder)
  return(gp)
}
checkInteger <- function(n) n == round(n)

getPrecisionRecall <- function(pred, labels, weights=NULL, ...) {
  ## Used currently by CollateABCTrainingData.R
  stopifnot(is.logical(labels))
  stopifnot(length(pred) == length(labels))
  
  if (is.null(weights)) weights <- rep(1, length(pred))
  
  if (is.logical(pred)) {
    recall <- sum(weights[pred & labels], na.rm=T)/sum(weights[labels], na.rm=T)
    precision <- sum(weights[pred & labels], na.rm=T)/sum(weights[pred], na.rm=T)
    result <- data.frame(
      Precision=precision,
      Recall=recall,
      nPred=sum(weights[pred],na.rm=T),
      nTruePred=sum(weights[pred & labels],na.rm=T),
      nTrue=sum(weights[labels],na.rm=T),
      alpha=TRUE,
      ...
    )
  } else if (all(checkInteger(pred), na.rm=T)) {
    ## Assume this is a ranking and that lower is better
    result <- do.call(rbind, lapply(na.omit(sort(unique(pred))), function(cutoff) {
      tmp <- getPrecisionRecall(pred <= cutoff, labels, weights, ...)
      tmp$alpha <- cutoff
      return(tmp)
    }))
  } else if (is.numeric(pred)) {
    result <- do.call(rbind, lapply(sort(unique(pred)), function(cutoff) {
      tmp <- getPrecisionRecall(pred >= cutoff, labels, weights, ...)
      tmp$alpha <- cutoff
      return(tmp)
    }))
  } else {
    stop("Class of 'pred' not supported: ", class(pred))
  }
  return(result)
}



























#################################################################
################################################################# 
## OLD
#################################################################
#################################################################


if (FALSE) {
  annotateVariants <- function(df, variant.gr, promoters.gr) {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(VariantAnnotation)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    
    library(org.Hs.eg.db)
    
    codingVariants <- locateVariants(variant.gr, txdb, CodingVariants())
    #symbols <- select(org.Hs.eg.db, keys=unique(na.omit(codingVariants$GENEID)), keytype="ENTREZID",
    #                 columns="SYMBOL")
    #df$Coding <- NA
    #df$Coding[codingVariants$QUERYID] <- sapply(codingVariants$GENEID, function(id) {
    #  index <- which(symbols$ENTREZID == id)
    #  if (length(index) == 1) return(symbols$SYMBOL[index])
    #  else return(NA)
    #})
    
    df$Coding <- 1:nrow(df) %in% codingVariants$QUERYID  ## Note that this only means it overlaps CDS, not that it's missense
    
    spliceVariants <- locateVariants(variant.gr, txdb, SpliceSiteVariants())
    #symbols <- select(org.Hs.eg.db, keys=unique(na.omit(spliceVariants$GENEID)), keytype="ENTREZID",
    #                 columns="SYMBOL")
    #df$SpliceSite <- NA
    #df$SpliceSite[spliceVariants$QUERYID] <- sapply(spliceVariants$GENEID, function(id) {
    #  index <- which(symbols$ENTREZID == id)
    #  if (length(index) == 1) return(symbols$SYMBOL[index])
    #  else return(NA)
    #})
    
    df$SpliceSite <- 1:nrow(df) %in% spliceVariants$QUERYID
    ixn <- findOverlaps(variant.gr, promoters.gr)
    df$Promoter <- 1:nrow(df) %in% queryHits(ixn)
    return(df)
  }
  
  annotateVariantsWithPartition <- function(variant.list, partition.file="/seq/lincRNA/RAP/GWAS/170704_Partitioning/170819-All/PartitionCombined.bed", tmp.file=".tmp.partition", sizes="/seq/lincRNA/data/hg19/sizes") {
    writeBed(with(variant.list, data.frame(chr=chr, start=position-1, end=position+1)), file=tmp.file)
    system(paste0("source /broad/software/scripts/useuse; use BEDTools; bedtools intersect -wa -b ", tmp.file, " -a ", partition.file, " > ", tmp.file, ".overlap"))
    partition <- unique(readBed(paste0(tmp.file,".overlap")))
    system(paste0("rm ",tmp.file, " ", tmp.file, ".overlap"))
    return(annotateVariantsWithPartitionHelper(variant.list, partition))
  }
  
  annotateVariantsWithPartitionHelper <- function(variant.list, annotation.partition) {
    ixn <- findOverlaps(GRangesFromBed(with(variant.list, data.frame(chr=chr, start=position-1, end=position))), GRangesFromBed(annotation.partition))
    stopifnot(length(ixn) == nrow(variant.list))
    stopifnot(!any(duplicated(queryHits(ixn))))
    variant.list$Partition <- annotation.partition$name[subjectHits(ixn)]
    return(variant.list)
  }
  
  getPermutedVariantsWithPartition <- function(n.perm, variant.list, annotation.partition, dir="permutations/") {
    do.call(rbind, lapply(1:n.perm, function(i) {
      p <- readBed(paste0(dir,"/",i,"/permuted.variants.bed"))
      p <- merge(p, variant.list[,c("variant","PosteriorProb")], by.x='name', by.y='variant')
      p$position <- p$start
      p <- annotateVariantsWithPartition(p, annotation.partition)
      p
    }))
  }
  
  
  collectGeneExpression <- function(neighborhood.dir, cell.types, col="pA.RNA.RPKM") {
    result <- lapply(cell.types, function(ct) {
      file <- paste0(neighborhood.dir,"/",ct,"/GeneList.txt")
      print(file)
      if (file.exists(file)) {
        g <- read.delim(file)
        expr <- g[,col]
        names(expr) <- g[,"name"]
        return(expr)
      } else {
        return(NA)
      }
    })
    num.genes <- max(as.numeric(names(table(unlist(lapply(result, length))))))
    for (i in 1:length(result)) 
      if (length(result[[i]]) == 1)
        result[[i]] <- rep(NA, num.genes)
    result <- do.call(cbind, result)
    colnames(result) <- cell.types
    return(result)
  }
  
  
  makeCredibleSets <- function(df, genes, Source=NA) {
    if ("Partition" %in% colnames(df)) {
      res <- do.call(rbind, tapply(1:nrow(df), df$CredibleSet, function(i) {
        with(df[i,], data.frame(chr=as.character(as.matrix(chr[1])), start=min(position), end=max(position),CredibleSet=CredibleSet[1], nSNP=length(position), AnyCoding=any(Partition == "CDS"), AnyPromoter=any(Partition == "TSS-100bp"), AnySpliceSite=any(Partition == "SpliceSite"), Disease=Disease[1], Source=Source, BestSNP=variant[which.max(PosteriorProb)[1]]))   
      }))
    } else {
      res <- do.call(rbind, tapply(1:nrow(df), df$CredibleSet, function(i) {
        with(df[i,], data.frame(chr=as.character(as.matrix(chr[1])), start=min(position), end=max(position),CredibleSet=CredibleSet[1], nSNP=length(position), AnyCoding=any(Coding), AnyPromoter=any(Promoter), AnySpliceSite=any(SpliceSite), Disease=Disease[1], Source=Source, LocusID=LocusID[1], BestSNP=variant[which.max(PosteriorProb)[1]]))   
      }))
    }
    res <- merge(res, unique(data.frame(BestSNP=df$variant, BestSNPPos=df$position)), all.x=TRUE)
    res <- res[,c(2:ncol(res),1)]
    
    coding.genes <- subset(genes, !grepl("NR_", name))
    ixn <- as.data.frame(nearest(GRangesFromBed(with(res, data.frame(chr=chr, start=BestSNPPos, end=BestSNPPos))),
                                 GRangesFromBed(coding.genes), select="all", ignore.strand=TRUE))
    res$BestSNPNearestGene <- sapply(1:nrow(res), function(i) paste(unique(coding.genes$symbol[subset(ixn, queryHits == i)$subjectHits]), collapse=","))
    return(res)
  }
  
  
  writeCredibleSetBed <- function(all.cs, variant.list, file) {
    ## Write credible set for viewing ... with variants as "exons"
    bed <- NULL
    for (i in 1:nrow(all.cs)) {
      nSNP <- all.cs$nSNP[i]
      dz <- all.cs$Disease[i]
      cs <- all.cs$CredibleSet[i]
      variants <- subset(variant.list, CredibleSet == cs)
      variants <- variants[order(variants$position),]
      bed.entry <- with(variants, data.frame(
        chr=chr[1], start=min(position)-1, end=max(position), name=paste(dz,cs,sep="_"),
        score=nSNP, strand=".", thickStart=min(position)-1, thickEnd=max(position), itemRgb="0,0,0", 
        blockCount=nSNP, blockSizes=paste0(rep("1,",nSNP), collapse=""), 
        blockStarts=paste0(sort(position)-all.cs$start[i],collapse=",")))
      if (nrow(bed.entry) != 1) browser()
      if (is.null(bed)) bed <- unfactor(bed.entry)
      else bed <- rbind(bed, bed.entry)
    }
    writeBed(bed, file)
    return(bed)
  }
  
  writeVariantBed <- function(variant.list, file, name.only=FALSE) {
    if (name.only) 
      tmp <- with(variant.list, data.frame(chr=chr, start=position-1, end=position, name=variant))
    else
      tmp <- with(variant.list, data.frame(chr=chr, start=position-1, end=position, name=paste0(variant,"-",Disease), score=PosteriorProb, strand="+"))
    tmp <- tmp[order(tmp$chr, tmp$start),]
    writeBed(format(tmp, scientific=F), file)
    return(tmp)
  }
  
  
  plotCredibleSetProperties <- function(all.cs, filter.cs=NULL) {
    hist(table(all.cs$LocusID), col='gray', xlab="Number of Diseases per Locus", ylab="Count", main="All Credible Sets: locus = within 500Kb")
    hist(all.cs$nSNP, col='gray', breaks=100, xlim=c(0,50), xlab="Number of SNPs in Credible Set", main="All Credible Sets")
    hist(subset(all.cs, !AnyCoding & !AnySpliceSite)$nSNP, col='gray', breaks=100, xlim=c(0,50), xlab="Number of SNPs in Credible Set", main="Credible Sets w/o Coding or Splice Site Variants")
    plotCumulativeDistribution(subset(all.cs, AnyCoding | AnySpliceSite)$nSNP, 
                               subset(all.cs, !AnyCoding & !AnySpliceSite)$nSNP, "Coding or Splice Site", "Neither",
                               xlim=c(0,50), ylab="Cumulative fraction", xlab="# SNPs in Credible Set")
    plotCumulativeDistribution(subset(all.cs, Source == "Farh2015" & !AnyCoding & !AnySpliceSite)$nSNP, 
                               subset(all.cs, Source != "Farh2015" & !AnyCoding & !AnySpliceSite)$nSNP, "Farh2015", "Other",
                               xlim=c(0,50), ylab="Cumulative fraction", xlab="# SNPs in Credible Set", main="Credible Sets w/o Coding or Splice site variants")
    plot(ecdf(subset(all.cs, !AnyCoding & !AnySpliceSite)$nSNP), xlim=c(0,50), ylab="Cumulative frequency", xlab="# SNPs in Credible Set")
    
    if (!is.null(filter.cs)) {
      hist(filter.cs$nSNP, col='gray', breaks=100, xlim=c(0,50), xlab="Number of SNPs in Credible Set", main="Filtered Credible Sets")
      plotCumulativeDistribution(subset(filter.cs, Source == "Farh2015")$nSNP, 
                                 subset(filter.cs, Source != "Farh2015")$nSNP, "Farh2015", "Other",
                                 xlim=c(0,50), ylab="Cumulative fraction", xlab="# SNPs in Credible Set", main="Filtered Credible Sets w/o Coding or ss variants")
      plot(ecdf(filter.cs$nSNP), xlim=c(0,50), ylab="Cumulative frequency", xlab="# SNPs in Credible Set", main="Filtered Credible Sets")
    }
  }
  
  
  getBestSNP <- function(df) {
    tmp <- unlist(by(df, df$CredibleSet, function(x) as.character(as.matrix(x$variant[which.max(x$pvalue)])), simplify=TRUE))
    tmp[as.character(as.matrix(df$CredibleSet))]    
  }
  
  getCredibleSetLoci <- function(all.cs, buffer=500000, overlap.sets=NULL) {
    ## e.g. overlap.sets = list(c("Crohns_disease","IBD","Ulcerative_colitis"))
    ## Find those with associations with more than one disease
    
    all.cs.gr <- GRangesFromBed(all.cs[,c("chr","start","end")])
    ixn <- findOverlaps(all.cs.gr, all.cs.gr, maxgap=500000)
    to.keep <- unique(unlist(mapply(function(i,j) {
      di <- all.cs$Disease[i]
      dj <- all.cs$Disease[j]
      
      if ( (i != j) & (di != dj) & !any(unlist(lapply(overlap.sets, function(overlap.set) all(c(di,dj) %in% overlap.set)))) )
        return(c(i,j))
    }, queryHits(ixn), subjectHits(ixn))))
    
    loci <- reduce(all.cs.gr, min.gapwidth=500000, drop.empty.ranges=TRUE, ignore.strand=TRUE) 
    return(loci)
  }
  
  assignLocusId <- function(all.cs, loci) {
    ixn <- findOverlaps(GRangesFromBed(all.cs[,c("chr","start","end")]), loci)
    all.cs$LocusID <- subjectHits(ixn)
    return(all.cs)
  }
  
  
  
  filterPredictions <- function(pred.flat, filter.cs, tss.cutoff=0.1, require.gex=FALSE) {
    ## Implements a different cutoff for distal enhancers versus distal promoters
    all.flat <- subset(pred.flat, ( (!is.na(GeneExpression) | !require.gex) & ABC.Score >= opt$cutoffLenient & ((class != "tss" & class != "promoter") | isOwnTSS)) | (Score.Fraction >= tss.cutoff & (class == "tss" | class == "promoter")))
    filter.flat <- subset(all.flat, CredibleSet %in% as.character(as.matrix(filter.cs$CredibleSet)))
    return(list(all=all.flat, filter=filter.flat))
  }
  
  
  getPredictionsForVariant <- function(variant, pred) {
    do.call(rbind, lapply(pred, function(cell.pred) {
      curr <- subset(cell.pred, QueryRegionName == variant)
      if (nrow(curr) > 0)
        curr[,common.col]
    }))
  }
  
  annotateCredibleSet <- function(name, variants, pred, fraction.cutoff=opt$cutoff, abs.cutoff=opt$absCutoff) {
    variant.overlaps <- do.call(rbind, lapply(as.character(as.matrix(variants)), function(variant) getPredictionsForVariant(variant, pred)))
    if (!is.null(variant.overlaps)) { 
      variant.overlaps$CredibleSet <- name
      result <- subset(variant.overlaps, (ABC.Score >= fraction.cutoff) | (Score >= abs.cutoff))
      if (nrow(result) == 0) return(NULL)
      return(result)
    } 
    return(NULL)
  }
  
  
  annotateCredibleSets <- function(variants, cs, pred) {
    names <- as.character(as.matrix(cs$CredibleSet))
    cs.annot <- lapply(names, function(cs.name) {
      print(subset(variants, CredibleSet == cs.name)$variant)
      annotateCredibleSet(cs.name, subset(variants, CredibleSet == cs.name)$variant, pred)
    })
    names(cs.annot) <- names
    cs.annot <- cs.annot[!unlist(lapply(cs.annot, function(curr) is.null(curr) | nrow(curr) == 0))]
    return(cs.annot)
  }
  
  
  annotateCredibleSetFlat <- function(variants, cs, pred) {
    names <- as.character(as.matrix(cs$CredibleSet))
    cs.annot <- lapply(names, function(cs.name) {
      print(subset(variants, CredibleSet == cs.name)$variant)
      annotateCredibleSet(cs.name, subset(variants, CredibleSet == cs.name)$variant, pred)
    })
    names(cs.annot) <- names
    cs.annot <- cs.annot[!unlist(lapply(cs.annot, function(curr) is.null(curr) | nrow(curr) == 0))]
    return(cs.annot)
  }
  
  
  
  getCredibleSetVariants <- function(cs, fraction.cutoff=NULL, abs.cutoff=NULL) {
    cs$QueryRegionName <- as.character(as.matrix(cs$QueryRegionName))
    if (length(cs) > 0) {
      if (!is.null(fraction.cutoff)) cs <- subset(cs, ABC.Score >= fraction.cutoff)
      if (!is.null(abs.cutoff)) cs <- subset(cs, Score >= abs.cutoff)
    }
    if (length(cs) > 0) table(cs[,c("QueryRegionName","CellType")])
  }
  
  getCredibleSetGenes <- function(cs, fraction.cutoff=NULL, abs.cutoff=NULL) {
    ## Number of enhancers regulating each gene
    if (length(cs) > 0) {
      if (!is.null(fraction.cutoff)) cs <- subset(cs, ABC.Score >= fraction.cutoff)
      if (!is.null(abs.cutoff)) cs <- subset(cs, Score >= abs.cutoff)
    }
    if (length(cs) > 0) with(cs[!duplicated(cs[,c("TargetGene","name","CellType")]),], table(as.matrix(TargetGene), CellType))
  }
  
  getCredibleSetElements <- function(cs, fraction.cutoff=NULL, abs.cutoff=NULL) {
    ## Tally of affected elements per cell type
    if (length(cs) > 0) {
      if (!is.null(fraction.cutoff)) cs <- subset(cs, ABC.Score >= fraction.cutoff)
      if (!is.null(abs.cutoff)) cs <- subset(cs, Score >= abs.cutoff)
    }
    if (length(cs) > 0) with(cs[!duplicated(cs[,c("name","CellType")]),], table(as.matrix(name), CellType))
  }
  
  getCredibleSetElementsPerCellType <- function(cs, fraction.cutoff=NULL, abs.cutoff=NULL) {
    ## Number of affected elements per cell type
    if (length(cs) > 0) {
      if (!is.null(fraction.cutoff)) cs <- subset(cs, ABC.Score >= fraction.cutoff)
      if (!is.null(abs.cutoff)) cs <- subset(cs, Score >= abs.cutoff)
    }
    if (length(cs) > 0) apply(getCredibleSetElements(cs), 2, sum)
  }
  
  getCredibleSetCellTypes <- function(cs, fraction.cutoff=NULL, abs.cutoff=NULL) {
    #if (length(cs) > 0) 
    if (length(cs) > 0) {
      if (!is.null(fraction.cutoff)) cs <- subset(cs, ABC.Score >= fraction.cutoff)
      if (!is.null(abs.cutoff)) cs <- subset(cs, Score >= abs.cutoff)
    }
    if (length(cs) > 0) unique(cs$CellType)
  }
  
  
  addUniqueCountsPerCS <- function(all.flat, all.cs, col) {
    count.col <- paste0('n',col,'s')
    cells.per.cs <- with(all.flat, tapply(get(col), CredibleSet, function(z) length(unique(z))))
    tmp <- data.frame(CredibleSet = factor(names(cells.per.cs), levels=levels(all.cs$CredibleSet)))
    tmp[,count.col] <- cells.per.cs
    all.cs <- merge(all.cs, tmp, all.x=T)
    all.cs[,count.col][is.na(all.cs[,count.col])] <- 0
    return(all.cs)  
  }
  
  getUniqueCountsFromOverlap <- function(all.flat, ids, count.feature, by.feature) {
    ids <- as.character(as.matrix(ids))
    count.per.ct <- with(all.flat, tapply(get(count.feature), get(by.feature), function(z) length(unique(z))))
    result <- rep(0, length(ids))
    names(result) <- ids
    for (id in ids)
      if (id %in% names(count.per.ct))
        result[names(result) == id] <- count.per.ct[id]
    result[is.na(result)] <- 0
    return(result)
  }
  
  getCredibleSetCounts <- function(cs.plot, all.flat, control.types) {
    ## Number of cell types per association, variant or locus
    cs.plot$hasExpressedPromoterVariant <- cs.plot$CredibleSet %in% unique(subset(all.flat, isOwnTSS)$CredibleSet)
    cs.plot$nCellTypes <- getUniqueCountsFromOverlap(all.flat, cs.plot$CredibleSet, "CellType", "CredibleSet")
    cs.plot$nControlCellTypes <- getUniqueCountsFromOverlap(subset(all.flat, CellType %in% control.types), cs.plot$CredibleSet, "CellType", "CredibleSet")
    cs.plot$nImmuneCellTypes <- getUniqueCountsFromOverlap(subset(all.flat, !(CellType %in% control.types)), cs.plot$CredibleSet, "CellType", "CredibleSet")
    cs.plot$nTargetGenes <- getUniqueCountsFromOverlap(all.flat, cs.plot$CredibleSet, "TargetGene", "CredibleSet")
    return(cs.plot)
  }
  
  getVariantCounts <- function(all.flat, variants, control.types) {
    variants$nCellTypes <- getUniqueCountsFromOverlap(all.flat, variants$variant, "CellType", "QueryRegionName")
    variants$nControlCellTypes <- getUniqueCountsFromOverlap(subset(all.flat, CellType %in% control.types), variants$variant, "CellType", "QueryRegionName")
    variants$nImmuneCellTypes <- getUniqueCountsFromOverlap(subset(all.flat, !(CellType %in% control.types)), variants$variant, "CellType", "QueryRegionName")
    variants$nTargetGenes <- getUniqueCountsFromOverlap(all.flat, variants$variant, "TargetGene", "QueryRegionName")
    return(variants)
  }
  
  getCellTypeCounts <- function(all.cell.types, all.flat, control.types, variants, cs.plot) {
    cell.types <- data.frame(CellType=all.cell.types, IsControl=all.cell.types %in% control.types)
    cell.types$nCredibleSets <- getUniqueCountsFromOverlap(all.flat, cell.types$CellType, "CredibleSet", "CellType")
    ## do it in subset where we only look at variants that don't overlap many cell types
    cell.types$nCredibleSetsMinusUbiq <- getUniqueCountsFromOverlap(subset(all.flat, QueryRegionName %in% subset(variants, nCellTypes <= 17)$variant), cell.types$CellType, "CredibleSet", "CellType")
    cell.types$nCredibleSetsMinusControl <- getUniqueCountsFromOverlap(subset(all.flat, QueryRegionName %in% subset(variants, nControlCellTypes == 0)$variant), cell.types$CellType, "CredibleSet", "CellType")
    cell.types$nCredibleSetsMinusK562 <- getUniqueCountsFromOverlap(subset(all.flat, CellType != "K562"), cell.types$CellType, "CredibleSet", "CellType")
    cell.types$nSmallCredibleSets <- getUniqueCountsFromOverlap(subset(all.flat, CredibleSet %in% subset(cs.plot, nSNP <= 10)$CredibleSet), cell.types$CellType, "CredibleSet", "CellType")
    cell.types <- cell.types[order(cell.types$IsControl),]
    return(cell.types)
  }
  
  
  plotCellTypeCounts <- function(cell.types, ...) {
    par(mar=c(5,20,5,3))
    p <- cell.types$nCredibleSets; names(p) <- cell.types$CellType
    colors <- c('red','gray')[as.numeric(cell.types$IsControl[order(p)])+1]
    p <- p[order(p)]
    par(mar=c(5.1,6.1,4.1,2.1))
    b <- barplot(p, horiz=TRUE, yaxt='n', col=colors, xlab="# Credible Sets", ...)
    axis(2, b, labels=names(p), cex=0.5, las=2)
    
    par(mar=c(5,10,5,3))
    p <- cell.types$nCredibleSetsMinusUbiq; names(p) <- cell.types$CellType
    colors <- c('red','gray')[as.numeric(cell.types$IsControl[order(p)])+1]
    p <- p[order(p)]
    par(mar=c(5.1,6.1,4.1,2.1))
    b <- barplot(p, horiz=TRUE, yaxt='n', col=colors, xlab="# Credible Sets (-ubiq variants)", ...)
    axis(2, b, labels=names(p), cex=0.5, las=2)
    
    par(mar=c(5,10,5,3))
    p <- cell.types$nSmallCredibleSets; names(p) <- cell.types$CellType
    colors <- c('red','gray')[as.numeric(cell.types$IsControl[order(p)])+1]
    p <- p[order(p)]
    par(mar=c(5.1,6.1,4.1,2.1))
    b <- barplot(p, horiz=TRUE, yaxt='n', col=colors, xlab="# Credible sets (<=10 variants)", ...)
    axis(2, b, labels=names(p), cex=0.5, las=2)
  }
  
  
  plotEnhancerOverlapStats <- function(cs.plot, cell.types, variants, all.flat, ...) {
    ## Show CDF for cell types per association, versus permuted, broken down by immune cell types vs control cell types
    b <- barplot(getSummaryCSCellTypeCounts(cs.plot), main="# Credible Sets with EP Predictions\nin cell type category", xlab="All CS", horiz=TRUE, yaxt='n')
    axis(2, b, labels=c("All","Immune","Control","Both","ControlOnly"), las=2, cex=0.7)
    b <- barplot(getSummaryCSCellTypeCounts(subset(cs.plot, nSNP <= 10)), main="# Credible Sets with EP Predictions\nin cell type category", xlab="CS with <= 10 variants", horiz=TRUE, yaxt='n')
    axis(2, b, labels=c("All","Immune","Control","Both","ControlOnly"), las=2, cex=0.7)
    
    ## Number of cell types per association, variant or locus
    hist(cs.plot$nCellTypes, col='gray', breaks=0:max(cs.plot$nCellTypes), main="Credible Sets with EP Predictions", xlab="Cell types per credible set", cex.lab=0.25)
    for (disease in unique(cs.plot$Disease)) {
      hist(subset(cs.plot, Disease == disease)$nCellTypes, col='gray', breaks=0:max(cs.plot$nCellTypes), main=paste0("Credible Sets with EP Predictions\n(",disease,")"), xlab="Cell types per credible set")
    }
    #hist(subset(cs.plot, CredibleSet %in% filter.cs$CredibleSet)$nCellTypes, col='gray', breaks=0:max(cs.plot$nCellTypes), main="Filtered Credible Sets", xlab="Cell types per credible set")
    if (!is.null(variants)) 
      hist(variants$nCellTypes, col='gray', breaks=0:max(variants$nCellTypes,na.rm=T), main="CS Variants with EP Predictions", xlab="Cell types per variant")
    
    ## Number of genes per association, variant, or locus
    hist(cs.plot$nTargetGenes, col='gray', breaks=0:max(cs.plot$nTargetGenes), main="Credible Sets with EP Predictions", xlab="Genes per credible set")
    #hist(subset(cs.plot, CredibleSet %in% filter.cs$CredibleSet)$nTargetGenes, col='gray', breaks=0:max(cs.plot$nTargetGenes), main="Filtered Credible Sets", xlab="Genes per credible set")
    if (!is.null(variants)) 
      hist(variants$nTargetGenes, col='gray', breaks=0:max(variants$nTargetGenes,na.rm=T), main="Variants with EP Predictions", xlab="Genes per variant")
    
    ## Relationship between # of diseases and # of immune cell types per association / variant
    
    ## Number of overlaps with immune cells vs number of overlaps with matched other cell type
    ## (enrichment in immune cell overlaps)
    plotCellTypeCounts(cell.types, main="CS per cell type (all dz)")
    for (disease in unique(cs.plot$Disease)) {
      dz.cell.types <- getCellTypeCounts(as.character(as.matrix(cell.types$CellType)), subset(all.flat, CredibleSet %in% subset(cs.plot, Disease == disease)$CredibleSet), control.types, variants, cs.plot)
      plotCellTypeCounts(dz.cell.types, main=paste0("CS per cell type (",disease,")"))
    }
    
    
    ## Pathway enrichment of predicted genes vs closest genes
    
    ## Highlight examples where closest gene differs from predictions
    
    ## Properties of G-E connections:
    ##   Genes per enhancer as function of cell type
    ##   Variants that regulate different genes across cell types?
    ##   This could be a set of functions that we use for analyzing genome-wide predictions as well
  }
  
  
  
  getSummaryCSCellTypeCounts <- function(cs.plot) {
    cs.plot$nBothCellTypes <- with(cs.plot, as.numeric(nImmuneCellTypes > 0 & nControlCellTypes > 0))
    cs.plot$nControlOnly <- with(cs.plot, as.numeric(nImmuneCellTypes == 0 & nControlCellTypes > 0))
    apply(cs.plot[,c("start","nImmuneCellTypes","nControlCellTypes","nBothCellTypes","nControlOnly")], 2, function(z) sum(z > 0))
  }
  
  
  plotCellTypeEnrichmentBarplots <- function(cs.plot, cell.types, permute.cs, permute.cell.types, n.perm, ...) {
    p <- cell.types$nCredibleSets / (permute.cell.types$nCredibleSets/n.perm + 0.1); names(p) <- cell.types$CellType
    colors <- c('red','gray')[as.numeric(cell.types$IsControl[order(p)])+1]
    p <- p[order(p)]
    par(mar=c(5.1,8.1,4.1,2.1))
    b <- barplot(p, horiz=TRUE, yaxt='n', col=colors, xlab="Enrichment (Observed / Expected)", main="Credible Sets per Cell Type", ...)
    axis(2, b, labels=names(p), cex.axis=0.8, las=2)
    abline(v=1, lty=2, col='gray')
    
    p <- cell.types$nSmallCredibleSets / (permute.cell.types$nSmallCredibleSets/n.perm + 0.1); names(p) <- cell.types$CellType
    colors <- c('red','gray')[as.numeric(cell.types$IsControl[order(p)])+1]
    p <- p[order(p)]
    par(mar=c(5.1,8.1,4.1,2.1))
    b <- barplot(p, horiz=TRUE, yaxt='n', col=colors, xlab="Enrichment (Observed / Expected)", main="Credible Sets (n<=10) per Cell Type", ...)
    axis(2, b, labels=names(p), cex.axis=0.8, las=2)
    abline(v=1, lty=2, col='gray')
    
    
    p <- -log10(mapply(function(real, perm) pbinom(real, nrow(cs.plot), perm/nrow(permute.cs), lower.tail=F), cell.types$nCredibleSets, permute.cell.types$nCredibleSets))
    names(p) <- cell.types$CellType
    colors <- c('red','gray')[as.numeric(cell.types$IsControl[order(p)])+1]
    p <- p[order(p)]
    par(mar=c(5.1,8.1,4.1,2.1))
    b <- barplot(p, horiz=TRUE, yaxt='n', col=colors, xlab="-log10 P-val (Observed / Expected)", main="Credible Sets per Cell Type", ...)
    axis(2, b, labels=names(p), cex.axis=0.8, las=2)
    abline(v=3, lty=2, col='gray')
    
    tryCatch({
      p <- -log10(mapply(function(real, perm) pbinom(real, with(cs.plot, sum(nSNP<=10)), perm/with(permute.cs, sum(nSNP<=10)), lower.tail=F), cell.types$nSmallCredibleSets, permute.cell.types$nSmallCredibleSets))
      names(p) <- cell.types$CellType
      colors <- c('red','gray')[as.numeric(cell.types$IsControl[order(p)])+1]
      p <- p[order(p)]
      par(mar=c(5.1,8.1,4.1,2.1))
      b <- barplot(p, horiz=TRUE, yaxt='n', col=colors, xlab="-log10 P-val (Observed / Expected)", main="Credible Sets (n<=10) per Cell Type", ...)
      axis(2, b, labels=names(p), cex.axis=0.8, las=2)
      abline(v=3, lty=2, col='gray')
    }, error = function(e) print("Skipping")
    )
  }
  
  
  plotEnhancerOverlapEnrichment <- function(cs.plot, cell.types, permute.cs, permute.cell.types, n.perm=10, ...) {
    ## Show CDF for cell types per association, versus permuted, broken down by immune cell types vs control cell types
    par(mar=c(5.1,8.1,4.1,4.1))
    tp <- rbind(getSummaryCSCellTypeCounts(cs.plot), getSummaryCSCellTypeCounts(permute.cs) / n.perm)
    b <- barplot(tp, main="# Credible Sets with EP Predictions\nin cell type category", xlab="All CS", horiz=T, yaxt='n', beside=T)
    text(x = tp, y = b, label = format(tp,digits=1), pos = 4, cex = 1, col = "black")
    axis(2, selectEveryNth(b,2), labels=c("All","Immune","Control","Both","ControlOnly"), las=2, cex=0.7)
    legend("topright", c("Expected (permuted)","Observed"), col=c("gray","black"), pch=19)
    tp <- rbind(getSummaryCSCellTypeCounts(subset(cs.plot, nSNP <= 10)), getSummaryCSCellTypeCounts(subset(permute.cs, nSNP <= 10)) / n.perm)
    b <- barplot(tp, main="# Credible Sets with EP Predictions\nin cell type category", xlab="CS with <= 10 variants", horiz=TRUE, yaxt='n', beside=T)
    axis(2, selectEveryNth(b,2), labels=c("All","Immune","Control","Both","ControlOnly"), las=2, cex=0.7)
    legend("topright", c("Expected (permuted)","Observed"), col=c("gray","black"), pch=19)
    text(x = tp, y = b, label = format(tp,digits=1), pos = 4, cex = 1, col = "black")
    
    plotCellTypeEnrichmentBarplots(cs.plot, cell.types, permute.cs, permute.cell.types, n.perm, ...)
  }
  
  
  plotGeneByCellTypeHeatmap <- function(cs.name, filter.flat, variant.list, genes, cell.list, buffer=1000000, write.matrix=NULL, ...) {
    cs.variants <- subset(variant.list, CredibleSet == cs.name)
    curr.pred <- subset(filter.flat, CredibleSet == cs.name)
    if (nrow(curr.pred) == 0) return()
    curr.pred$CellType <- as.character(as.matrix(curr.pred$CellType))
    
    ## Plot for all genes within 1 Mb of variants or out to furthest gene
    #plot.genes <- unique(subset(genes, chr==as.character(as.matrix(cs.variants$chr[1])) & 
    #  end >= min(min(cs.variants$position)-buffer,min(curr.pred$TargetGeneTSS)) & 
    #  start <= max(max(cs.variants$position)+buffer,max(curr.pred$TargetGeneTSS)))$symbol)
    plot.genes <- subset(genes, symbol %in% as.character(as.matrix(unique(curr.pred$TargetGene))))
    plot.genes <- plot.genes[!duplicated(plot.genes$symbol),]
    plot.genes <- plot.genes[order(plot.genes$start),]
    plot.genes <- plot.genes$symbol
    
    tab <- matrix(0, nrow=length(cell.list)+1, ncol=length(plot.genes)+1)
    rownames(tab) <- c(cell.list,""); colnames(tab) <- c(plot.genes,"")
    for (i in 1:nrow(curr.pred)) {
      target.gene <- as.character(as.matrix(curr.pred[i,"TargetGene"]))
      if (!(target.gene %in% plot.genes)) {
        cat(paste0("Found gene that doesn't match: ", target.gene, "\n"))
      } else {
        tryCatch({
          tab[curr.pred[i,"CellType"], target.gene] <- curr.pred[i,"ABC.Score"]*100
        }, error = function(e) { print("Caught in plotGeneByCellTypeHeatmap"); browser() })
      }
    }
    
    if (ncol(tab) < 2) {
      print(paste0("Skipping ",cs.name," because only 1 gene is predicted"))
    } else {
      plotCellTypeHeatmap(tab, ...)
      if (!is.null(write.matrix)) {
        write.table(tab, file=write.matrix, sep='\t', quote=F)
      }
    }
  }
  
  
  plotVariantByCellTypeHeatmap <- function(cs.name, filter.flat, variant.list, cell.list, buffer=1000000, ...) {
    cs.variants <- subset(variant.list, CredibleSet == cs.name)
    cs.variants <- cs.variants[order(cs.variants$position),]
    curr.pred <- subset(filter.flat, CredibleSet == cs.name)
    
    tab <- matrix(0, nrow=length(cell.list)+1, ncol=nrow(cs.variants)+1)
    rownames(tab) <- c(cell.list,""); colnames(tab) <- c(as.character(as.matrix(cs.variants$variant)),"")
    for (i in 1:nrow(curr.pred)) {
      variant <- as.character(as.matrix(curr.pred[i,"QueryRegionName"]))
      if (!(variant %in% colnames(tab))) {
        print(paste0("Skipping ", variant, " because variant is not found"))
      } else {
        tab[curr.pred[i,"CellType"], variant] <- curr.pred[i,"ABC.Score"]*100
      }
    }
    
    if (ncol(tab) < 2) {
      print(paste0("Skipping ",cs.name," because only 1 variant is predicted"))
    } else {
      plotCellTypeHeatmap(tab, ...)
    }
  }
  
  
  plotVariantByGeneHeatmap <- function(cs.name, filter.flat, variant.list, genes, buffer=1000000, ...) {
    cs.variants <- subset(variant.list, CredibleSet == cs.name)
    cs.variants <- cs.variants[order(cs.variants$position),]
    curr.pred <- subset(filter.flat, CredibleSet == cs.name)
    
    plot.genes <- subset(genes, symbol %in% as.character(as.matrix(unique(curr.pred$TargetGene))))
    plot.genes <- plot.genes[!duplicated(plot.genes$symbol),]
    plot.genes <- plot.genes[order(plot.genes$start),]
    plot.genes <- plot.genes$symbol
    
    tab <- matrix(0, nrow=length(plot.genes)+1, ncol=nrow(cs.variants)+1)
    rownames(tab) <- c(plot.genes,""); colnames(tab) <- c(as.character(as.matrix(cs.variants$variant)),"")
    for (i in 1:nrow(curr.pred)) {
      variant <- as.character(as.matrix(curr.pred[i,"QueryRegionName"]))
      target.gene <- as.character(as.matrix(curr.pred[i,"TargetGene"]))
      if (!(target.gene %in% plot.genes)) {
        print(paste0("Found gene that doesn't match: ", target.gene, "\n"))
      } else if (!(variant %in% colnames(tab))) {
        print(paste0("Skipping ", variant, " because variant is not found"))
      } else {
        tab[target.gene, variant] <- curr.pred[i,"ABC.Score"]*100
      }
    }
    
    if (ncol(tab) < 2 | nrow(tab) < 2) {
      print(paste0("Skipping ",cs.name," because only 1 variant is predicted"))
    } else {
      plotCellTypeHeatmap(tab, ...)
    }
  }
  
  plotCredibleSetCellTypeHeatmap <- function(all.flat, return.raw=FALSE, histogram=TRUE, ...) {
    all.flat <- unfactor(all.flat)
    tab <- with(all.flat, table(CredibleSet, CellType))
    tab[tab > 0] <- 1
    tab <- as.matrix(tab)
    
    # add gene list to credible set names
    tmp <- unique(all.flat[,c("CredibleSet","TargetGene")])
    cs.gene.map <- with(tmp, tapply(TargetGene, CredibleSet, paste0, collapse=','))
    rownames(tab) <- paste0(rownames(tab),"  ",cs.gene.map[rownames(tab)])
    
    par(mar=c(12,0,4,8))
    tab.sorted <- plotCredibleSetCellTypeHeatmapHelper(tab, histogram=histogram, breaks=c(-1, 0.5, 1), col=c('white','red'), ...)
    
    if (return.raw) return(tab)
    else return(tab.sorted)
  }
  
  plotCredibleSetCellTypeHeatmapHelper <- function(tab, histogram=TRUE, margins=c(16,16), ...) {
    #tab <- matrix(0, nrow=length(unique(all.flat$CellType)), ncol=length(unique(all.flat$CredibleSet)))
    ## SUggest providing col (colors) and breaks (vector)
    library(gplots)
    
    res <- heatmap.2(tab, trace='none', key=F, margins=margins, cex.main=0.8, cexCol=1, cexRow=1, ...)
    
    tab <- tab[res$rowInd, res$colInd]
    
    if (histogram) {
      par(mar=c(3,1,1,3))
      h.row <- apply(tab, 1, sum)
      h.col <- apply(tab, 2, sum)
      barplot(h.col, axes=T, ylim=c(0,max(h.col)), space=0, col='gray')
      barplot(h.row, axes=T, xlim=c(0,max(h.row)), space=0, col='gray', horiz=T)
    }
    return(tab)
  }
  
  
  plotCellTypeHeatmap <- function(tab, ...) {
    ## Plot a heatmap with genes or variants on X-axis and cell types on Y-axis
    require(gplots)
    require(RColorBrewer)
    par(mar=c(8.1,4.1,4.1,8.1))
    palette <- colorRampPalette(c("white","red"))(n = 12)
    heatmap.2(tab, Rowv=F, Colv=F, trace='none', dendrogram='none', col=palette, key=F, 
              breaks=c(0:10,30,100), margins=c(6,12), cex.main=0.8, cexCol=1, cexRow=1, ...)
  }
  
  getGenesInLocus <- function(genes, locus, buffer=200000) {
    ixn <- findOverlaps(genes, locus, maxgap=buffer)  
    return(unique(genes[queryHits(ixn)]$name))
  }
  
  getGeneCellTypeMatrix <- function(gene.name, filter.cs, all.flat, loci, cell.lines, gex, genes=NULL, buffer=200000, include=NULL, locus=NULL) {
    require(reshape2)
    if (!is.null(gene.name)) {
      credible.sets <- subset(filter.cs, CredibleSet %in% subset(all.flat, TargetGene == gene.name)$CredibleSet)
      locus <- subjectHits(findOverlaps(GRangesFromBed(credible.sets[,c("chr","start","end")]), loci))[1]
    } else if (is.null(locus)) {
      print("Must specify locus or gene.name")
      return(NULL)
    }
    credible.sets <- subset(filter.cs, LocusID == locus)
    curr.pred <- subset(all.flat, CredibleSet %in% credible.sets$CredibleSet & CellType %in% cell.lines)
    
    if (is.null(genes)) {
      predictions <- dcast(curr.pred[,c("CellType","TargetGene","ABC.Score")], CellType ~ TargetGene, mean, value.var="Score.Fraction")
      rownames(predictions) <- predictions$CellType
      return(list(predictions=t(predictions[,-1]), gex=gex[colnames(predictions),predictions$CellType]))
    } else {
      if (!nrow(curr.pred)) return(NULL)
      nearby.genes <- intersect(unique(c(getGenesInLocus(GRangesFromBed(genes), loci[locus], buffer=buffer), as.character(as.matrix(curr.pred$TargetGene)))), rownames(gex))
      if (!is.null(include))
        nearby.genes <- c(nearby.genes, include)
      predictions <- matrix(NA, nrow=length(nearby.genes), ncol=length(cell.lines))
      rownames(predictions) <- nearby.genes
      colnames(predictions) <- cell.lines
      for (i in 1:nrow(curr.pred)) {
        predictions[as.character(as.matrix(curr.pred$TargetGene[i])), 
                    as.character(as.matrix(curr.pred$CellType[i]))] <- curr.pred$ABC.Score[i]
      }
      #browser()
      return(list(predictions=predictions, gex=gex[rownames(predictions), colnames(predictions)], CredibleSets=credible.sets$CredibleSet))
    }
  }
  
  getGeneAnnotationDf <- function(m, all.flat, cell.lines, genes.uniq) {
    nearby.genes <- rownames(m$predictions)
    
    ## Reformat as data frame and include annotations
    g <- subset(genes.uniq, name %in% nearby.genes)
    g$name <- as.character(as.matrix(g$name))
    g$TargetGene <- g$name
    common <- intersect(g$name, rownames(m$gex))
    for (cell.type in cell.lines) {
      g[,paste0("TPM.",cell.type)] <- m$gex[g$name,cell.type]
    }
    for (cell.type in cell.lines) {
      g[,paste0("ABC.Score.",cell.type)] <- m$predictions[g$name,cell.type]
    }
    
    g$CredibleSets <- NA
    g$variants <- NA
    g$Diseases <- NA
    for (i in 1:nrow(g)) {
      g$CredibleSets[i] <- paste0(unique(subset(all.flat, TargetGene == g$TargetGene[i])$CredibleSet), collapse=',')
      g$variants[i] <- paste0(unique(subset(all.flat, TargetGene == g$TargetGene[i])$QueryRegionName), collapse=',')
      g$Diseases[i] <- paste0(unique(subset(all.flat, TargetGene == g$TargetGene[i])$Disease), collapse=',')
      g$CellTypes[i] <- paste0(unique(subset(all.flat, TargetGene == g$TargetGene[i])$CellType), collapse=',')
    }
    return(g)
  }
  
  
  getVariantAnnotationDf <- function(variant.list, all.flat) {
    all.flat$QueryRegionName <- factor(as.character(as.matrix(all.flat$QueryRegionName)), levels=levels(variant.list$variant))
    v <- variant.list
    v$variants <- NA
    v$Diseases <- NA
    for (i in 1:nrow(v)) {
      v$CredibleSets[i] <- paste0(unique(subset(all.flat, QueryRegionName == v$variant[i])$CredibleSet), collapse=',')
      v$variants[i] <- paste0(unique(subset(all.flat, QueryRegionName == v$variant[i])$QueryRegionName), collapse=',')
      v$Diseases[i] <- paste0(unique(subset(all.flat, QueryRegionName == v$variant[i])$Disease), collapse=',')
      v$CellTypes[i] <- paste0(unique(subset(all.flat, QueryRegionName == v$variant[i])$CellType), collapse=',')
    }
    return(v)
  }
  
  
  getCredibleSetDf <- function(all.cs, all.flat, control.types) {
    ## Columns:
    ## chr
    ## start
    ## end
    ## name
    ## Disease
    ## nSNP
    ## Coding (variant:gene pairs) .. for now, variant list
    ## Promoter (variant:gene pairs) .. for now, variant list
    ## Splice Site (variant:gene pairs) .. for now, variant list
    ## UTR (variant:gene pairs) .. skip for now
    ## CTCF Motif (variant list) .. skip for now
    ## BestSNPposition (SNP with highest p-value / Bayes posterior proability)
    ## BestSNP (name of SNP)
    ## BestSNPNearestGene (Gene closest to best SNP)
    ## LocusID (for linking credible sets in same locus)
    ## Num Strong Roadmap DHS overlaps .. skip for now
    ## Num any Roadmap DHS Overlaps .. skip for now
    ## Variant:CellType:DHS value for Roadmap DHS overlaps .. skip for now
    
    ## Variants with Immune-cell ABC predictions
    ## Genes with Immune-cell ABC predictions
    ## Cell types with Immune-cell ABC predictions
    ## Variants with Control-cell ABC predictions
    ## Genes with Control-cell ABC predictions
    ## Cell types with Control-cell ABC predictions
    
    coding <- unlist(by(variant.list, variant.list$CredibleSet, function(x) {
      if (any(x$Coding)) paste0(x$variant[x$Coding], collapse=',')
      else ""
    }, simplify=FALSE))
    all.cs$Variants.Coding <- coding[as.character(as.matrix(all.cs$CredibleSet))]
    
    ss <- unlist(by(variant.list, variant.list$CredibleSet, function(x) {
      if (any(x$SpliceSite)) paste0(x$variant[x$SpliceSite], collapse=',')
      else ""
    }, simplify=FALSE))
    all.cs$Variants.SpliceSite <- ss[as.character(as.matrix(all.cs$CredibleSet))]
    
    promoter <- unlist(by(variant.list, variant.list$CredibleSet, function(x) {
      if (any(x$Promoter)) paste0(x$variant[x$Promoter], collapse=',')
      else ""
    }, simplify=FALSE))
    all.cs$Variants.Promoter <- promoter[as.character(as.matrix(all.cs$CredibleSet))]
    
    immune.flat <- subset(all.flat, !(CellType %in% control.types))
    control.flat <- subset(all.flat, CellType %in% control.types)
    all.cs$Variants.ABC.Immune <- sapply(all.cs$CredibleSet, function(cs) {
      curr <- subset(immune.flat, CredibleSet == cs)
      if (nrow(curr) == 0) return("")
      return(paste0(unique(curr$QueryRegionName), collapse=','))
    })
    all.cs$Genes.ABC.Immune <- sapply(all.cs$CredibleSet, function(cs) {
      curr <- subset(immune.flat, CredibleSet == cs)
      if (nrow(curr) == 0) return("")
      return(paste0(unique(curr$TargetGene), collapse=','))
    })
    all.cs$CellTypes.ABC.Immune <- sapply(all.cs$CredibleSet, function(cs) {
      curr <- subset(immune.flat, CredibleSet == cs)
      if (nrow(curr) == 0) return("")
      return(paste0(unique(curr$CellType), collapse=','))
    })
    
    all.cs$Variants.ABC.Control <- sapply(all.cs$CredibleSet, function(cs) {
      curr <- subset(control.flat, CredibleSet == cs)
      if (nrow(curr) == 0) return("")
      return(paste0(unique(curr$QueryRegionName), collapse=','))
    })
    all.cs$Genes.ABC.Control <- sapply(all.cs$CredibleSet, function(cs) {
      curr <- subset(control.flat, CredibleSet == cs)
      if (nrow(curr) == 0) return("")
      return(paste0(unique(curr$TargetGene), collapse=','))
    })
    all.cs$CellTypes.ABC.Control <- sapply(all.cs$CredibleSet, function(cs) {
      curr <- subset(control.flat, CredibleSet == cs)
      if (nrow(curr) == 0) return("")
      return(paste0(unique(curr$CellType), collapse=','))
    })
    
    all.cs$nVariants.ABC.Immune <- sapply(all.cs$Variants.ABC.Immune, function(v) if (is.na(v)) 0 else length(strsplit(v,",")[[1]]))
    all.cs$nVariants.ABC.Control <- sapply(all.cs$Variants.ABC.Control, function(v) if (is.na(v)) 0 else length(strsplit(v,",")[[1]]))
    
    all.cs$nGenes.ABC.Immune <- sapply(all.cs$Genes.ABC.Immune, function(v) if (is.na(v)) 0 else length(strsplit(v,",")[[1]]))
    all.cs$nGenes.ABC.Control <- sapply(all.cs$Genes.ABC.Control, function(v) if (is.na(v)) 0 else length(strsplit(v,",")[[1]]))
    
    all.cs$nCellTypes.ABC.Immune <- sapply(all.cs$CellTypes.ABC.Immune, function(v) if (is.na(v)) 0 else length(strsplit(v,",")[[1]]))
    all.cs$nCellTypes.ABC.Control <- sapply(all.cs$CellTypes.ABC.Control, function(v) if (is.na(v)) 0 else length(strsplit(v,",")[[1]]))
    
    return(all.cs)
  }
  
  
  
  # Function to pull flanking sequence. Defaults to +/- 10 bp
  getSequenceWithSNP <- function(region.chr="chr1", snp.position=114433946, region.start=114433946-10, region.end=114433946+10, alleles=c("G","A")) {
    ## Temporary method -- doesn't deal with indels.  Use GATK method instead
    require('BSgenome.Hsapiens.UCSC.hg19')
    seq  <- as.character(getSeq(Hsapiens,region.chr,region.start,region.end))
    paste0(substr(seq,0,snp.position-region.start),alleles,substr(seq,snp.position-region.start+2,nchar(seq)))
  }
  
  
  
  
  loadVariantOverlapMatrix <- function(variant.list, dir, cell.types, suffix=".overlaps.txt", score.col=2) {
    ## Load in overlaps between variants and DHS sites, generated in log.sh
    ## JE 3/5/18
    
    variants <- levels(variant.list$variant)
    
    dhs.overlap <- matrix(0, nrow=length(variants), ncol=length(cell.types)) 
    rownames(dhs.overlap) <- variants
    colnames(dhs.overlap) <- cell.types
    
    ## Loop through cell types and annotate matrix with the designated score column
    for (cell in cell.types) {
      file <- paste0(dir, "/", cell, suffix)
      if (file.exists(file)) {
        curr <- read.delim(file, header=F)
        
        ## Collapse by variant to the highest score
        variant.scores <- tapply(curr[,score.col], curr[,1], max)
        variant.scores.df <- data.frame(variant=factor(names(variant.scores), levels=variants), score=variant.scores)
        
        ## Add to matrix - uses factor index (see variants and variant.scores.df above) to make it quick
        dhs.overlap[variant.scores.df$variant,cell] <- variant.scores.df$score
      } else {
        dhs.overlap <- dhs.overlap[,-which(colnames(dhs.overlap) == cell)]
      }
    }
    
    dhs.overlap <- data.frame(variant=levels(variant.list$variant), dhs.overlap)
    
    return(dhs.overlap)
  }
  
  plotCredibleSetCellTypeHeatmapFromVariantOverlap <- function(variant.list, overlap.matrix, n.breaks=100, binarize=FALSE, ...) {
    require(RColorBrewer)
    
    overlap.matrix <- merge(variant.list[,c("variant","CredibleSet")], overlap.matrix, by="variant")
    
    ## For each credible set, select strongest value for each cell type
    tab <- do.call(rbind, by(overlap.matrix[,-1:-2], overlap.matrix$CredibleSet, function(x) apply(x, 2, max), simplify=FALSE))
    #rownames(tab) <- unique(overlap.matrix$CredibleSet)
    #colnames(tab) <- colnames(overlap.matrix)[-1:-2]
    #browser()
    
    
    par(mar=c(12,0,4,8))  
    
    breaks <- quantile(as.numeric(tab)[as.numeric(tab) > 0], seq(0,1,1/n.breaks))
    palette <- colorRampPalette(c("white","red"))(n = length(breaks)-1)
    
    if (binarize) {
      tab[tab > 0] <- 1
      breaks <- c(0,0.5,1)
      palette <- c("white","red")
    }
    tab.sorted <- plotCredibleSetCellTypeHeatmapHelper(tab, histogram=TRUE, breaks=breaks, col=palette, margins=c(16,16), ...)
    return(tab.sorted)
  }
  
  
}
