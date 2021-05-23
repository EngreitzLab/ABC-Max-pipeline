## Jesse Engreitz
## February 12, 2020
## Quick script to output summary statistics about E-G connections

library(dplyr)
library(tidyr)
source("/seq/lincRNA/RAP/Promoters/LanderLab-EP-Prediction/src/libs/JuicerUtilities.R")
library(GenomicRanges)

ABC="AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"

x <- read.delim(gzfile(ABC))


distal <- subset(x, class != "promoter")

totalEperG <- distal %>% group_by(TargetGene) %>% summarise(
	TotalEnhancers=n(), 
	AverageEnhancerDistance=mean(distance),
	TotalBases=sum(width(reduce(GRangesFromBed(data.frame(chr=chr, start=start, end=end))))))
cellTypesPerG <- x %>% group_by(TargetGene) %>% summarise(CellTypesWithPrediction=length(unique(CellType)))

totalDPperG <- subset(x, class == "promoter" & isSelfPromoter == "False") %>% group_by(TargetGene) %>% summarise(TotalDistalPromoters=n())

all(totalEperG$TargetGene == cellTypesPerG$TargetGene)
result <- merge(totalEperG, cellTypesPerG, all.x=TRUE, all.y=TRUE)
result <- merge(result, totalDPperG, all.x=TRUE)
result[is.na(result)] <- 0
result$EnhancersPerCellType <- with(result, TotalEnhancers / CellTypesWithPrediction)
result$EnhancerBasesPerCellType <- with(result, TotalBases / CellTypesWithPrediction)
result$DistalPromotersPerCellType <- with(result, TotalDistalPromoters / CellTypesWithPrediction)

write.tab(result, file="GeneEnhancerStats.ForABCPaperV3.tsv")


## Add in Roadmap metrics
result <- read.delim("GeneEnhancerStats.ForABCPaperV3.tsv", stringsAsFactors=F, check.names=F)
roadmap <- read.delim("EnhancersPerGenePerCelltypeRoadmapABC.tsv", stringsAsFactors=F, check.names=F)

m <- merge(result, roadmap %>% select(-which(grepl("ABC",colnames(roadmap)))), all.x=TRUE)

write.tab(m, file="GeneEnhancerStats.ForABCPaperV3.WithRoadmap.tsv")
