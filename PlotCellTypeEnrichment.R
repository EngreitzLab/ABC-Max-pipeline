# Jesse Engreitz
# Plot IBD gene x cell type enrichment table
# use R-3.4


suppressPackageStartupMessages(library("optparse"))

## To do -- make this more customizable
option.list <- list(
  make_option("--params", type="character", help="Disease params file with cell type enrichment stats and others (output by CellTypeEnrichment.R"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option("--cellTypes", type="character", help="Cell type annotation file", default="ABC-GWAS/data/CellTypes.Annotated.ABCPaper.txt"),
  make_option("--cellTypeEnrichments", type="character", default="plots/CellTypeEnrichment.tsv", help="File containing merged cell type enrichments across traits")
  make_option("--codeDir", type="character", default="/seq/lincRNA/RAP/GWAS/200406_ABCPaper/LanderLab-EP-Prediction/"))
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(devtools))

source(paste0(opt$codeDir, "/src/libs/JuicerUtilities.R"))
source(paste0(opt$codeDir, "/src/libs/VariantPrediction/CredibleSetTools.R"))

setwd(opt$outdir)
save.image(file=paste0(opt$outdir, "/IBD.RData"))

###################################################
## Load in data relevant to IBD and format for plotting

params <- read.delim(opt$params, stringsAsFactors=F)
params <- subset(params, IncludeInPaper)

ibddir <- subset(params, Disease == "IBD")$ABCOverlap

catOrder <- c("myeloid","Bcell","Tcell","hematopoietic","fibroblast","epithelial","other")
#catOrder <- c("myeloid","Bcell","Tcell","hematopoietic","epithelial","other")
catColors <- c("green","orange","blue","purple","pink","brown","gray"); names(catColors) <- catOrder

cell.type.annot.all <- read.delim(opt$cellTypes, stringsAsFactors=F)
cellEnrichment <- read.delim("plots/CellTypeEnrichment.tsv", check.names=F, stringsAsFactors=F)

mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15))

########
## Plot cell type enrichments for IBD
{
  ibdEnhancerList <- cell.type.annot.all
  enrichPlot <- merge(cellEnrichment, cell.type.annot.all[,c("CellType","Categorical.IBDTissueAnnotations2")], by="CellType")
  ibdEnhancerList <- merge(ibdEnhancerList, subset(enrichPlot, Disease == "IBD")[,c("CellType","LDSC.Enrichment","FM.vsGenomeEnrichment","FM.vsGenomeSignificant")], by="CellType")
  ibdEnhancerList$CellCat <- ibdEnhancerList$Categorical.IBDTissueAnnotations2
  ibdEnhancerList$CellType <- ordered(ibdEnhancerList$CellType, levels=rev(cell.type.annot.all[order(cell.type.annot.all$Categorical.IBDTissueAnnotations2),"CellType"]))

  formatEnrichmentBarPlot <- function(p) p + geom_boxplot() + geom_jitter(position=position_jitter(0.2), size=3) + ylim(0, 22) + mytheme + geom_hline(yintercept=1, linetype="dashed", color="gray") + scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[-5]) + ylab("Enrichment\n(IBD variants / all variants)")
  pdf(file=paste0(opt$outdir, "/CellTypeEnrichment.IBD.pdf"), width=15, height=4)
  p1 <- ggplot(ibdEnhancerList, aes(x=CellCat, y=FM.vsGenomeEnrichment, color=CellCat)) + ggtitle("ABC Enhancers") 
  p1 <- p1 %>% formatEnrichmentBarPlot()
  p2 <- ggplot(ibdEnhancerList, aes(x=CellCat, y=EnhancerList.vsGenomeEnrichment, color=CellCat)) + ggtitle("Candidate Enhancers") 
  p2 <- p2 %>% formatEnrichmentBarPlot()
  p3 <- ggplot(ibdEnhancerList, aes(x=CellCat, y=AllPeaks.vsGenomeEnrichment, color=CellCat)) + ggtitle("All accessible regions") 
  p3 <- p3 %>% formatEnrichmentBarPlot()
  p4 <- ggplot(ibdEnhancerList, aes(x=CellCat, y=AllPeaksMinusABC.vsGenomeEnrichment, color=CellCat)) + ggtitle("Non-ABC accessible regions") 
  p4 <- p4 %>% formatEnrichmentBarPlot()
  print(plot_grid(p1,p4,p2, labels=c('A','B','C'), nrow=1))
  ggsave(paste0(opt$outdir, "/CellTypeEnrichment.IBD.eps"), width=15, height=4)
  dev.off()
}

save.image(file=paste0(opt$outdir, "/IBD.RData"))


