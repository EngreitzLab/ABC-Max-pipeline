# Jesse Engreitz
# Plot IBD gene x cell type enrichment table
# use R-3.4


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))

## To do -- make this more customizable
option.list <- list(
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option("--cellTypes", type="character", help="Cell type annotation file", default="Test_data/CellTypes.Annotated.ABCPaper.txt"),
  make_option("--cellTypeEnrichments", type="character", default="plots/CellTypeEnrichment.tsv", help="File containing merged cell type enrichments across traits"),
  make_option("--codeDir", type="character", default="ABC-Max-pipeline/"),
  make_option("--pred", type="character", default="ABC"))
opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(cowplot))

source(paste0(opt$codeDir, "/Utilities/JuicerUtilities.R"))
source(paste0(opt$codeDir, "/Utilities/CredibleSetTools.R"))

setwd(opt$outdir)
save.image(file=paste0(opt$outdir, "/IBD.RData"))

###################################################
## Load in data relevant to IBD and format for plotting

catOrder <- c("myeloid","Bcell","Tcell","hematopoietic","fibroblast","epithelial","other")
#catOrder <- c("myeloid","Bcell","Tcell","hematopoietic","epithelial","other")
catColors <- c("green","orange","blue","purple","pink","brown","gray"); names(catColors) <- catOrder

cell.type.annot.all <- read.delim(opt$cellTypes, sep = "\t", header=TRUE, stringsAsFactors=F)
cellEnrichment <- read.delim(opt$cellTypeEnrichments, sep = "\t", header=TRUE, check.names=F, stringsAsFactors=F, row.names=NULL)

mytheme <- theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(size = 13), axis.title = element_text(size = 15))

########
## Plot cell type enrichments for IBD
{
  ibdEnhancerList <- cell.type.annot.all
  enrichPlot <- merge(cellEnrichment, cell.type.annot.all[,c("CellType","Categorical.IBDTissueAnnotations2")], by="CellType")
  ibdEnhancerList <- merge(ibdEnhancerList, subset(enrichPlot, Disease == "IBD")[,c("CellType","enrichment")], by="CellType")
  ibdEnhancerList$CellCat <- ibdEnhancerList$Categorical.IBDTissueAnnotations2
  ibdEnhancerList$CellType <- ordered(ibdEnhancerList$CellType, levels=rev(cell.type.annot.all[order(cell.type.annot.all$Categorical.IBDTissueAnnotations2),"CellType"]))

  formatEnrichmentBarPlot <- function(p) p + geom_boxplot() + geom_jitter(position=position_jitter(0.2), size=3) + ylim(0, 22) + mytheme + geom_hline(yintercept=1, linetype="dashed", color="gray") + scale_color_manual(values=brewer.pal(n = 8, name = "Dark2")[-5]) + ylab("Enrichment\n(IBD variants / all variants)")
  pdf(file=paste0(opt$outdir, "/CellTypeEnrichment.", opt$pred, ".pdf"), width=5, height=4)
  p1 <- ggplot(ibdEnhancerList, aes(x=CellCat, y=enrichment, color=CellCat)) + ggtitle("ABC Enhancers") 
  p1 <- p1 %>% formatEnrichmentBarPlot()
  print(p1)
  ggsave(paste0(opt$outdir, "/CellTypeEnrichment.", opt$pred, ".eps"), width=5, height=4)
  dev.off()
}

save.image(file=paste0(opt$outdir, "/IBD.RData"))


