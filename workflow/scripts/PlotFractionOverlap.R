# Jesse Engreitz
# Plot IBD gene x cell type enrichment table
# use R-3.4


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))

## To do -- make this more customizable
option_list <- list(
		    make_option("--outdir", type="character", default="test", help="Output directory"),
	            make_option("--outPdf", type="character", help="Output PDF file for enrichment"),
        	    make_option("--outEps", type="character", help="Output EPS file for enrichment"),
		    make_option(c("--cellTypes"), type="character", default="Test_data/CellTypes.Annotated.ABCPaper.txt", help="Cell type annotation file"),
	    	    make_option(c("--cellTypeEnrichments"), type="character", default="plots/CellTypeEnrichment.tsv", help= "File containing merged cell type enrichments across traits"),
  		    make_option(c("--codeDir"), type="character", default="ABC-Max-pipeline/", help="code directory"),
		    make_option(c("--entry"), type="character", default="enrichment", help="feature to plot"),    
		    make_option(c("--trait"), type="character", default="IBD", help="trait name")
		    )

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(cowplot))

source(paste0(opt$codeDir, "/JuicerUtilities.R"))
source(paste0(opt$codeDir, "/CredibleSetTools.R"))

setwd(opt$outdir)
save.image(file=paste0(opt$outdir, "/PlotCellTypeEnrichment.RData"))

###################################################
## Load in data relevant to IBD and format for plotting

cellEnrichment <- read.delim(opt$cellTypeEnrichments, sep = "\t", header=TRUE, check.names=F, stringsAsFactors=F, row.names=NULL)
cellEnrichment <- transform(cellEnrichment, n.FractionOverlap = n / total)
cellEnrichment <- transform(cellEnrichment, n.noPromoters.FractionOverlap = n.NoPromoters / total.NoPromoters)

mytheme <- theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(size = 13), axis.title = element_text(size = 15))

cell.type.annot.all <- read.delim(opt$cellTypes, sep = "\t", header=TRUE, stringsAsFactors=F)
########
## Plot cell type enrichments for IBD
{
  ibdEnhancerList <- cell.type.annot.all
  enrichPlot <- merge(cellEnrichment, cell.type.annot.all[,c("CellType","Categorical.IBDTissueAnnotations2")], by="CellType")
  catOrder <- sort(unique(cell.type.annot.all$Categorical.IBDTissueAnnotations2))
  n <- length(unique(catOrder))
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  catColors <- c("green","orange","blue","purple","pink","brown","gray", sample(color, n-7))
  if (opt$entry == 'enrichment.NoPromoters'){
	ibdEnhancerList <- merge(ibdEnhancerList, subset(enrichPlot, Disease == opt$trait)[,c("CellType", "n.noPromoters.FractionOverlap")], by="CellType")
	
#	  names(enrichPlot)[names(enrichPlot) == 'enrichment.NoPromoters'] <- 'enrichment'
  } else {
        ibdEnhancerList <- merge(ibdEnhancerList, subset(enrichPlot, Disease == opt$trait)[,c("CellType", "n.FractionOverlap")], by="CellType")
  }
  ibdEnhancerList$CellCat <- ordered(ibdEnhancerList$Categorical.IBDTissueAnnotations2, levels=catOrder)
  ibdEnhancerList$CellType <- ordered(ibdEnhancerList$CellType, levels=rev(cell.type.annot.all[order(cell.type.annot.all$Categorical.IBDTissueAnnotations2),"CellType"]))
  formatEnrichmentBarPlot <- function(p) p + geom_boxplot() + geom_jitter(position=position_jitter(0.2), size=3) + ylim(0, 0.20) + mytheme + geom_hline(yintercept=1, linetype="dashed", color="gray") + scale_color_manual(values=catColors) + scale_color_hue(l=50, c=90) + ylab(paste0("Fraction of \n", opt$trait, " variants overlap Enhancers"))
  
  pdf(file=opt$outPdf, width=5, height=4)
  if (opt$entry == 'enrichment.NoPromoters'){
	  p1 <- ggplot(ibdEnhancerList, aes(x=CellCat, y=n.noPromoters.FractionOverlap, color=CellCat)) + ggtitle("Enhancers")
} else {
  	  p1 <- ggplot(ibdEnhancerList, aes(x=CellCat, y=n.FractionOverlap, color=CellCat)) + ggtitle("Enhancers") 
  }
  p1 <- p1 %>% formatEnrichmentBarPlot()
  print(p1)
  
  ggsave(opt$outEps, width=5, height=4)

  dev.off()
}

# For debugging: save.image(file=paste0(opt$outdir, "/PlotCellTypeEnrichment.RData"))


