# Plot precision-recall at identifying known IBD genes


suppressPackageStartupMessages(library(optparse))


option_list <- list(
        make_option("--genePredTable", help="Input gene predictions file (for one or more traits)"),
        make_option("--outPdf", type="character", help="Output PDF for precision-recall plot"),
        make_option("--knownGenes", type="character", help="File with columns corresponding to lists of genes with which to evaluate predictors. Columns to use must start with 'GeneList.'"),
        make_option("--codeDir", type="character", default="ABC-Max-pipeline/", help="code directory"),
        make_option("--knownGeneMaxDistance", type="numeric", default=1000000, help="Maximum distance from credible set to TSS to count a known gene")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))

#######################################################################
## Check inputs

checkInputs <- function(opt) {
  names = strsplit(opt$genePredTable, " ") %>% unlist()
  print(names)
  for (k in 1:length(names)){
      if (!file.exists(names[k])) stop(paste0("Gene prediction file does not exist: ", names[k]))
  }
  if (!file.exists(opt$knownGenes)) stop(paste0("Known genes file does not exist: ", opt$knownGenes))
  if (!file.exists(opt$codeDir)) stop(paste0("Pipeline code directory does not exist: ", opt$codeDir))
}
checkInputs(opt)

source(paste0(opt$codeDir, "/JuicerUtilities.R"))
source(paste0(opt$codeDir, "/CredibleSetTools.R"))

saveProgress <- function() save.image(file=paste0(opt$outPdf,".RData"))
saveProgress()

print("Open files")
#######################################################################
## Load gene prediction files 

names = strsplit(opt$genePredTable, " ") %>% unlist()
# aggregate data 
gp <- read.delim(names[1], check.names=F, stringsAsFactors=F, comment.char='#')
write.table(gp, file="gp.tsv")
# merge cols
mergecols <- colnames(gp[, 1:10])
print(mergecols)
print("reading...")
if (length(names) > 1) {
	for (i in 2:length(names)){
		temp = read.delim(file=names[i], check.names=F, stringsAsFactors=F, comment.char='#')
		gp = merge(gp, temp, by=mergecols)
	}
}
knownGenes <- read.delim(opt$knownGenes, header=T, check.names=F, stringsAsFactors=F, comment.char='#')
predictors <- colnames(gp)[greplany(c("GeneScore.","GenePrediction.","GenePredictionMax."), colnames(gp))]
print(knownGenes)

#######################################################################
## Functions for plotting PR curves

mytheme <- theme_classic() + theme(
  axis.text = element_text(size = 13), 
  axis.title = element_text(size = 15))

####################################################################
## PLOT:

pdf(file=opt$outPdf, width=6, height=4, onefile=T)
for (gl in colnames(knownGenes))
  gp.plot <- doOneKnownGeneList(gl, gp, predictors, maxKnownGenes=length(knownGenes))
  nRecall <- sum(gp.plot$knownGene)
  pr <- getPRTable(gp.plot, predictors)
  singlePR  <- getPRPlot(pr, baseline=getPrecisionBaseline(gp), xlab=paste0("Recall (n=",nRecall,")"))
dev.off()



