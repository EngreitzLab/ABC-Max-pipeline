suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("GenomicRanges"))

option.list <- list(
	make_option("--allFlat", type="character", help="get all.flat.tsv files generated from AnnotateCredibleSets.R"),
	make_option("--codeDir", type="character", default="/oak/stanford/groups/akundaje/kmualim/ABC-Max-pipeline/workflow/scripts/"),
	make_option("--datadir", type="character", help="Data directory"),
	make_option("--traitTable", type="character", help="name of trait parameter file"), 
	make_option("--trait", type="character", help="trait", default="IBD"),
	make_option(c("-o", "--outdir"), type="character", help="name for pdf output file"))

# example of these scripts 
# geneFeatures -- generated from log.SummaryProperties.R to get genePredictions == for all diseases, get directory of ABCOverlap
# QueryRegionName, Disease, PosteriorProb, CredibleSet, Coding, CellType
# From GenePredictions.all.tsv
# "Disease","CredibleSet","TargetGene","ConnectionStrengthRank.Binary.AnyDisease_FMOverlap_Enriched","DistanceRank","POPS.Rank"
# params == params for traits 
# necessary columns : Disease, ABCOverlapColumn
# cellTypeEnrich == invidivual celltype enrichment files 
# necessary columns : Disease, CellType, LDSC.Significant, FM.vsGenomeSignificant 


opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)

source(paste0(opt$codeDir, "/JuicerUtilities.R"))
source(paste0(opt$codeDir, "/CredibleSetTools.R"))

names = strsplit(opt$allFlat, " ") %>% unlist()
datadir = opt$datadir
trait = opt$trait

trait_params = read.delim(file=opt$traitTable, sep="\t", header=TRUE, stringsAsFactors=FALSE)

all.flat.full = read.csv(file=names[1], sep='\t', header=TRUE, stringsAsFactors = FALSE)
scorecol <- all.flat.full %>% select(contains("Score"))
all.flat.full.subset <- all.flat.full %>% select(QueryRegionName,chr,start,end,TargetGene,CellType,variant.chr,variant.start,variant.end,TargetGeneIsCoding,CredibleSet,Disease,PosteriorProb,Coding,SpliceSite,Promoter,LocusID)
all.flat.full.subset$Score <- scorecol 

if(length(names)>1){
	for (i in 2:length(names)){
		tryCatch(
		{
			temp = read.delim(file=names[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
			scorecol <- temp %>% select(contains("Score"))
			temp.subset <- temp %>% select(QueryRegionName,chr,start,end,TargetGene,CellType,variant.chr,variant.start,variant.end,TargetGeneIsCoding,CredibleSet,Disease,PosteriorProb,Coding,SpliceSite,Promoter,LocusID)
			temp.subset$Score <- scorecol
			all.flat.full.subset = rbind(all.flat.full.subset, temp.subset)
		}, 
		error=function(cond) {
			message(paste("File does not exist::", i))
		}, 
		finally={
			message(paste("Processed:",i))
		})
	}
}

all.flat.full <- all.flat.full.subset
all.flat.full <- all.flat.full %>% select(QueryRegionName,Disease,PosteriorProb,CredibleSet,Coding,CellType) %>% unique()
all.flat.full$CellType <- factor(all.flat.full$CellType)
all.flat.v <- all.flat.full %>% select(-CellType) %>% unique()
all.v <- readAllVariants(trait_params, datadir)

# seems specific to IBD 
#final.dzs <- setdiff(subset(params, IncludeInPaper)$Disease,c("CD","UC"))
#final.cs <- lapply(final.dzs, function(dz) subset(all.cs[[dz]], !AnyCoding & !AnySpliceSite & AnyPp10)); names(final.cs) <- final.dzs
final.flat.full <- NULL #merge(all.flat.full, do.call(rbind, lapply(final.cs, function(cs) cs %>% select(Disease, CredibleSet))))
final.flat.v <- NULL #final.flat.full %>% select(-CellType) %>% unique()

ppBreaks <- c(0,0.0001,0.001,seq(0.01,1,0.01))
dzs <- unique(all.flat.v$Disease); dzs <- dzs[!(dzs %in% c("IBD+UC","IBD+CD"))]
#enrichVsPp <- getAllVariantEnrichment(ppBreaks, all.flat.full=all.flat.full, all.flat.v=all.flat.v, final.flat.v=final.flat.v, vlist=all.v)
#enrichVsPp.byDisease <- do.call(rbind, lapply(dzs, function(dz) getAllVariantEnrichment(ppBreaks, all.flat.full=all.flat.full, all.flat.v=all.flat.v, final.flat.v=final.flat.v, vlist=all.v, dz=dz)))
#enrichVsPp <- rbind(enrichVsPp, enrichVsPp.byDisease)
#write.table(enrichVsPp, file=paste0(opt$outdir, "/EnrichmentVsPosteriorProb.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
enrichVsPp.noncoding <- getAllVariantEnrichment(ppBreaks, all.flat.full=all.flat.full, all.flat.v=all.flat.v, final.flat.v=final.flat.v, vlist=all.v, noncodingOnly=TRUE)
enrichVsPp.noncoding <- rbind(enrichVsPp.noncoding, do.call(rbind, lapply(dzs, function(dz) getAllVariantEnrichment(ppBreaks, all.flat.full=all.flat.full, all.flat.v=all.flat.v, final.flat.v=final.flat.v, vlist=all.v, dz=dz, noncodingOnly=TRUE))))
write.table(enrichVsPp.noncoding, file=paste0(opt$outdir, "/EnrichmentVsPosteriorProb.nonCoding.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# 04122021 KSM: set this as an input instead of hardcoded 
bloodRelatedTraits <- c("Baso","Eosino","Hb","LOY","Lym","MCH","MCHC","MCV","Mono","Neutro","RBC","WBC")
#enrichVsPp.noncoding.bloodTraits <- getAllVariantEnrichment(c(0.1,0.95), all.flat.full=all.flat.full, all.flat.v=all.flat.v, final.flat.v=final.flat.v, vlist=all.v, noncodingOnly=TRUE, dz=bloodRelatedTraits)

# set theme 
mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15))

# PLOT 
{
  pdf(file=paste0(opt$outdir, "/EnrichmentVsPosteriorProb.pdf"), width=5, height=5)
#  p <- ggplot(enrichVsPp, aes(x=log10(PosteriorProb), y=enrichment, group=Disease)) + geom_line(color='gray') + xlab("PIP (log10)") + ylab("Enrichment") + ylim(1,NA) + ggtitle("ABC Enhancers, all traits") + mytheme
#  p <- p + geom_line(data=subset(enrichVsPp,Disease=="All"), aes(x=log10(PosteriorProb), y=enrichment), color='black')
#  p <- p + geom_line(data=subset(enrichVsPp,Disease==trait), aes(x=log10(PosteriorProb), y=enrichment), color='red')
#  p <- p + geom_hline(yintercept=1, linetype='dashed', color='black')
#  p <- p + theme(legend.position="none")
#  print(p)
#
#  p <- ggplot(enrichVsPp, aes(x=log10(PosteriorProb), y=fraction, group=Disease)) + geom_line(color='gray') + xlab("PIP (log10)") + ylab("Fraction variants") + ylim(0,1) + ggtitle("ABC Enhancers, all traits") + mytheme
#  p <- p + geom_line(data=subset(enrichVsPp,Disease=="All"), aes(x=log10(PosteriorProb), y=fraction), color='black')
#  p <- p + geom_line(data=subset(enrichVsPp,Disease==trait), aes(x=log10(PosteriorProb), y=fraction), color='red')
#  p <- p + geom_hline(yintercept=0.07568, linetype='dashed', color='black')
#  p <- p + theme(legend.position="none")
#  print(p)

  p <- ggplot(subset(enrichVsPp.noncoding, nTotal >= 5 & PosteriorProb <= 0.95), aes(x=log10(PosteriorProb), y=fraction*100, group=Disease)) + geom_line(color='gray') + xlab("PIP (log10)") + ylab("% noncoding variants") + ylim(0,100) + ggtitle("ABC Enhancers, all traits") + mytheme
  p <- p + geom_line(data=subset(enrichVsPp.noncoding,Disease%in%bloodRelatedTraits & nTotal >= 5 & PosteriorProb <= 0.95), aes(x=log10(PosteriorProb), y=fraction*100), color='red')
  p <- p + geom_line(data=subset(enrichVsPp.noncoding,Disease=="All" & PosteriorProb <= 0.95), aes(x=log10(PosteriorProb), y=fraction*100), color='black')
  p <- p + geom_line(data=subset(enrichVsPp.noncoding,Disease==trait & nTotal >= 5 & PosteriorProb <= 0.95), aes(x=log10(PosteriorProb), y=fraction*100), color='blue')
  p <- p + geom_hline(yintercept=7.5, linetype='dashed', color='black')
  p <- p + theme(legend.position="none")
  print(p)
  dev.off()
}
