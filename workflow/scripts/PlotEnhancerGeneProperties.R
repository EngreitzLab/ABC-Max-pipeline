library(ggplot2)
library(dplyr)
library(tidyr)
library(optparse)

######################################################

option_list <- list(
		    make_option(c("--predFile"), type="character", default=NA, help="Enhancer Gene Links File"), 
		    make_option(c("--threshold"), type="character", default=NA, help="Value to threshold significant and non-significant links"),
		    make_option(c("--prefix_outfile"), type="character", default=NA, help="prefix for output files"),
		    make_option(c("--outDir"), type="character", default=NA, help="out directory")
		    )

opt <- parse_args(OptionParser(option_list=option_list))
predFiles = strsplit(opt$predFile, " ") %>% unlist()
outDir = opt$outdir

# aggregate data
enr.all = read.csv(file=predFiles[1], sep='\t', header=TRUE, stringsAsFactors = FALSE)
enr.all$predictionSet = names[1]

for (i in 2:length(names)){
	temp = read.csv(file=enrichmentTables[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
	temp$predictionSet = names[i]
	enr.all = rbind(enr.all, temp)
}
if (!is.na(opt$filterCellTypes)){
	cells = read.csv(file=opt$filterCellTypes, sep='\t', header=TRUE, stringsAsFactors = FALSE)
	enr.all = enr.all[enr.all$CellType %in% cells$cell, ]
}

cdf = ggplot(enr.all, aes(enrichment.NoPromoters, col=predictionSet)) + stat_ecdf(geom = "step") + ylab('Cumulative fraction') + xlab('Enrichment (GWAS variants/all common variants)') + theme_minimal() + scale_color_discrete(name='Prediction set') + theme(text = element_text(size = rel(4)), legend.text=element_text(size=rel(2.5)))
d = ggplot(enr.all,aes(x=enrichment.NoPromoters, col=predictionSet)) + geom_density() + theme_minimal() + ylab('Density') + xlab('Enrichment (GWAS variants/all common variants)') + scale_color_discrete(name='Prediction set') + theme(text = element_text(size = rel(4)), legend.text=element_text(size=rel(2.5)))

write.table(enr.all, file=paste0(opt$outDir, "EnrichmentDensityPlots.tsv"))
pdf(file=opt$outPdf, width=7, height=5); print(cdf); dev.off()
pdf(file=opt$outDensity, width=7, height=5); print(d); dev.off()
		      
