library(ggplot2)
library(dplyr)
library(tidyr)
library(optparse)

######################################################

option_list <- list(
		    make_option(c("--names"), type="character", default=NA, help="list of prediction names"), 
		    make_option(c("--tables"), type="character", default=NA, help="list of enrichment tables"),
		    make_option(c("--outPdf"), type="character", default=NA, help="output file for CDF plot"),
		    make_option(c("--outDensity"), type="character", default=NA, help="output file for Aggregate Density Plot"),
		    make_option(c("--outDir"), type="character", default=NA, help="out directory"),
		    make_option(c("--filterCellTypes"), type="character", default=NA, help="Celltypes to use to plot enrichment")
		    )

opt <- parse_args(OptionParser(option_list=option_list))
names = strsplit(opt$names, " ") %>% unlist()
print(names)
enrichmentTables = strsplit(opt$tables, " ") %>% unlist()
outDir = opt$outdir

# aggregate data
enr.all = read.csv(file=enrichmentTables[1], sep='\t', header=TRUE, stringsAsFactors = FALSE)
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
		      
