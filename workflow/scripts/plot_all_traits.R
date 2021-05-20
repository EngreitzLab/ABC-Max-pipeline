suppressPackageStartupMessages({
	library(plyr)
	library(dplyr)
	library(ggplot2)
	library(tidyr)
	library(plotly)
	library(stringr)
	library(gtools)
	library(viridis)
	library(DT)
	library(optparse)
})

main <- function() {
	option_list <- list(
			    make_option(c("--enrichmentTables"), type="character", default=NA, help="output filepath"),
			    make_option(c("--outDir"), type="character", default=NA, help="enrichment table"),
			    make_option(c("--pred"), type="character", default=NA, help= "number of base pairs per enhancer set"),
			    make_option(c("--outfile"), type="character", default=NA, help="")
			    )
	opt <- parse_args(OptionParser(option_list=option_list))
	enrichmentTables = (opt$enrichmentTables) %>% strsplit(" ") %>% unlist()
	pred = (opt$pred)
	outDir = (opt$outDir)
	outfile = (opt$outfile)
	cellEnrichment <- read.delim(enrichmentTables[1], sep = "\t", header=TRUE, check.names=F, stringsAsFactors=F, row.names=NULL)
	cellEnrichment <- transform(cellEnrichment, n.FractionOverlap = n / total)
	cellEnrichment <- transform(cellEnrichment, n.noPromoters.FractionOverlap = n.NoPromoters / total.NoPromoters)
	
	if (length(enrichmentTables) > 1){
		for (i in 2:length(enrichmentTables)){
			tryCatch(
			{
				temp = read.delim(file=enrichmentTables[i], sep='\t', header=TRUE, stringsAsFactors = FALSE)
				temp <- transform(temp, n.FractionOverlap = n / total)
				temp <- transform(temp, n.noPromoters.FractionOverlap = n.NoPromoters / total.NoPromoters)
				cellEnrichment = rbind(cellEnrichment, temp)
			},
			error=function(cond) {
				message(paste("File does not exist::", i))
			},
			finally={
				message(paste("Processed:",i))
			})
		}
	}

	write.table(cellEnrichment, file=paste0(outDir, pred, "_aggregateTraitEnrichment.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
	
	cellEnrichmentSignificant <- cellEnrichment %>% filter(p.NoPromoters < 0.05)
	numCellTypesVal <- cellEnrichmentSignificant %>% group_by(Disease) %>% summarize(numCellTypes = length(unique(CellType)))
	dz <- unique(cellEnrichmentSignificant$Disease)
	numCellTypesDF <- data.frame(n=numCellTypesVal$numCellTypes, Disease=dz)
	write.table(numCellTypesDF, file=outfile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

main()
