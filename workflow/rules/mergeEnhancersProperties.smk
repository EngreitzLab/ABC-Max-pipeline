rule mergeEnhancerRegions: 
	input:
		predFile = lambda wildcard: config["predDir"]+str(preds_config_file.loc[wildcard.pred, "predFile"])
	params:
		cellType = lambda wildcard: str(preds_config_file.loc[wildcard.pred, "hasCellType"]), 
		outDir = os.path.join(config["outDir"], "{pred}")
	output:
		mergedPredFile = os.path.join(config["outDir"], "{pred}/{pred}.mergedEnhancerRegions.tsv.gz"),
	message: "Merging enhancer regions in Prediction File"
	run:
		shell(
			"""
			# TODO: find an alternative to deal with pipefail
			set +o pipefail;
			
			# Merge Enhancer Regions across Biosamples 
			# Answers the question: 
			# For a given region of the genome, how many genes does it regulate across ALL Biosamples 
			# Given an enhancer (in a particular biosample), in how many other biosamples are there overlapping enhancers called?
			zcat {input.predFile} | csvtk cut -t -f chr,start,end,TargetGene,CellType | sed 1d | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,4,5,5 -o collapse,count_distinct,collapse,count_distinct | gzip > {output.mergedPredFile}
 			
			""")

rule generateCountMetrics:
	input:
		mergedPredFile = os.path.join(config["outDir"], "{pred}/{pred}.mergedEnhancerRegions.tsv.gz")
	params:
		metrics = os.path.join(config["codeDir"], "grabMetrics.py"),
		genes = config["promoterActivityRef"],
		outDir = os.path.join(config["outDir"], "{pred}")
	output:
		numGenes = os.path.join(config["outDir"], "{pred}/{pred}.metrics.numGenes.tsv"),
		numEGCounts = os.path.join(config["outDir"], "{pred}/{pred}.metrics.numEGCounts.tsv"),
		numBiosamplesCounts = os.path.join(config["outDir"], "{pred}/{pred}.metrics.numBiosamplesCounts.tsv"),
		totalUniquebp = os.path.join(config["outDir"], "{pred}/{pred}.metrics.uniquebp.tsv")
	run:
		shell(
			"""
			# Generates output files, each answering a global statistic about the set of enhancer-gene connections
			# numGenes : Number of Genes that regulate an enhancer
			# numEGCounts : Number of Enhancer-Gene Connections
			# numBiosamplesCounts : For a given enhancer, in how many other biosamples does this region also regulate a gene?
			# totalUniquebp: How many total unique basepairs constitute all enhancers that regulate a gene?
			python {params.metrics} \
			--predFile {input.mergedPredFile} \
			--totalUniqueBp {output.totalUniquebp} \
			--genes {params.genes} \
			--numGenes {output.numGenes} \
			--numBiosamplesCounts {output.numBiosamplesCounts} \
			--numEGCounts {output.numEGCounts} \
			--outdir {params.outDir}
			""")
