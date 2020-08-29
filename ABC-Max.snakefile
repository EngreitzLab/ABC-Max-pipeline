# coding: utf-8

# This snakefile contains the rules for conducting the ABC-Max analysis as in
# [citation].

configfile: "ABC-Max.config.json"

from os.path import join

## For any set of variants and predictions, overlap

# Compute fraction of variants overlapping predictions in each cell type?
# Compute fraction of non-coding variants overlapping predictions in each cell type?
# Compute fraction of non-coding variants overlapping predictions in any cell type?

## Credible set analysis: AnnotateCredibleSets.R

# For any set of variants and predictions, and for every trait and cell type category, annotate
# and run the enrichment analysis.

# For every trait and credible set, create heatmaps.
# For every trait and credible set, create gene prioritization tables (ABC-Max)

# Same analysis by cell type?

#def getCSnames(wildcards):
#	return ["a list of credible set names depending on given wildcards"]
#
#def getTraits(wildcards):
#	return ["a list of dz traits depending on given wildcards"]

# Gathering all the outputs for all sets of predictions and variants
outputSet = set()
predictions = config["predictions"]
traits = config["traits"]

# Outputs specific for predictions
for pred in predictions:
	# Do stuff
	OverlapAllSNPs = os.path.join([pred]["outDir"][0], pred, ".OverlapAllSNPs.tsv.gz")
	OverlapCounts = os.path.join([pred]["outDir"][0], pred, ".OverlapCounts.tsv")
	OverlapCountsAllNoncoding = os.path.join([pred]["outDir"][0], pred, ".OverlapCounts.AllNoncoding.tsv")
	outputSet.add(OverlapAllSNPs, OverlapCounts, OverlapCountsAllNoncoding) 

# Outputs specific for variant sets
for trait in traits:
	bedFile = os.path.join(config[trait]["outDir"][0], trait, ".bed")
	bedgraph = os.path.join(config[trait]["outDir"][0], trait, ".bedgraph")
	outputSet.add(bedFile, bedgraph)

# Outputs for predictions and variant sets
for pred in predictions:
	for trait in traits:
		overlap = os.path.join(config[trait]["outDir"][0], trait, ".", pred, ".tsv.gz")
		outputSet.add(overlap)

rule all:
	input:
		outputSet

# ToDo: add ABC preprocessing steps

rule computeBackgroundOverlap:
	input:
		predFile = lambda wildcard: config[wildcard.pred]["predFile"][0],
		allVariants = config["allVariants"],
		chrSizes = config["chrSizes"],
		CDS = config["CDS"]
	output:
		overallOverlap = os.path.join(lambda wildcard: config[wildcard.pred]["outDir"][0], "{pred}.OverlapAllSNPs.tsv.gz"),
		overallOverlapCounts = os.path.join(lambda wildcard: config[wildcard.pred]["outDir"][0], "{pred}.OverlapCounts.tsv"),
		noncodingOverlap = os.path.join(lambda wildcard: config[wildcard.pred]["outDir"][0], "{pred}.OverlapCounts.AllNoncoding.tsv")
	log: os.path.join(config["logDir"][0], ".log")
	params:
	message: "Overlapping background variants with predictions: {wildcards.pred}"
	run:
		shell(
			"""
			# Intersect background list of variants (e.g. all 1000 genomes variants) with predicted enhancers 
			# to compute background rate at which common variants overlap enhancers overall
			
			# Compute fraction of variants overlapping predictions in each cell type
			zcat {input.predFile} | cut -f 1-3,8,9 | sed 1d | sort | uniq | bedtools sort -i stdin -faidx {input.chrSizes} | \
			bedtools intersect -sorted -g {input.chrSizes} -a {input.allVariants} -b stdin -wa -wb | gzip > {output.overallOverlap}
			zcat {output.overallOverlap} | cut -f 7-8 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > {output.overallOverlapCounts}

			# Compute fraction of noncoding variants overlapping predictions in any cell type
			zcat {output.overallOverlap} | bedtools intersect -v -a stdin -b {input.CDS} | cut -f 1-3,8 | sort | uniq | cut -f 4 | 	sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > {output.noncodingOverlap}
			""")

rule overlapVariants:
	input:
		predFile = lambda wildcard: config[wildcard.pred]["predFile"][0],
		varList = lambda wildcard: config[wildcard.trait]["varList"][0]
	output:
		varBed = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.bed"),
		varBedgraph = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.bedgraph"),
		overlap = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.{pred}.tsv.gz")
	log: os.path.join(config["logDir"][0], ".log")
	params:
		outDir = lambda wildcard: config[wildcard.trait]["outDir"][0],
		paramsFile = config["params"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers"
	run:
		shell(
			"""
			# Creating the bed file
			cut -f1,2,3 {input.varList} | sed '1d' | awk -F'\t' '$1 = $1 FS $2-1' > {output.varBed}

			# Ensure that variants are sorted for bedtools -sorted overlap algorithm
			cut -f 1-4 {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
			zcat {input.predFile} | head -1 | awk '{ print $0 "\tvariant.chr\tvariant.start\tvariant.end\tQueryRegionName" }' | gzip > {output.overlap};

			# Intersecting variants with predictions
			zcat {input.predFile} | sed 1d | bedtools intersect -sorted -g {params.chrSizes} -b {output.varBedgraph} -a stdin -wb | gzip >> {output.overlap}

			""")

rule overlapVariantsByCelltype:
	input:
		allPred = lambda wildcard: config[wildcard.pred]["predFile"][0],
		celltypeAnnotation = lambda wildcard: config[wildcard.pred]["celltypeAnnotation"][0],
		#predOverlapFile = lambda wildcard: config[wildcard.pred]["predFile"][0],
		#varBed = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.bed"),
		varBedgraph = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.bedgraph"),
	output:
		overlap = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.{pred}.CellType.tsv.gz")
	log: os.path.join(config["logDir"][0], ".log")
	params:
		outDir = lambda wildcard: config[wildcard.trait]["outDir"][0],
		chrSizes = config["chrSizes"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers by cell type"
	run:
		shell(
			"""
			#  Looping through cell types
			zcat {input.allPred} | head -1 | awk '{ print $0 "\tvariant.chr\tvariant.start\tvariant.end\tQueryRegionName" }' | gzip > {output.overlap} \
			cat {params.paramsFile} | while read cellType bedFile; do
				cat $bedFile | cut -f 1-3 | uniq | awk -v ct=$cellType -F $'\t' '{
				print $1 "\t" $2 "\t" $3 "\t" ct "-" NR "\tintergenic\t2\tGATA1\t48644981\tNaN\t0.9\tTrue\t100000\tFalse\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\t0.03\t0.5\t0.03\t" ct
			  }' | bedtools intersect -sorted -g {input.chrSizes} -b {input.varBedgraph} -a stdin -wb
			done | bedtools sort -i stdin -faidx {input.chrSizes} | gzip >> {output.overlap}
			""")

rule annotateVariants:
	output: 
		#join(outDir, "/data/all.flat.tsv"),
		#join(outDir, "/data/filter.flat.tsv"),
		#join(outDir, "/data/pp10.flat.tsv"),
		#join(outDir, "/data/all.flat.tsv"),
		## For every trait in the variant list
		#join(outDir, "Enrichment.CellType.vsPP.", "{trait}", ".tsv")
		## For every cell type category
		#join(outDir, "Enrichment.CellType.vsPP.", "{trait}", ".", "{celltype}", ".tsv"),
		#join(outDir, "CellTypes.Annotated.txt"),
		## For every credible set
		#join(outDir, "CredibleSetHeatmaps/", "{cs.name}" ,".tsv"),
		#join(outDir, "CredibleSetHeatmaps.pdf")
		# ...
	log: os.path.join(config["logDir"][0], ".log")
	params:
		cellTypeTable = "",
		geneLists = "",
		codeDir = "",
		outDir = "",
		PP = 0.1
		#CTRLPP = 0.01
	message: "{wildcards.method}, {wildcards.trait}"
	run:
		shell(
			"""
			R CMD BATCH AnnotateCredibleSets.R \
			--variants {input.varList} \
			--predictionFile {input.predictionFile} \
			--outbase {params.outDir} \
			--credibleSets {input.csList} \
			--codeDir {params.codeDir} \
			--posteriorProb {params.PP} \
			--cellTypeTable {params.cellTypeTable} \
			--geneLists {params.geneLists}"
			""")