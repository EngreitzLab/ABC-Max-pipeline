# coding: utf-8

# This snakefile contains the rules for conducting the ABC-Max analysis as in
# [citation].

configfile: "ABC-Max.config.json"

from os.path import join

#def getCSnames(wildcards):
#	return ["a list of credible set names depending on given wildcards"]
#
#def getTraits(wildcards):
#	return ["a list of dz traits depending on given wildcards"]

# Gathering all the outputs for all sets of predictions and variants
outputSet = set()

# Outputs specific for predictions
for pred in config["predictions"]:
	OverlapAllSNPs = os.path.join([pred]["outDir"][0], pred, ".OverlapAllSNPs.tsv.gz")
	OverlapCounts = os.path.join([pred]["outDir"][0], pred, ".OverlapCounts.tsv")
	OverlapCountsAllNoncoding = os.path.join([pred]["outDir"][0], pred, ".OverlapCounts.AllNoncoding.tsv")
	outputSet.add(OverlapAllSNPs, OverlapCounts, OverlapCountsAllNoncoding) 

# Outputs specific for variant sets
for trait in config["traits"]:
	bedFile = os.path.join(config[trait]["outDir"][0], trait, ".bed")
	bedgraph = os.path.join(config[trait]["outDir"][0], trait, ".bedgraph")

	outputSet.add(bedFile, bedgraph)

# Outputs for predictions and variant sets
for pred in config["predictions"]:
	for trait in config["traits"]:
		overlap = os.path.join(config[trait]["outDir"][0], trait, ".", pred, ".tsv.gz")
		outputSet.add(overlap)

# Preprocessing ABC?
if "ABC" in config["predictions"]:
	outputSet.add(config["ABC"]["shrunkPredFile"][0])

rule all:
	input:
		outputSet

# Move this to another snakefile?
rule preprocessABC:
	input:
		rawPredFile = config["ABC"]["rawPredFile"][0]
	output:
		ABC015predFile = config["ABC"]["ABC015PredFile"][0],
		shrunkPredFile = config["ABC"]["shrunkPredFile"][0]
	params:
		chrSizes = config["chrSizes"]
	message: ""
	run:
		shell(
			"""
			# Removing promoter elements, retaining predictions with an ABC score >0.15
			zcat {input.rawPredFile} | sed 1d | awk '(($5 != "promoter" && $21>=0.015) || ($5 == "promoter" && $21>=0.1)) && $21 != "NaN" && int($12) <= 2000000'| bedtools sort -i stdin -faidx {params.chrSizes} | gzip > {input.ABC015predFile}

			# Shrinking elements from 500bp to 200bp
			zcat {input.ABC015predFile} | head -1 | gzip > {input.shrunkPredFile}
			zcat {input.ABC015predFile} | sed 1d | bedtools slop -b -150 -g {params.chrSizes} | gzip >> {input.shrunkPredFile}
			""")

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
			# Intersecting a background list of variants with predicted enhancers 
			# to compute background rate at which common variants overlap enhancers
			# overall
			
			# Compute fraction of variants overlapping predictions in each cell type
			# Finding the relevant columns
			chrCol=zcat {input.predFile} | awk -v RS='\t' '/chr/{print NR; exit}';
			startCol=zcat {input.predFile} | awk -v RS='\t' '/start/{print NR; exit}';
			endCol=zcat {input.predFile} | awk -v RS='\t' '/end/{print NR; exit}';
			CelltypeCol=zcat {input.predFile} | awk -v RS='\t' '/Celltype/{print NR; exit}';

			zcat {input.predFile} | cut -f $chrCol,$startCol,$endCol,$CelltypeCol | sed 1d | sort | uniq | bedtools sort -i stdin -faidx {input.chrSizes} | \
			bedtools intersect -sorted -g {input.chrSizes} -a {input.allVariants} -b stdin -wa -wb | gzip > {output.overallOverlap}
			
			# Getting the cell type column and counting
			zcat {output.overallOverlap} | cut -f 7 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > {output.overallOverlapCounts}

			# Compute fraction of noncoding variants overlapping predictions in any cell type
			zcat {output.overallOverlap} | bedtools intersect -v -a stdin -b {input.CDS} | cut -f 1-3,7 | sort | uniq | cut -f 4 | 	sort | uniq -c | sed 's/^ *//' | tr ' ' '\t' > {output.noncodingOverlap}
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
		varFilterCol = varFilterCol,
		varFilterThreshold = varFilterThreshold,
		sigvarList = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.sig.varList.tsv"),
		outDir = lambda wildcard: config[wildcard.trait]["outDir"][0],
		paramsFile = config["params"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers"
	run:
		if wildcard.varFilterCol is not None:
			shell(
				"""
				## Subsetting the variant list based on significance
				# Finding the score column
				scoreCol=$(awk -v RS='\t' '/{params.varFilterCol}/{print NR; exit}' {input.varList});

				# Filtering to retain variants exceeding the threshold
				awk '{{ if ($($scoreCol) >= {params.varFilterThreshold}) {{ print }} }}' {input.varList}  > {params.sigvarList};
	
				# Creating the bed file
				# Finding and cutting chr, pos, and var columns
				chrCol=zcat {params.sigvarList} | awk -v RS='\t' '/chr/{print NR; exit}';
				posCol=zcat {params.sigvarList} | awk -v RS='\t' '/position/{print NR; exit}';
				varCol=zcat {params.sigvarList} | awk -v RS='\t' '/variant/{print NR; exit}';
				cut -f $chrCol,$posCol,$varCol {params.sigvarList} | sed '1d' | awk -F'\t' '$1 = $1 FS $2-1' > {output.varBed};
	
				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cut -f 1-4 {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				zcat {input.predFile} | head -1 | awk '{ print $0 "\tvariant.chr\tvariant.start\tvariant.end\tQueryRegionName" }' | gzip > {output.overlap};
	
				# Intersecting variants with predictions
				zcat {input.predFile} | sed 1d | bedtools intersect -sorted -g {params.chrSizes} -b {output.varBedgraph} -a stdin -wb | gzip >> {output.overlap}
	
				""")
		else:
			shell(
				"""
				# Creating the bed file for all variants
				# Finding and cutting chr, pos, and var columns
				chrCol=zcat {params.sigvarList} | awk -v RS='\t' '/chr/{print NR; exit}';
				posCol=zcat {params.sigvarList} | awk -v RS='\t' '/position/{print NR; exit}';
				varCol=zcat {params.sigvarList} | awk -v RS='\t' '/variant/{print NR; exit}';
				cut -f $chrCol,$posCol,$varCol {input.varList} | sed '1d' | awk -F'\t' '$1 = $1 FS $2-1' > {output.varBed};
	
				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cut -f 1-4 {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				zcat {input.predFile} | head -1 | awk '{ print $0 "\tvariant.chr\tvariant.start\tvariant.end\tQueryRegionName" }' | gzip > {output.overlap};
	
				# Intersecting variants with predictions
				zcat {input.predFile} | sed 1d | bedtools intersect -sorted -g {params.chrSizes} -b {output.varBedgraph} -a stdin -wb | gzip >> {output.overlap}
	
				""")


rule overlapVariantsByCelltype:
	input:
		allPred = lambda wildcard: config[wildcard.pred]["predFile"][0],
		#celltypeAnnotation = lambda wildcard: config[wildcard.pred]["celltypeAnnotation"][0],
		varBedgraph = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.bedgraph"),
	output:
		overlap = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.{pred}.CellType.tsv.gz")
	log: os.path.join(config["logDir"][0], ".log")
	params:
		#outDir = lambda wildcard: config[wildcard.trait]["outDir"][0],
		# A file with paths to BED files containing predicted enhancers per cell type
		# Can be tested once all ABC outputs have been moved over
		enhancerListBeds = lambda wildcard: config[wildcard.pred]["enhancerListBeds"][0],
		chrSizes = config["chrSizes"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers by cell type"
	run:
		shell(
			"""
			# 


			# Creating an empty file with additional columns
			zcat {input.allPred} | head -1 | awk '{ print $0 "\tvariant.chr\tvariant.start\tvariant.end\tQueryRegionName" }' | gzip > {output.overlap} \
			
			# Overlapping elements (neighborhoods, all peaks, or all peaks without ABC)
			# Looping through cell types
			cat {params.paramsFile} | while read cellType bedFile; do
				cat $bedFile | cut -f 1-3 | uniq | awk -v ct=$cellType -F $'\t' '{{
				print $1 "\t" $2 "\t" $3 "\t" ct "-" NR "\tintergenic\t2\tGATA1\t48644981\tNaN\t0.9\tTrue\t100000\tFalse\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\t0.03\t0.5\t0.03\t" ct
			  }}' | bedtools intersect -sorted -g {input.chrSizes} -b {input.varBedgraph} -a stdin -wb
			done | bedtools sort -i stdin -faidx {input.chrSizes} | gzip >> {output.overlap}
			"""
			)

rule annotateVariants:
	input:
		varList = lambda wildcard: config[wildcard.trait]["varList"][0],
		predOverlapFile = os.path.join(lambda wildcard: config[wildcard.trait]["outDir"][0], "{trait}.{pred}.tsv.gz"),
		bgVars = config["bgVariants"]
	output:
		# ...
	log: os.path.join(config["logDir"][0], ".log")
	params:
		cellTypeTable = lambda wildcard: config[wildcard.pred]["celltypeAnnotation"][0],
		codeDir = config["codeDir"],
		outDir = lambda wildcard: config[wildcard.trait]["outDir"][0],
		scoreCol = lambda wildcard: config[wildcard.trait]["varFilterCol"][0],
		scoreType = lambda wildcard: config[wildcard.trait]["varScoreType"][0],
		scoreThreshold = lambda wildcard: config[wildcard.trait]["varFilterThreshold"][0],
		ctrlThreshold = lambda wildcard: config[wildcard.trait]["varCtrlThreshold"][0]
	message: "Annotating {wildcards.trait} variants with {wildcards.method} predictions"
	run:
		# If using ABC predictions, plotting some additional features
		if pred=="ABC":
			shell(
				"""
				R CMD BATCH AnnotateCredibleSets.R \
				--variants {input.varList} \
				--predictionFile {input.predOverlapFile} \
				--outbase {params.outDir} \
				--credibleSets {input.csList} \
				--codeDir {params.codeDir} \
				--cellTypeTable {params.cellTypeTable} \
				--geneLists {params.geneLists}"
				""")
		else:
			shell(
				"""
				R CMD BATCH AnnotateCredibleSets.R \
				--variants {input.varList} \
				--predictionFile {input.predOverlapFile} \
				--isABC FALSE \
				--outbase {params.outDir} \
				--credibleSets {input.csList} \
				--codeDir {params.codeDir} \
				--variantScoreCol {params.scoreCol} \
				--scoreType {params.scoreType} \
				--variantScoreThreshold {params.scoreThreshold} \
				--variantCtrlScoreThreshold {params.ctrlThrehold} \
				--cellTypeTable {params.cellTypeTable} \
				--geneLists {params.geneLists}"
				""")



