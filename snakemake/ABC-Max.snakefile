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

# TODO: automatically create new directories for each trait+prediction set combo
# TODO: overlap by celltype

# Preprocessing ABC?
# TODO: change shrunkPredFile to predFile
# TODO: move over the raw ABC output
if "ABC" in config["predictions"]:
	outputSet.add(config["ABC"]["shrunkPredFile"][0])

rule all:
	input:
		#outputSet,
		expand(os.path.join(config["outDir"], "{pred}.OverlapAllSNPs.tsv.gz"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}.OverlapCounts.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}.OverlapCounts.AllNoncoding.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{trait}.bed"), trait=config["traits"]),
		expand(os.path.join(config["outDir"], "{trait}.bedgraph"), trait=config["traits"]),
		expand(os.path.join(config["outDir"], "{trait}.{pred}.tsv.gz"), trait=config["traits"], pred=config["predictions"])
		#expand(os.path.join(config["outDir"], "{trait}.{pred}.CellType.tsv.gz"), trait=config["traits"], pred=config["predictions"])

# TODO: Move this to another snakefile?
rule preprocessABC:
	input:
		rawPredFile = config["ABC"]["rawPredFile"][0]
	output:
		ABC015predFile = config["ABC"]["ABC015PredFile"][0],
		shrunkPredFile = config["ABC"]["shrunkPredFile"][0]
	params:
		chrSizes = config["chrSizes"]
	message: "Preprocessing ABC predictions"
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
		allVariants = config["bgVariants"],
		chrSizes = config["chrSizes"],
		CDS = config["CDS"]
	output:
		overallOverlap = os.path.join(config["outDir"], "{pred}.OverlapAllSNPs.tsv.gz"),
		overallOverlapCounts = os.path.join(config["outDir"], "{pred}.OverlapCounts.tsv"),
		noncodingOverlap = os.path.join(config["outDir"], "{pred}.OverlapCounts.AllNoncoding.tsv")
	log: os.path.join(config["logDir"], "{pred}.bgoverlap.log")
	params:
	message: "Overlapping background variants with predictions: {wildcards.pred}"
	run:
		shell(
			"""
			#set +o pipefail;
			# Intersecting a background list of variants with predicted enhancers 
			# to compute background rate at which common variants overlap enhancers
			# overall

			# TODO: if cell type column is not provided, use all predictions, or require Celltype col?
			
			# Compute fraction of variants overlapping predictions in each cell type
			# Finding the relevant columns
			zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType | sed 1d | sort -k 1,1 -k 2,2n | uniq | bedtools sort -i stdin -faidx {input.chrSizes} | \
			bedtools intersect -sorted -g {input.chrSizes} -a {input.allVariants} -b stdin -wa -wb | gzip > {output.overallOverlap};
			
			# Getting the cell type column and counting
			zcat {output.overallOverlap} | cut -f 7 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.overallOverlapCounts};

			# Compute fraction of noncoding variants overlapping predictions in any cell type
			zcat {output.overallOverlap} | bedtools intersect -v -a stdin -b {input.CDS} | cut -f 1-3,7 | sort | uniq | cut -f 4 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.noncodingOverlap}
			""")

rule createVarFiles:
	input:
		varList = lambda wildcard: config[wildcard.trait]["varList"][0]
	output:
		varBed = os.path.join(config["outDir"], "{trait}.bed"),
		varBedgraph = os.path.join(config["outDir"], "{trait}.bedgraph"),
		sigvarList = os.path.join(config["outDir"], "{trait}.sig.varList.tsv")
	log: os.path.join(config["logDir"], "{trait}.createbed.log")
	params:
		varFilterCol = lambda wildcard: config[wildcard.trait]["varFilterCol"][0],
		varFilterThreshold = lambda wildcard: config[wildcard.trait]["varFilterThreshold"][0],
		outDir = lambda wildcard: config["outDir"],
		chrSizes = config["chrSizes"]
		#paramsFile = config["params"]
	message: "Creating variant BED files"
	run:
		if "{params.varFilterCol}" is not None:
			shell(
				"""
				## Subsetting the variant list based on significance
				# Finding the score column
				#scoreCol=$(awk -v RS='\\t' '/{params.varFilterCol}/{{print NR; exit}}' {input.varList});

				# Filtering to retain variants exceeding the threshold
				#awk '{{ if ($($scoreCol) >= {params.varFilterThreshold}) {{ print }} }}' {input.varList}  > {output.sigvarList};

				cat {input.varList} | csvtk -t filter -f "{params.varFilterCol}>={params.varFilterThreshold}" > {output.sigvarList};
	
				# Creating the bed file
				# Finding and cutting chr, position, and variant columns
				# TODO: do not require start and stop, only position?
				cat {output.sigvarList} | csvtk cut -t -f chr,position,variant | sed '1d' | awk -F "\\t" "\$1 = \$1 FS \$2-1 FS \$2 FS \$3 FS" | cut -f1-4 > {output.varBed};

				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cat {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				""")
		else:
			shell(
				"""
				# Creating the bed file for all variants
				# Finding and cutting chr, pos, and var columns
				cat {input.varList} | csvtk cut -t -f chr,position,variant | sed '1d' | awk -F "\\t" "\$1 = \$1 FS \$2-1 FS \$2 FS \$3 FS" | cut -f1-4 > {output.varBed};

				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cat {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				""")
	

rule overlapVariants:
	input:
		predFile = lambda wildcard: config[wildcard.pred]["predFile"][0],
		varList = lambda wildcard: config[wildcard.trait]["varList"][0],
		varBed = os.path.join(config["outDir"], "{trait}.bed"),
		varBedgraph = os.path.join(config["outDir"], "{trait}.bedgraph")
	output:
		overlap = os.path.join(config["outDir"], "{trait}.{pred}.tsv.gz")
	log: os.path.join(config["logDir"], "{trait}.{pred}.overlap.log")
	params:
		chrSizes = config["chrSizes"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers"
	run:
		shell(
			"""
			zcat {input.predFile} | head -1 | awk '{ print $0 "\\tvariant.chr\\tvariant.start\\tvariant.end\\tQueryRegionName" }' | gzip > {output.overlap};
	
			# Intersecting variants with predictions
			zcat {input.predFile} | sed 1d | bedtools intersect -sorted -g {params.chrSizes} -b {input.varBedgraph} -a stdin -wb | gzip >> {output.overlap}
	
			""")


rule overlapVariantsByCelltype:
	input:
		allPred = lambda wildcard: config[wildcard.pred]["predFile"][0],
		#celltypeAnnotation = lambda wildcard: config[wildcard.pred]["celltypeAnnotation"][0],
		varBedgraph = os.path.join(config["outDir"], "{trait}.bedgraph")
	output:
		overlap = os.path.join(config["outDir"], "{trait}.{pred}.CellType.tsv.gz")
	log: os.path.join(config["logDir"], "{trait}.{pred}.overlap.celltype.log")
	params:
		# A file with paths to BED files containing predicted enhancers per cell type
		enhancerListBeds = lambda wildcard: config[wildcard.pred]["enhancerListBeds"][0],
		chrSizes = config["chrSizes"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers by cell type"
	run:
		shell(
			"""
			# TODO: Automatically create the celltype bed files from the prediction file
			# TODO: generalize for any set of predictions

			# Creating an empty file with additional columns
			zcat {input.allPred} | head -1 | awk '{{ print $0 "\tvariant.chr\tvariant.start\tvariant.end\tQueryRegionName" }}' | gzip > {output.overlap} \
			
			# Looping through cell types
			cat {params.enhancerListBeds} | while read cellType bedFile; do
				cat $bedFile | cut -f 1-3 | uniq | awk -v ct=$cellType -F $'\\t' '{{
				print $1 "\t" $2 "\t" $3 "\t" ct "-" NR "\tintergenic\t2\tGATA1\t48644981\tNaN\t0.9\tTrue\t100000\tFalse\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\t0.03\t0.5\t0.03\t" ct
			  }}' | bedtools intersect -sorted -g {input.chrSizes} -b {input.varBedgraph} -a stdin -wb
			done | bedtools sort -i stdin -faidx {input.chrSizes} | gzip >> {output.overlap}
			"""
			)

rule annotateVariants:
	input:
		varList = lambda wildcard: config[wildcard.trait]["varList"][0],
		predOverlapFile = os.path.join(config["outDir"], "{trait}.{pred}.tsv.gz"),
		bgVars = config["bgVariants"]
	output:
		# ...
	log: os.path.join(config["logDir"], "{trait}.{pred}.annotate.log")
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



