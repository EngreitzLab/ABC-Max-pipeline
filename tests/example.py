#! bin/python3

import pandas as pd
import os, sys
from multiprocessing import Pool

# data contains the information of the trait, the variant list, and the variant creidble set 
data = pd.read_csv(sys.argv[1], sep="\t", header=None)
predictions = ["Granja2019", "ABC", "BLUEPRINT_201217-ChromHMM"]

outdir="/oak/stanford/groups/akundaje/kmualim/test_code/"

inputs = []
for die, variant, sets in zip(data[0], data[1], data[2]):
    inputs.append((die, variant, sets))

def runAnnotate(input):
    # trait 
    # like IBD, etc
    die = input[0]
    # variant file 
    variant = input[1]
    # credible set file 
    sets = input[2]
    pred = "Granja2019"
    os.system("Rscript AnnotateCredibleSets.R --outbase tests/{}/{}/ --variants {} --credibleSets {} --codeDir /oak/stanford/groups/akundaje/kmualim/github/ABC-Max-pipeline/Utilities/ --predictionFile /oak/stanford/groups/akundaje/kmualim/tests/{}/{}/{}.{}.tsv.gz --bgOverlap /oak/stanford/groups/akundaje/kmualim/AltPreds/{}/{}.OverlapAllSNPs.tsv.gz --trait {} --cellType FALSE --TargetGene TRUE --TargetGeneTSS FALSE".format(pred, die, variant, sets, pred, die, die, pred, pred, pred, die))
    pred = "BLUEPRINT_201217-ChromHMM"
    os.system("Rscript AnnotateCredibleSets.R --outbase tests/{}/{}/ --variants {} --credibleSets {} --codeDir /oak/stanford/groups/akundaje/kmualim/github/ABC-Max-pipeline/Utilities/ --predictionFile /oak/stanford/groups/akundaje/kmualim/tests/{}/{}/{}.{}.tsv.gz --bgOverlap /oak/stanford/groups/akundaje/kmualim/AltPreds/{}/{}.OverlapAllSNPs.tsv.gz --trait {} --cellType TRUE --TargetGene FALSE --TargetGeneTSS FALSE".format(pred, die, variant, sets, pred, die, die, pred, pred, pred, die))
    pred = "ABC"
    os.system("Rscript AnnotateCredibleSets.R --outbase tests/{}/{}/ --variants {} --credibleSets {} --codeDir /oak/stanford/groups/akundaje/kmualim/github/ABC-Max-pipeline/Utilities/ --predictionFile /oak/stanford/groups/akundaje/kmualim/tests/{}/{}/{}.{}.tsv.gz --bgOverlap /oak/stanford/groups/akundaje/kmualim/AltPreds/{}/{}.OverlapAllSNPs.tsv.gz --trait {}".format(pred, die, variant, sets, pred, die, die, pred, pred, pred, die))


for input in inputs:
    runAnnotate(input)
