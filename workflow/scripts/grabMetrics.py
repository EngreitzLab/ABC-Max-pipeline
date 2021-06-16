import pandas as pd 
import numpy as np
import sys, os 
import argparse
import time
from collections import Counter
from multiprocessing import Pool

# TODO: condense the number of files that is output from here
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--predFile', required=True, help="merged Enhancer Prediction file"), 
    parser.add_argument('--totalUniqueBp', required=True, help="metrics file to store # of total unique bp"),
    parser.add_argument('--genes', required=True, help="File with gene bodies"),
    parser.add_argument('--numGenes', required=True, help="metrics file to store # genes with a predicted enhancer in another biosample"),
    parser.add_argument('--numEGCounts', required=True, help="num EnhancerGene counts txtfile"),
    parser.add_argument('--numBiosamplesCounts', required=True, help="numBiosamples count textfile"),
    parser.add_argument('--outdir', required=True, help="outDir")
    args = parser.parse_args()
    return args

def grepGenes(inputs):
    predFile = inputs[0]
    gene = inputs[1]
    outFile = inputs[2]
    outGene = inputs[3]
    outdir = inputs[4]
    os.system('zcat {} | grep "{}" > {}/temp{}'.format(args.predFile, gene, outdir, gene))
    if os.stat("{}/temp{}".format(outdir, gene)).st_size != 0:
        temp = pd.read_csv("{}/temp{}".format(outdir,gene), sep="\t", names=['chr', 'start', 'end', 'Genes', 'numGenes', 'Biosample', 'numBiosample'])
        entry = np.sum(temp['end'].astype('int')-temp['start'].astype('int'))
        with open(outFile, "a") as f:
            f.write(str(entry))
            f.write("\n")
        f.close()

        if len(temp['numBiosample'].loc[temp['numBiosample']>1]) > 0:
            with open(outGene, "a") as f1:
                f1.write(str(gene))
                f1.write("\n")
            f1.close()
        os.system("rm {}/temp{}".format(outdir, gene))

def grabMetrics(args):
    # read in prediction file 
    data = pd.read_csv(args.predFile, sep="\t", names=['chr', 'start', 'end', 'Genes', 'numGenes', 'Biosample', 'numBiosample'])
    genes = pd.read_csv(args.genes, sep="\t")
    
    numGenes = 0
    print("Grabbing total number of unique basepairs")
    # Number of total unique basepairs contained in enhancers for each gene 
    starttime = time.time()
    inputs = []
    for gene in genes['name'].drop_duplicates():
        inputs.append((args.predFile, gene, args.totalUniqueBp, args.numGenes, args.outdir))

    with Pool(30) as p:
        p.map(grepGenes, inputs)
    print("Done")
    endtime = time.time()
    print("Amount of time taken to loop through 24k genes: {}".format(str(endtime-starttime)))

    # Given an enhancer-gene connection in a particular biosample, in how many other biosamples is there an overlapping enhancer connected to the same gene?
    print("Grabbing number of biosamples with same EG Links")
    data['numBiosamplesConnectedToSameGene'] = 0
    starttime1 = time.time()
    concatNumBiosamples = []
    numEGenes = []
    numEGVals = []
    for index, row in data.iterrows():
        genes = str(row['Genes']).split(",")
        if len(genes) > 0:
            uniq_genes = list(set(genes))
            if len(uniq_genes) > 1:
                genes_to_use = [i for i in uniq_genes if str(i) not in numEGenes]
                if not genes_to_use:
                    numEGVals = [len(data[['chr', 'start', 'end']].loc[data['Genes']==gene].drop_duplicates()) for gene in genes_to_use] 
                    with open(args.numEGCounts, "a") as f:
                        for i,j in zip(genes_to_use, numEGVals):
                            f.write(str(i))
                            f.write("\t")
                            f.write(str(j))
                            f.write("\n")
                        f.close()
                for i in uniq_genes:
                    numEGenes.append(i)
            else:
                if uniq_genes[0] not in numEGenes:
                    numEGenes.append(uniq_genes[0])
                    numEGVals = len(data[['chr', 'start', 'end']].loc[data['Genes']==uniq_genes[0]])
                    with open(args.numEGCounts, "a") as f:
                        f.write(str(uniq_genes[0]))
                        f.write("\t")
                        f.write(str(numEGVals))
                        f.write("\n")
                        f.close()
            biosamples = row['Biosample']
            genes_collections = Counter(genes)
            val = [int(v) for k,v in sorted(genes_collections.items()) if v>1]
            concatNumBiosamples.append(val)
            if val:
               # getting the number of similar genes can be a proxy to the number of biosamples
                data.loc[index, 'numBiosamplesConnectedToSameGene'] = str(val)
            with open(args.numBiosamplesCounts, "a") as f:
                for i in val:
                    f.write(str(i))
                    f.write("\n")
                f.close()
                    
        else:
            print("Only one biosample detected: {}".format(index))
    print("Done!")
    endtime1 = time.time()
    print("Finished grabbing all metrics:{}".format(str(endtime1-starttime1)))
if __name__=="__main__":
    args = parse_args()
    grabMetrics(args)
