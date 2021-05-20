import os,sys
import pandas as pd
import glob
import numpy as np
from subprocess import check_call, check_output, PIPE, Popen, getoutput, CalledProcessError

def grabStatistics(arr_values):
    """
    Takes in numpy array and generates mean, median, std values of that array
    """
    mean = arr_values.mean()
    median = arr_values.median()
    std = arr_values.std()
    return mean, median, std

def GetUniqueBpConnections(prediction):
    """
    Merge enhancer regions in prediction files 
    1. Calculate total number of unique bp EG connections 
    """
    os.system("cat {} | sed '1d' | bedtools merge -i stdin > {}_merged.tsv".format(prediction, prediction))

def checkFileExists(predFile):
    if os.path.isfile(predFile):
        print("File exists : {}".format(predFile))

def EnhancerGenePredictionProperties(predFile, expressedGenes, outdir):
    """
    Takes as input: 
    predFile : prediction file name 
    expressedGenes : list of expressed genes in that biosample 
    outdir : Directory to save output files

    This function answers questions regarding E-G link properties: 
    1. Number of enhancer in this biosample 
    2. Size of enhancers in this biosample 
    3. Distance of Enhancer-Gene Links in this biosample 
    4. Number of enhancer-gene connections per gene in this biosample 
    5. Number of genes per Enhancer in this biosample
    6. Get Unique Basepair Connections  
    7. Median values of 2-5
    """
    # Read in prediction file 
    prediction = pd.read_csv(predFile, sep="\t")

    EnhancerPerGene = prediction.groupby(['TargetGene']).size()
    prediction['SizeEnhancers'] = prediction['end'] - prediction['start'] 
    EnhancerRegions = EnhancerPerGene.drop_duplicates(['chr', 'start', 'end'])

    # Get Number and Size of Enhancers 
    NumEnhancers = len(EnhancerRegions)
    SizeEnhancers = list(EnhancerRegions['SizeEnhancers'])
    MedianSizeEnhancers = np.median(EnhancerRegions['SizeEnhancers'])

    # Grab Number of Enhancers Per Gene
    GeneMean, GeneMedian, GeneStdev = grabStatistics(EnhancerPerGene)
    # Grab Number of Enhancers Per Expressed Gene 
    prediction_expressed = prediction.loc[prediction['TargetGene'].isin(list(expressedGenes))]
    prediction_expressed.to_csv("{}_expressed.csv".format(predFile), sep="\t")
    ExpressedEnhancerPerGene = prediction_expressed.groupby(['TargetGene']).size()
    ExpressedGeneMean, ExpressedGeneMedian, ExpressedGeneStdev = grabStatistics(ExpressedEnhancerPerGene)
    # Get Merged Enhancer Gene Regions 
    GetUniqueBpConnections(predFile)
    GetUniqueBpConnections("{}_expressed.csv".format(predFile))

    # Check that files have been created 
    checkFileExists(predFile)
    checkFileExists("{}_expressed.csv".format(predFile))

    # Grab Number of genes per enhancers
    NumGenesPerEnhancer = prediction_df[['chr', 'start', 'end']].groupby(['chr', 'start', 'end']).size()
    mean_genes_per_enhancer, median_genes_per_enhancer,stdev_genes_per_enhancer = grabStatistics(NumGenesPerEnhancer) 
    
    # Grab Number of Enhancer-Gene Pairs Per Chromsome
    enhancergeneperchrom = prediction_df.groupby(['chr']).size()
    mean_enhancergeneperchrom, median_enhancergeneperchrom ,stdev_enhancergeneperchrom = grabStatistics(enhancergeneperchrom) 
    
    # Enhancer-Gene Distancee
    distance = np.array(prediction_df['distance'])
    thquantile = np.percentile(distance, 10)
    testthquantile = np.percentile(distance, 90)

    # Genes 
    genes = prediction['TargetGene'].drop_duplicates()

    # Get number of unique bp in EG connections and E-ExpressedGene connections 
    data = pd.read_csv("{}_merged.tsv".format(predFile), sep="\t", header=None)
    data['bp'] = data[2] - data[1]
    dataExpressed = pd.read_csv("{}_mergedExpressed.tsv".format(predFile), sep="\t", header=None)
    numUniqueBp = np.sum(data['bp'])
    numUniqueBpExpressed = np.sum(dataExpressed['bp'])

    # Making property dataframe 
    property_df = pd.DataFrame()
    property_df['NumEnhancers'] = NumEnhancers
    property_df['SizeEnhancers'] = [SizeEnhancers]
    property_df['MedianSizeEnhancers'] = MedianSizeEnhancers
    property_df['NumEnhancersPerGene'] = [EnhancerPerGene]
    property_df['MedianNumEnhancersPerGene'] = GeneMedian
    property_df['NumEnhancersPerExpressedGene'] = [ExpressedEnhancerPerGene]
    property_df['MedianNumEnhancersPerExpressedGene'] = ExpressedGeneMedian
    property_df['NumGenesPerEnhancer'] = [NumGenesPerEnhancer]
    property_df['MedianNumGenesPerEnhancer'] = median_genes_per_enhancer
    property_df['MedianEGPerChrom'] = median_enhancergeneperchrom
    property_df['NumUniqueBpInEGLinks'] = numUniqueBp
    property_df['NumUniqueBpInEExpressedGLinks'] = numUniqueBpExpressed
    property_df['Distance'] = [distance]
    property_df['MedianEGDist'] = np.median(distance)
    property_df.to_csv("{}_EnhancerGeneProperties.tsv", sep="\t", index=False)

    return genes
