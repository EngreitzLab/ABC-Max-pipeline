import pandas as pd
import numpy as np
import argparse 
import os
import seaborn as sns
import matplotlib.pyplot as plt

def get_parameters():
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description='Plotting aggregate metrics for Enrichment',
                                                             formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--traits', required=True, help="List of traits to plot enrichment scores")
    parser.add_argument('--predictor_of_choice', default="ABC", help="predictor of choice to compare with other predictors, if applicable")
    parser.add_argument('--data_outdir', required=True, help="Directory holding enrichment files")
    parser.add_argument('--cellTypeTable', default="/oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/ABC-GWAS/data/CellTypes.Annotated.ABCPaper.txt", help="Table with annotations of cell types, with columns 'CellType', 'Categorical.*', 'Binary.*' for plotting enrichments")
    parser.add_argument('--outdir', required=True, help="Directory to save plots in")
    parser.add_argument('--predictors', required=True, help="List of predictors to compare")
    return parser

def get_predict_argument_parser():
    parser = get_parameters()
    return parser

def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()
    trait = pd.read_csv(args.traits, sep="\t", header=None)
    predictor = pd.read_csv(args.predictors, sep="\t", header=None)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # plot single enrichment plots 
    # plot_single_enrichment_plots(list(trait[0]), list(predictor[0]), args)
    if len(predictor[0]) > 1:
        # plot single comparison barplots 
        plot_single_comparison_barplot(list(trait[0]), list(predictor[0]), args)

        # plot_aggregate_cdf
        plot_aggregate_cdf(list(trait[0]), list(predictor[0]), args)

def plot_single_enrichment_plots(traits, predictors, args):
    for trait in traits:
        for pred in predictors:
            data_outdir = os.path.join(args.data_outdir, pred, trait)
            os.system("Rscript PlotCellTypeEnrichment.R -o {} --cellTypes {} --cellTypeEnrichments {}/enrichment/Enrichment.CellType.vsScore.{}.tsv --codeDir /oak/stanford/groups/akundaje/kmualim/github/ABC-Max-pipeline/ --trait {}".format(args.outdir, args.cellTypeTable, data_outdir, trait, trait))

def plot_single_comparison_barplot(traits, predictors, args):
    for trait in traits:
        concat = None
        data = pd.read_csv("{}/{}/{}/enrichment/Enrichment.CellType.vsScore.{}.tsv".format(args.data_outdir, args.predictor_of_choice, trait, trait), sep="\t")
        data['predictions'] = args.predictor_of_choice
        if concat is None:
            concat = data
        else:
            concat = pd.concat([concat, data])
        label="other predictions"
        for pred in predictors:
            if pred != args.predictor_of_choice:
                data = pd.read_csv("{}/{}/{}/enrichment/Enrichment.CellType.vsScore.{}.tsv".format(args.data_outdir, pred, trait, trait), sep="\t")
                data['predictions'] = label
                concat = pd.concat([concat, data])
    ax = sns.boxplot(x="predictions", y="enrichment", data=concat)
    ax.set_xticklabels(ax.get_xticklabels(),rotation=30)
    plt.gcf().set_size_inches(10, 10)
    plt.title("Enrichment of {} variants of {} and other predictors".format(args.predictor_of_choice, trait))
    if not os.path.exists(os.path.join(args.outdir, trait)):
        os.makedirs(os.path.join(args.outdir, trait))
    plt.savefig("{}/{}/{}_enrichment_barplots_{}_against_all_predictions.pdf".format(args.outdir, trait, trait, args.predictor_of_choice), format='pdf')
    plt.clf()

def plot_aggregate_cdf(traits, predictors, args):
    for trait in traits:
        concat = None
        for pred in predictors:
            try:
                data = pd.read_csv("{}/{}/{}/enrichment/Enrichment.CellType.vsScore.{}.tsv".format(args.data_outdir, pred, trait, trait), sep="\t")
                data['predictions'] = pred
                if concat is None:
                    concat = data
                else:
                    concat = pd.concat([concat, data])
            except:
                print("PRED: {}".format(pred))
                print("trait: {}".format(trait))
    sns.ecdfplot(data=concat, x="enrichment", hue='predictions', stat='proportion')    
    plt.ylabel("Cumulative fraction")
    plt.title("Enrichment of {} variants intersected with {} Enhancers".format(trait, pred))
    plt.gcf().set_size_inches(15, 10)
    if not os.path.exists(os.path.join(args.outdir, trait)):
        os.makedirs(os.path.join(args.outdir, trait))
    plt.savefig("{}/{}/{}_across_all_predictions.pdf".format(args.outdir, trait, trait), format='pdf')
    plt.clf()
    ax = sns.boxplot(x="predictions", y="enrichment", data=concat)
    ax.set_xticklabels(ax.get_xticklabels(),rotation=30)
    plt.gcf().set_size_inches(10, 10)
    plt.title("Enrichment of {} variants intersected with {} Enhancers".format(trait, args.predictor_of_choice))
    plt.savefig("{}/{}/{}_barplot_enrichment_across_all_predictions.pdf".format(args.outdir, trait, trait), format='pdf')
    plt.clf()

if __name__=="__main__":
    main()
