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
    parser.add_argument('--data_outdir', default=".", help="Directory holding enrichment files")
    parser.add_argument('--outdir', required=True, help="Directory to save plots in")
    parser.add_argument('--predictors', nargs='+', required=True, help="List of predictors to compare")
    return parser

def get_predict_argument_parser():
    parser = get_parameters()
    return parser

def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()
    trait = args.traits 
    predictor = args.predictors 
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if len(predictor)> 1:
        # plot single comparison barplots 
        plot_single_comparison_barplot(trait, predictor, args)
        # plot_aggregate_cdf
        plot_aggregate_cdf(trait, predictor, args)

def plot_single_comparison_barplot(trait, predictors, args):
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

def plot_aggregate_cdf(trait, predictors, args):
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
