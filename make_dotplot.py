#!/usr/bin/env python3
#Copyright © 2020 Giulia I. Corsi
import argparse as ap
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import seaborn as sb
import os
import numpy as np

sb.set_style("whitegrid")

PALETTE = {'start':0.2, 'end':0.8} # used to smooth the palette and avoid very dark or very light colors that are at the edges of palettes
def parse_args():
    parser = ap.ArgumentParser()
    parser.add_argument('-e', '--enrichment_results', help='Path to table with enrichment results, in csv format ',
                        type=str, required=True)
    parser.add_argument('-o', help='Path to output folder', type=str, required = True)
    parser.add_argument('--GSEA', help='Source of the GSEA analysis', choices=['StringApp', 'Other'],
                        default='StringApp')
    parser.add_argument('-c_n', help='Name of the column containing the number of background genes', type=str, required=False)
    parser.add_argument('-c_nb', help='Name of the column containing the number of genes assigned to each term',
                        type=str, required=False)
    parser.add_argument('-c_c', help='Name of the column containing a category, used to group enriched terms', type=str, required=False)
    parser.add_argument('-c_d', help='Name of the column containing the term description', type=str, required=False)
    parser.add_argument('-c_f', help='Name of the column containing the FDR value', type=str, required=False)
    parser.add_argument('-n', help='Maximum number of terms to plot', type=int, default=20)
    parser.add_argument('--ratio_min',
                        help='Consider only terms for which the gene-ratio term is above the given threshold',
                        type=float, default=0.00)
    parser.add_argument('--palette', help='Color palette used', default='viridis')
    parser.add_argument('--n_bins', help = 'Number of bins for grouping gene counts', default=4, type=int, choices = [2,3,4])
    parser.add_argument('--groups', help = 'List of categories to be grouped together', default=['KEGG Pathways','Reactome Pathways'], nargs='+')
    parser.add_argument('--increase_height', help = 'Integer, increases the height of the plot', type = int, default=0)
    parser.add_argument('--increase_width', help = 'Integer, increases the width of the plot', type = int, default=0)
    args = parser.parse_args()

    if args.GSEA == 'Other':
        if not args.c_n or not args.c_nb or not args.c_d:  # c_c is optional
            parser.error(
                'If the GSEA source is unknown, names of the columns to be used in the table have to be provided (args.c_*)')
    if args.GSEA == 'StringApp':
        args.c_n = '# genes'
        args.c_nb = '# background genes'
        args.c_d = 'description'
        args.c_c = 'category'
        args.c_f = 'FDR value'
    if not os.path.exists(args.o):
        os.mkdir(args.o)
    return args


def parse_data(path_enrichment: str, c_n: str, c_nb: str, ratio_min:float, c_f: str, c_d: str, c_c: str):
    df = pd.read_csv(path_enrichment)
    df['ratio'] = df[c_n] / df[c_nb]
    if ratio_min:
        df = df[df['ratio'] > ratio_min]
    df = df[[c_n, c_f, c_d, c_c, 'ratio']]
    df.columns = ['Gene Count', 'FDR', 'Description', 'Category', 'Gene Ratio']
    return df


def dot_plot_enrichment(df: pd.DataFrame, title: str, out_folder: str, c_d: str, c_r: str, c_n: str, n_max: int,
                        c_f: str, palette:str, n_bins:int, height_add:int, width_add :int):
    def get_exponent(x):
        str_x = '%.2E' % x
        str_x = str_x[str_x.find('E'):]
        return str_x

    def get_color_bins(ser: pd.Series):
        min_exp = int(get_exponent(ser.min())[1:])
        max_exp = int(get_exponent(ser.max())[1:])
        exponents = ['E%.2i'%x for x in range(min_exp, max_exp+1,1) ]
        color_bins = {x: id for id, x in enumerate(exponents)}
        return color_bins


    df.reindex()
    df.sort_values(c_f)
    df = df.iloc[:n_max].copy()
    SIZES = {0:20,1:120,2:250,3:400}
    width_xfigure = width_add+5+df[c_d].apply(lambda x: len(x)).max()/15
    height_yfigure = height_add+len(df)/3.5
    plt.figure(figsize=(width_xfigure, height_yfigure))
    df['count_bins'] = (pd.cut(df[c_n], min(n_bins, len(df[c_n].unique()))).apply(lambda x: x.right)).astype(int)
    if len(df)>=5: #len(df['count_bins'].unique()) == n_bins and
        size_bins =  {x: SIZES[i] for i,x in enumerate(sorted(df['count_bins'].unique()))}
        df['count_bins_sizes'] = df['count_bins'].apply(lambda x:size_bins[x])
        print(df)
        print(size_bins)
        color_bins = get_color_bins(df[c_f])
        df[c_f] = df[c_f].apply(lambda x: get_exponent(x))
        #color_bins = {x:id for id, x in enumerate(df[c_f].unique())}
        df['color_bins']=df[c_f].apply(lambda x: color_bins[x])
        newcmp = ListedColormap(cm.get_cmap(palette, 256)(np.linspace(PALETTE['start'], PALETTE['end'], 256)))
        dot_plot = plt.scatter(x=df[c_r], y=range(len(df)), s=df['count_bins_sizes'], c=df['color_bins'], cmap=newcmp)
        plt.gca().invert_yaxis()
        plt.yticks(range(len(df[c_d])),df[c_d], fontsize=11)
        markers=[]
        lst_size_bins_counts = [0]
        lst_size_bins_counts.extend([x for x in sorted(size_bins.keys())])
        for i, bin in enumerate(sorted(size_bins.keys())):
            markers.append(plt.scatter([], [], c='grey', s=size_bins[bin], label=str(lst_size_bins_counts[i])+'-'+str(bin)))
        plt.legend(title = 'Gene counts', handles=markers,labelspacing=1, borderpad=1, fontsize=11)
        cbar = plt.colorbar(dot_plot, label = c_f, ticks = [0,int(len(color_bins)/2),len(color_bins)-1])
        cbar.set_ticklabels([df[c_f].unique()[0],df[c_f].unique()[int(len(df[c_f].unique())/2)],df[c_f].unique()[-1]])
        plt.xlabel('Gene Ratio', fontsize=11)
        plt.tight_layout()
        plt.savefig(os.path.join(out_folder,'%s_enrichment.svg'%title), format='svg')
    else:
        print('Figure for category %s is not produced because the number of terms is smaller than the number of gene count bins or because the number of terms is <5'%title)


if __name__ == '__main__':
    args = parse_args()
    df = parse_data(path_enrichment=args.enrichment_results, c_n=args.c_n, c_nb=args.c_nb, ratio_min=args.ratio_min,
                    c_f=args.c_f, c_d=args.c_d, c_c=args.c_c)
    if args.c_c:
        print(args.groups)
        if args.groups:
            category = '-'.join(args.groups)
            print(category)
            df['Category']=df['Category'].apply(lambda x: category if x in args.groups else x)
        for category, df_category in df.groupby('Category'):
            print(category)
            dot_plot_enrichment(df=df_category.copy(), title=category.replace(' ','_'), out_folder=args.o, c_d='Description',
                                c_r='Gene Ratio', c_n='Gene Count', n_max=args.n, c_f='FDR', palette=args.palette, n_bins=args.n_bins, height_add=args.increase_height, width_add = args.increase_width)
    else:
        dot_plot_enrichment(df=df, title='All', out_folder=args.o, c_d='Description', c_r='Gene Ratio',
                            c_n='Gene Count', n_max=args.n, c_f='FDR', palette=args.palette, n_bins= args.n_bins, height_add=args.increase_height, width_add = args.increase_width)
