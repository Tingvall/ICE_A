#!/usr/bin/env python

### PROCESS 12 MULTIPLE: UPSET PLOT FOR FACTOR BINDING IN ANCHOR POINTS - ALL INTERACTIONS WITH FACTOR AND INTERACTIONS FILTERES FOR GENES ###
import pandas as pd
import numpy as np
from upsetplot import plot
import matplotlib.pyplot as plt
import argparse

# PARSE ARGUMENTS
Description = 'UPSET PLOT FOR FACTOR BINDING IN ANCHOR POINTS'
Epilog = """Example usage: upset_plot.py <UPSET_PROMOTER> <UPSET_DISTAL> --prefix <PREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# Arguments
argParser.add_argument('UPSET_PROMOTER', help="Factor overlap in promoter regions.")
argParser.add_argument('UPSET_DISTAL', help="Factor overlap in distal regions.")
argParser.add_argument('DISTAL_PROMOTER', help="List of Distal-Promoter interactions for circos plot")
argParser.add_argument('--upset_promoter_g', dest='UPSET_PROMOTER_G', help="Factor overlap in promoter regions filtered by genes.")
argParser.add_argument('--upset_distal_g', dest='UPSET_DISTAL_G', help="Factor overlap in distal regions filtered by genes.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--circos_plot', dest="CIRCOS_PLOT", help="Specifies if Circos plot of peak overlap will be created.", choices=['true', 'false'])
argParser.add_argument('--filter_genes', dest="FILTER_GENES", help="Specifies if additional plot (Upset and/or Circos plots) should be created based on interactions filtered by provided gene list (default: false). This option requires that a gene list is provided with the argument --genes.", choices=['true', 'false'])
argParser.add_argument('--complete', dest="COMPLETE", help="If set to true, all available processes for the selected mode and provided inputs are run.", choices=['true', 'false'])

args = argParser.parse_args()

# DEFINE FUNCTION
def upset_plot(upset_promoter, upset_distal, distal_promoter, upset_promoter_g, upset_distal_g, prefix, circos_plot, filter_genes, complete):

  ### Loading and organizing data
  upset_promoter = pd.read_table(upset_promoter)
  upset_distal = pd.read_table(upset_distal)

  factor = pd.concat([upset_promoter, upset_distal])['Source'].unique()

  ## Upset PLOTS
  # Promoter all
  for f in factor:
      upset_promoter[f] = np.where(upset_promoter['Source'] == f, True, False)
  upset_promoter_anchor = upset_promoter.iloc[:,np.r_[2,4:len(upset_promoter.columns)]]
  upset_promoter_anchor = upset_promoter_anchor.groupby(upset_promoter_anchor.columns[0]).max()
  upset_promoter_anchor_group = upset_promoter_anchor.groupby(list(factor)).size().to_frame('size')
  plot(upset_promoter_anchor_group['size'], sort_by="cardinality")
  plt.savefig('Upset_plot_Promoter_all.pdf')


  # Distal all
  for f in factor:
      upset_distal[f] = np.where(upset_distal['Source'] == f, True, False)
  upset_distal_anchor = upset_distal.iloc[:,np.r_[2,4:len(upset_distal.columns)]]
  upset_distal_anchor = upset_distal_anchor.groupby(upset_distal_anchor.columns[0]).max()
  upset_distal_anchor_group = upset_distal_anchor.groupby(list(factor)).size().to_frame('size')
  plot(upset_distal_anchor_group['size'], sort_by="cardinality")
  plt.savefig('Upset_plot_Distal_all.pdf')

  if filter_genes == 'true':
    upset_promoter_g = pd.read_table(upset_promoter_g)
    upset_distal_g = pd.read_table(upset_distal_g)

    # Promoter genes
    for f in factor:
        upset_promoter_g[f] = np.where(upset_promoter_g['Source'] == f, True, False)
    upset_promoter_g_anchor = upset_promoter_g.iloc[:,np.r_[2,4:len(upset_promoter_g.columns)]]
    upset_promoter_g_anchor = upset_promoter_g_anchor.groupby(upset_promoter_g_anchor.columns[0]).max()
    upset_promoter_g_anchor_group = upset_promoter_g_anchor.groupby(list(factor)).size().to_frame('size')
    plot(upset_promoter_g_anchor_group['size'], sort_by="cardinality")
    plt.savefig('Upset_plot_Promoter_genelist.pdf')

    # Distal genes
    for f in factor:
        upset_distal_g[f] = np.where(upset_distal_g['Source'] == f, True, False)
    upset_distal_g_anchor = upset_distal_g.iloc[:,np.r_[2,4:len(upset_distal_g.columns)]]
    upset_distal_g_anchor = upset_distal_g_anchor.groupby(upset_distal_g_anchor.columns[0]).max()
    upset_distal_g_anchor_group = upset_distal_g_anchor.groupby(list(factor)).size().to_frame('size')
    plot(upset_distal_g_anchor_group['size'], sort_by="cardinality")
    plt.savefig('Upset_plot_Distal_genelist.pdf')

  ### Preperations for circos PLOTS
  if (complete == 'true' or circos_plot == 'true'):
    upset_promoter = upset_promoter.iloc[:,np.r_[0,4:len(upset_promoter.columns)]]
    upset_promoter = upset_promoter.groupby(upset_promoter.columns[0]).max()
    upset_promoter['promoter_cat'] = 'Promoter'+upset_promoter.eq(True).dot('_'+upset_promoter.columns)
    upset_distal = upset_distal.iloc[:,np.r_[0,4:len(upset_distal.columns)]]
    upset_distal = upset_distal.groupby(upset_distal.columns[0]).max()
    upset_distal['distal_cat'] = 'Distal'+upset_distal.eq(True).dot('_'+upset_distal.columns)

    circos_f = upset_promoter.merge(upset_distal, left_index=True, right_index=True, how = 'outer')
    circos_f.fillna(value={'promoter_cat': 'Promoter_NoBinding', 'distal_cat': 'Distal_NoBinding'}, inplace=True)
    circos_f.fillna(False,inplace=True)

    #Filter for interactions in distal-promoter list (promoter-promoter excluded)
    distal_promoter = pd.read_table(distal_promoter, index_col=0)
    circos_f = circos_f[circos_f.index.isin(distal_promoter.index)]
    circos_f = circos_f.groupby(list(circos_f.columns)).size().to_frame('size').reset_index()
    circos_f.to_csv('Circos_peaks_' + prefix + '_interactions.txt', index=False, sep='\t' )

    if filter_genes == 'true':
      upset_promoter_g = upset_promoter_g.iloc[:,np.r_[0,4:len(upset_promoter_g.columns)]]
      upset_promoter_g = upset_promoter_g.groupby(upset_promoter_g.columns[0]).max()
      upset_promoter_g['promoter_cat'] = 'Promoter'+upset_promoter_g.eq(True).dot('_'+upset_promoter_g.columns)
      upset_distal_g = upset_distal_g.iloc[:,np.r_[0,4:len(upset_distal_g.columns)]]
      upset_distal_g = upset_distal_g.groupby(upset_distal_g.columns[0]).max()
      upset_distal_g['distal_cat'] = 'Distal'+upset_distal_g.eq(True).dot('_'+upset_distal_g.columns)

      circos_g = upset_promoter_g.merge(upset_distal_g, left_index=True, right_index=True, how = 'outer')
      circos_g.fillna(value={'promoter_cat': 'promoter_NoBinding', 'distal_cat': 'distal_NoBinding'}, inplace=True)
      circos_g.fillna(False,inplace=True)
      circos_g = circos_g[circos_g.index.isin(distal_promoter.index)]
      circos_g = circos_g.groupby(list(circos_g.columns)).size().to_frame('size').reset_index()
      circos_g.to_csv('Circos_genes_' + prefix + '_interactions.txt', index=False, sep='\t' )

# RUN FUNCTION
upset_plot(upset_promoter=args.UPSET_PROMOTER, upset_distal=args.UPSET_DISTAL, distal_promoter=args.DISTAL_PROMOTER,upset_promoter_g=args.UPSET_PROMOTER_G, upset_distal_g=args.UPSET_DISTAL_G, prefix=args.PREFIX, circos_plot=args.CIRCOS_PLOT, filter_genes=args.FILTER_GENES, complete=args.COMPLETE)
