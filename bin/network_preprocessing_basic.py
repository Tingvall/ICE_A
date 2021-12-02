#!/usr/bin/env python

### PROCESS 10 BASIC: CREATES NETWORK FILES FOR CYTOSCAPE ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'PREPROCESSING OF ANNOTATED INTERATION FOR CYTOSCAPE NETWORK VISUALIZATION'
Epilog = """Example usage: network_preprocessing_basic.py <INTERACTIONS_ANNO_AGG> <INTERACTIONS_ANNO> <GENES> --prefix <PREFIX> --sample <SAMPLE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# General arguments
argParser.add_argument('INTERACTIONS_ANNO_AGG', help="Annotated and aggregated interactions.")
argParser.add_argument('INTERACTIONS_ANNO', help="Annotated, not aggregated interactions")
argParser.add_argument('GENES', help="Text file specifying genes for filtering.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--sample', dest="SAMPLE", help="Name of sample.")
argParser.add_argument('--network_mode', dest="NETWORK_MODE", help="Defines mode network. Options are all (all interaction in the 2D-bed file), factor (all interaction with at least on peak overlap either anchor point) or genes (interactions associates with a gene list, provided by --genes)." , choices=['all', 'factor', 'genes'])
argParser.add_argument('--promoter_promoter', dest="PROMOTER_PROMOTER", help="If set to true, promoter-promoter interactions included in network (default: false).", choices=['true', 'false'])
argParser.add_argument('--complete', dest="COMPLETE", help="If set to true, all available processes for the selected mode and provided inputs are run.", choices=['true', 'false'])

args = argParser.parse_args()

# DEFINE FUNCTION
def network_preprocessing_basic(interactions_annotated, interactions_annotated_not_aggregated, genes, prefix, sample, network_mode, promoter_promoter, complete):

  #Loading input file
  anchors_peaks_anno = pd.read_table(interactions_annotated_not_aggregated, index_col=0)
  interactions_anno = pd.read_table(interactions_annotated, index_col=0)

  # Aggregating interaction file to only incude one row per interaction
  interactions_anno = interactions_anno.iloc[:,np.r_[0:5,7:15, 17:len(anchors_peaks_anno.columns)]]
  interactions_anno['Anchor1'] = interactions_anno["chr1"].map(str) +':'+ (interactions_anno["s1"]).map(str) +'-'+ interactions_anno["e1"].map(str)
  interactions_anno['Anchor2'] = interactions_anno["chr2"].map(str) +':'+ (interactions_anno["s2"]).map(str) +'-'+ interactions_anno["e2"].map(str)
  interactions_anno = pd.concat([interactions_anno['Anchor1'], interactions_anno.iloc[:,3:5], interactions_anno['Anchor2'],interactions_anno.iloc[:,8:(len(interactions_anno.columns)-2)]], axis=1)

  ### Creating edge table for cytoscape
  #Factor-Interaction
  Factor_Interaction = anchors_peaks_anno.copy(deep=True)
  Factor_Interaction.loc[Factor_Interaction.Overlap_1 == 1, 'Peak1'] = sample
  Factor_Interaction.loc[Factor_Interaction.Overlap_2 == 1, 'Peak2'] = "sample
  Factor_Interaction = Factor_Interaction[['chr1', 's1', 'e1','Gene_Name_1', 'Peak1','Peak1_ID','Peak1_score', 'chr2', 's2', 'e2',  'Gene_Name_2','Peak2','Peak2_ID','Peak2_score', 'TSS_1', 'TSS_2']]
  Factor_Interaction['Anchor1'] = Factor_Interaction['chr1'].map(str) +':'+ (Factor_Interaction['s1']).map(str) +'-'+ Factor_Interaction['e1'].map(str)
  Factor_Interaction['Anchor2'] = Factor_Interaction['chr2'].map(str) +':'+ (Factor_Interaction['s2']).map(str) +'-'+ Factor_Interaction['e2'].map(str)
  Factor_Interaction = Factor_Interaction.dropna(subset=['Peak1', 'Peak2'], thresh=1)

  #Factor-Distal
  Factor_Distal_1 = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 0) & (Factor_Interaction['TSS_2'] == 1), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Distal_1.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal_2 = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['TSS_2'] == 0), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Distal_2.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal = Factor_Distal_1.append(Factor_Distal_2)
  Factor_Distal['Edge_type'] = 'Factor-Distal'

  #Factor-Promoter
  Factor_Promoter_1 = Factor_Interaction.loc[Factor_Interaction['TSS_1'] == 1, ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Promoter_1.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_2 = Factor_Interaction.loc[Factor_Interaction['TSS_2'] == 1, ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Promoter_2.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter = Factor_Promoter_1.append(Factor_Promoter_2)
  Factor_Promoter['Edge_type'] = 'Factor-Promoter'

  #Distal-Promoter
  DP_1 = interactions_anno.loc[(interactions_anno['TSS_1'] == 0) & (interactions_anno['TSS_2'] == 1), ['Anchor1','Anchor2', 'Q-Value_Bias']]
  DP_1.columns = ['Source', 'Target', 'Edge_score']
  DP_2 = interactions_anno.loc[(interactions_anno['TSS_1'] == 1) & (interactions_anno['TSS_2'] == 0), ['Anchor2',  'Anchor1', 'Q-Value_Bias']]
  DP_2.columns = ['Source', 'Target', 'Edge_score']
  Distal_Promoter = DP_1.append(DP_2)
  Distal_Promoter['Edge_type'] = 'Distal-Promoter'
  Distal_Promoter['Edge_score'] = - np.log10(Distal_Promoter['Edge_score'])

  #Promoter-Promoter
  Promoter_Promoter = interactions_anno.loc[(interactions_anno['TSS_1']==1) & (interactions_anno['TSS_2']==1),:][['Anchor1', 'Anchor2', 'Q-Value_Bias']]
  Promoter_Promoter['Edge_type'] = 'Promoter-Promoter'
  Promoter_Promoter.columns = ['Source', 'Target', 'Edge_score', 'Edge_type']
  Promoter_Promoter['Edge_score'] = - np.log10(Promoter_Promoter['Edge_score'])

  #Promoter-Gene
  Promoter_Gene_1 = Factor_Interaction.loc[Factor_Interaction['TSS_1'] == 1, ['Anchor1', 'Gene_Name_1']].dropna(subset=['Gene_Name_1']).drop_duplicates()
  Promoter_Gene_1.columns = ['Source', 'Target']
  Promoter_Gene_2 = Factor_Interaction.loc[Factor_Interaction['TSS_2'] == 1, ['Anchor2',  'Gene_Name_2']].dropna(subset=['Gene_Name_2']).drop_duplicates()
  Promoter_Gene_2.columns = ['Source', 'Target']
  Promoter_Gene = Promoter_Gene_1.append(Promoter_Gene_2)
  Promoter_Gene['Edge_score'], Promoter_Gene['Edge_type'] = [1, 'Promoter-Gene']

  # Filtering of edges based on network mode
  if network_mode == "factor":
    #Filter edges based on factor
    Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target'])]
    Promoter_Promoter_filt_f = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter['Target'])]
    Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Target'])]

  elif network_mode == "genes":
    #Filter edges based on gene
    genes = pd.read_table(genes, header=None)
    Promoter_Gene_filt_g = Promoter_Gene[Promoter_Gene['Target'].isin(genes.iloc[:,0])]
    Promoter_Promoter_filt_g = Promoter_Promoter[Promoter_Promoter['Source'].isin(Promoter_Gene_filt_g['Source']) | Promoter_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Distal_Promoter_filt_g = Distal_Promoter[Distal_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Factor_Promoter_filt_g = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Factor_Distal_filt_g = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g['Source'])]

  ### Creating edge table for cytoscape
  if network_mode == "all":
    if promoter_promoter =="true":
      Edges = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Promoter, Promoter_Gene]).drop_duplicates()
    else:
      Edges = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Gene]).drop_duplicates()

  elif network_mode == "factor":
    if promoter_promoter =="true":
      Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
    else:
      Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()

  elif network_mode == "genes":
    if promoter_promoter =="true":
      Edges =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()
    else:
      Edges =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()

  Edges.to_csv('Network_Edges_' + prefix + '_interactions.txt', index=False, sep='\t' )


  ### Creating node table for cytoscape
  if network_mode == "all":
    #Specifying node type for all nodes
    Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
    Nodes.columns=['Node']
    if promoter_promoter =="true":
      Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Distal_Promoter['Target']) | Nodes['Node'].isin(Promoter_Promoter['Source']) | Nodes['Node'].isin(Promoter_Promoter['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))
    else:
      Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                     (np.where(Nodes['Node'].isin(Distal_Promoter['Target']) | Nodes['Node'].isin(Promoter_Promoter['Source']) , 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))

  elif network_mode == "factor":
    # Specifying node type for all nodes that are associated with factor binding
    Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
    Nodes.columns=['Node']
    if promoter_promoter =="true":
      Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))
    else:
      Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                   (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Promoter_Promoter['Source']) , 'Promoter',
                                      (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))

  elif network_mode == "genes":
    # Specifying node type for all nodes that are associated with selected genes
    Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
    Nodes.columns=['Node']
    if promoter_promoter =="true":
      Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_g['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_g['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))
    else:
      Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                     (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Target'])| Nodes['Node'].isin(Promoter_Promoter_filt_g['Source']) , 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))

  Nodes.to_csv('Network_Nodes_' + prefix + '_interactions.txt', index=False, sep='\t' )


# RUN FUNCTION
network_preprocessing_basic(interactions_annotated=args.INTERACTIONS_ANNO_AGG, interactions_annotated_not_aggregated=args.INTERACTIONS_ANNO, genes=args.GENES, prefix=args.PREFIX, sample=args.SAMPLE, network_mode=args.NETWORK_MODE, promoter_promoter=args.PROMOTER_PROMOTER, complete=args.COMPLETE)
