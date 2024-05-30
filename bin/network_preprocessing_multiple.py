#!/usr/bin/env python

### PROCESS 10 MULTIPLE: CREATES NETWORK FILES FOR CYTOSCAPE ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'PREPROCESSING OF ANNOTATED INTERATION FOR CYTOSCAPE NETWORK VISUALIZATION IN MULTIPLE MODE'
Epilog = """Example usage: network_preprocessing_multiple.py <INTERACTIONS_ANNO_AGG> <INTERACTIONS_ANNO> --genes <GENES> --prefix <PREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# General arguments
argParser.add_argument('INTERACTIONS_ANNO_AGG', help="Annotated and aggregated interactions.")
argParser.add_argument('INTERACTIONS_ANNO', help="Annotated, not aggregated interactions")
argParser.add_argument('--genes', dest='GENES', help="Text file specifying genes for filtering.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--network_mode', dest="NETWORK_MODE", help="Defines mode network. Options are a (all interaction in the 2D-bed file), f (all interaction with at least on peak overlap either anchor point),g (interactions associates with a gene list, provided by --genes) or fg (combiantion of option f & g)." , choices=['a', 'f', 'g', 'fg'])
argParser.add_argument('--promoter_promoter', dest="PROMOTER_PROMOTER", help="If set to true, promoter-promoter interactions included in network (default: false).", choices=['true', 'false'])
argParser.add_argument('--complete', dest="COMPLETE", help="If set to true, all available processes for the selected mode and provided inputs are run.", choices=['true', 'false'])
argParser.add_argument('--network_distal_only', dest="NETWORK_DISTAL_ONLY", help="If true, only distal factor binding are shown in the netork.", choices=['true', 'false'])

# Multiple mode specific arguments
argParser.add_argument('--upset_plot', dest="UPSET_PLOT", help="Specifies if Upset plot of peak overlap will be created.", choices=['true', 'false'])
argParser.add_argument('--circos_plot', dest="CIRCOS_PLOT", help="Specifies if Circos plot of peak overlap will be created.", choices=['true', 'false'])
argParser.add_argument('--filter_genes', dest="FILTER_GENES", help="Specifies if additional plot (Upset and/or Circos plots) should be created based on interactions filtered by provided gene list (default: false). This option requires that a gene list is provided with the argument --genes.", choices=['true', 'false'])
argParser.add_argument('--circos_use_promoters', dest="CIRCOS_USE_PROMOTERS", help="Specifies if TF overlap in promoters (defined based on promoter_start/end) should be used in circos plot in multiple mode when regions are specified.", choices=['true', 'false'])

args = argParser.parse_args()

# DEFINE FUNCTION
def network_preprocessing_multiple(interactions_annotated, interactions_annotated_not_aggregated, genes, prefix, network_mode, promoter_promoter, upset_plot, circos_plot, filter_genes, complete, network_distal_only, circos_use_promoters):

    #Loading input file
    anchors_peaks_anno = pd.read_table(interactions_annotated_not_aggregated, index_col=0)
    anchors_peaks_anno = anchors_peaks_anno.assign(Peak1_score = pd.to_numeric(anchors_peaks_anno['Peak1_score']))
    anchors_peaks_anno = anchors_peaks_anno.assign(Peak2_score = pd.to_numeric(anchors_peaks_anno['Peak2_score']))
    anchors_peaks_anno = anchors_peaks_anno.loc[((anchors_peaks_anno['Is_Promoter_1'] == 0) & (anchors_peaks_anno['Peak1_score'] == 1) & (anchors_peaks_anno['Is_Promoter_2'] == 1) & (anchors_peaks_anno['Peak2_score'] == 0)) | ((anchors_peaks_anno['Is_Promoter_1'] == 1) & (anchors_peaks_anno['Peak1_score'] == 0) & (anchors_peaks_anno['Is_Promoter_2'] == 0) & (anchors_peaks_anno['Peak2_score'] == 1)),:]

    interactions_anno = pd.read_table(interactions_annotated, index_col=0)
    interactions_anno =  interactions_anno[interactions_anno.index.isin(anchors_peaks_anno.index)]

    # Aggregating interaction file to only incude one row per interaction
    interactions_anno = interactions_anno.iloc[:,np.r_[0:5,6,10:15, 16, 20:len(anchors_peaks_anno.columns)]]
    interactions_anno['Anchor1'] = interactions_anno["chr1"].map(str) +':'+ (interactions_anno["s1"]).map(str) +'-'+ interactions_anno["e1"].map(str)
    interactions_anno['Anchor2'] = interactions_anno["chr2"].map(str) +':'+ (interactions_anno["s2"]).map(str) +'-'+ interactions_anno["e2"].map(str)
    interactions_anno = pd.concat([interactions_anno['Anchor1'], interactions_anno.iloc[:,3:6], interactions_anno['Anchor2'],interactions_anno.iloc[:,9:(len(interactions_anno.columns)-2)]], axis=1)

    # Factor-Interaction
    Factor_Interaction_all = anchors_peaks_anno[['chr1', 's1', 'e1','Gene_Name_1', 'Peak1','Peak1_ID', 'Peak1_score', 'chr2', 's2', 'e2',  'Gene_Name_2','Peak2','Peak2_ID','Peak2_score', 'Is_Promoter_1', 'Is_Promoter_2']]
    Factor_Interaction_all['Anchor1'] = Factor_Interaction_all['chr1'].map(str) +':'+ (Factor_Interaction_all['s1']).map(str) +'-'+ Factor_Interaction_all['e1'].map(str)
    Factor_Interaction_all['Anchor2'] = Factor_Interaction_all['chr2'].map(str) +':'+ (Factor_Interaction_all['s2']).map(str) +'-'+ Factor_Interaction_all['e2'].map(str)
    Factor_Interaction = Factor_Interaction_all.dropna(subset=['Peak1', 'Peak2'], thresh=1)

    #Factor-Distal
    if (circos_use_promoters == "true"):
        Factor_Distal_1 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 0) & (Factor_Interaction['Peak1_score'] == 1) & (Factor_Interaction['Is_Promoter_2'] == 1) & (Factor_Interaction['Peak2_score'] == 0), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).reset_index().drop_duplicates().set_index('Interaction')
        Factor_Distal_2 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['Peak1_score'] == 0) & (Factor_Interaction['Is_Promoter_2'] == 0) & (Factor_Interaction['Peak2_score'] == 1), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).reset_index().drop_duplicates().set_index('Interaction')
    else:
        Factor_Distal_1 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 0) & (Factor_Interaction['Is_Promoter_2'] == 1), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).reset_index().drop_duplicates().set_index('Interaction')
        Factor_Distal_2 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['Is_Promoter_2'] == 0), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).reset_index().drop_duplicates().set_index('Interaction')
    Factor_Distal_1.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_2.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal = Factor_Distal_1.append(Factor_Distal_2)
    Factor_Distal['Edge_type'] = 'Factor-Distal'

    #Factor-Promoter
    if (circos_use_promoters == "true"):
        Factor_Promoter_1 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['Peak1_score'] == 0) & (Factor_Interaction['Is_Promoter_2'] == 0) & (Factor_Interaction['Peak2_score'] == 1), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).reset_index().drop_duplicates().set_index('Interaction')
        Factor_Promoter_2 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 0) & (Factor_Interaction['Peak1_score'] == 1) & (Factor_Interaction['Is_Promoter_2'] == 1) & (Factor_Interaction['Peak2_score'] == 0), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).reset_index().drop_duplicates().set_index('Interaction')
    else:
        Factor_Promoter_1 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_1'] == 1, ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).reset_index().drop_duplicates().set_index('Interaction')
        Factor_Promoter_2 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_2'] == 1, ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).reset_index().drop_duplicates().set_index('Interaction')
    Factor_Promoter_1.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_2.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter = Factor_Promoter_1.append(Factor_Promoter_2)
    Factor_Promoter['Edge_type'] = 'Factor-Promoter'

    #Distal-Promoter
    DP_1 = interactions_anno.loc[(interactions_anno['Is_Promoter_1'] == 0) & (interactions_anno['Is_Promoter_2'] == 1), ['Anchor1','Anchor2', 'Interaction_score']]
    DP_2 = interactions_anno.loc[(interactions_anno['Is_Promoter_1'] == 1) & (interactions_anno['Is_Promoter_2'] == 0), ['Anchor2',  'Anchor1', 'Interaction_score']]
    DP_1.columns = ['Source', 'Target', 'Edge_score']
    DP_2.columns = ['Source', 'Target', 'Edge_score']
    Distal_Promoter = DP_1.append(DP_2)
    Distal_Promoter['Edge_type'] = 'Distal-Promoter'
    Distal_Promoter['Edge_score'] = - np.log10(Distal_Promoter['Edge_score'])
    Distal_Promoter.to_csv('Distal_promoter_for_circos.txt', index=True, sep='\t' )

    #Promoter-Promoter
    Promoter_Promoter = interactions_anno.loc[(interactions_anno['Is_Promoter_1']==1) & (interactions_anno['Is_Promoter_2']==1),:][['Anchor1', 'Anchor2', 'Interaction_score']]
    Promoter_Promoter['Edge_type'] = 'Promoter-Promoter'
    Promoter_Promoter.columns = ['Source', 'Target', 'Edge_score', 'Edge_type']
    Promoter_Promoter['Edge_score'] = - np.log10(Promoter_Promoter['Edge_score'])

    #Promoter-Gene
    Promoter_Gene_1 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_1'] == 1, ['Anchor1', 'Gene_Name_1']].dropna(subset=['Gene_Name_1']).reset_index().drop_duplicates().set_index('Interaction')
    Promoter_Gene_2 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_2'] == 1, ['Anchor2',  'Gene_Name_2']].dropna(subset=['Gene_Name_2']).reset_index().drop_duplicates().set_index('Interaction')
    Promoter_Gene_1.columns = ['Source', 'Target']
    Promoter_Gene_2.columns = ['Source', 'Target']
    Promoter_Gene = Promoter_Gene_1.append(Promoter_Gene_2)
    Promoter_Gene['Edge_score'], Promoter_Gene['Edge_type'] = [1, 'Promoter-Gene']

    # Filtering of edges based on network mode
    if (network_mode == "f" or network_mode == "fg"):
    #Filter edges based on factor
        Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target'])]
        if promoter_promoter =="true":
            Promoter_Promoter_filt_f = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter['Target'])]
            if network_distal_only=="true":
                Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Target'])]
            else:
                Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Target'])]
        else:
            if network_distal_only=="true":
                Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Target'])]
            else:
                Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Target'])]

        if network_mode == "fg":
            #Filter edges based on gene
            genes = pd.read_table(genes, header=None)
            if promoter_promoter =="true":
                if network_distal_only=="true":
                    Promoter_Gene_filt_fg = Promoter_Gene_filt_f[(Promoter_Gene_filt_f['Target'].isin(genes.iloc[:,0])) & (Promoter_Gene_filt_f['Source'].isin(Distal_Promoter_filt_f['Target']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Source']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Target']))]
                else:
                    Promoter_Gene_filt_fg = Promoter_Gene_filt_f[(Promoter_Gene_filt_f['Target'].isin(genes.iloc[:,0])) & (Promoter_Gene_filt_f['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene_filt_f['Source'].isin(Distal_Promoter_filt_f['Target']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Source']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Target']))]
            else:
                if network_distal_only=="true":
                    Promoter_Gene_filt_fg = Promoter_Gene_filt_f[(Promoter_Gene_filt_f['Target'].isin(genes.iloc[:,0])) & (Promoter_Gene_filt_f['Source'].isin(Distal_Promoter_filt_f['Target']))]
                else:
                    Promoter_Gene_filt_fg = Promoter_Gene_filt_f[(Promoter_Gene_filt_f['Target'].isin(genes.iloc[:,0])) & (Promoter_Gene_filt_f['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene_filt_f['Source'].isin(Distal_Promoter_filt_f['Target']))]
            Distal_Promoter_filt_fg = Distal_Promoter_filt_f[Distal_Promoter_filt_f['Target'].isin(Promoter_Gene_filt_fg['Source'])]
            Factor_Distal_filt_fg = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_fg['Source'])]
            Factor_Promoter_filt_fg = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_fg['Source'])]
            if promoter_promoter =="true":
                Promoter_Promoter_filt_fg = Promoter_Promoter_filt_f[(Promoter_Promoter_filt_f['Source'].isin(Promoter_Gene_filt_fg['Source']) & Promoter_Promoter_filt_f['Target'].isin(Factor_Promoter_filt_fg['Source'])) | (Promoter_Promoter_filt_f['Target'].isin(Promoter_Gene_filt_fg['Source']) & Promoter_Promoter_filt_f['Source'].isin(Factor_Promoter_filt_fg['Source']))]
                if network_distal_only=="true":
                    Promoter_Gene_filt_fg = Promoter_Gene_filt_fg[Promoter_Gene_filt_fg['Source'].isin(Distal_Promoter_filt_fg['Target']) | Promoter_Gene_filt_fg['Source'].isin(Promoter_Promoter_filt_fg['Source']) | Promoter_Gene_filt_fg['Source'].isin(Promoter_Promoter_filt_fg['Target'])]
                else:
                    Promoter_Gene_filt_fg = Promoter_Gene_filt_fg[Promoter_Gene_filt_fg['Source'].isin(Factor_Promoter_filt_fg['Target']) | Promoter_Gene_filt_fg['Source'].isin(Distal_Promoter_filt_fg['Target']) | Promoter_Gene_filt_fg['Source'].isin(Promoter_Promoter_filt_fg['Source']) | Promoter_Gene_filt_fg['Source'].isin(Promoter_Promoter_filt_fg['Target'])]

    elif network_mode == "g":
        #Filter edges based on gene
        genes = pd.read_table(genes, header=None)
        Promoter_Gene_filt_g = Promoter_Gene[Promoter_Gene['Target'].isin(genes.iloc[:,0])]
        if promoter_promoter =="true":
            Promoter_Promoter_filt_g = Promoter_Promoter[Promoter_Promoter['Source'].isin(Promoter_Gene_filt_g['Source']) | Promoter_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
        Distal_Promoter_filt_g = Distal_Promoter[Distal_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
        Factor_Promoter_filt_g = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
        Factor_Distal_filt_g = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g['Source'])]

    ### Creating edge table for cytoscape
    if network_mode == "a":
        if promoter_promoter =="true":
            Edges = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Promoter, Promoter_Gene]).drop_duplicates()
        else:
            Edges = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Gene]).drop_duplicates()

    elif network_mode == "f":
        if promoter_promoter =="true":
            if network_distal_only=="true":
                Edges =  Factor_Distal.append([Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
            else:
                Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
        else:
            if network_distal_only=="true":
                Edges =  Factor_Distal.append([Distal_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
            else:
                Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()

    elif network_mode == "g":
        if promoter_promoter =="true":
          Edges =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()
        else:
          Edges =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()

    elif network_mode == "fg":
        if promoter_promoter =="true":
            if network_distal_only=="true":
                Factor_Promoter_filt_fg_2 = Factor_Promoter_filt_fg[Factor_Promoter_filt_fg['Target'].isin(Promoter_Promoter_filt_fg["Target"]) | Factor_Promoter_filt_fg['Target'].isin(Promoter_Promoter_filt_fg["Source"])]
                Edges =  Factor_Distal_filt_fg.append([Factor_Promoter_filt_fg_2, Distal_Promoter_filt_fg, Promoter_Promoter_filt_fg, Promoter_Gene_filt_fg]).drop_duplicates()
            else:
                Edges =  Factor_Distal_filt_fg.append([Factor_Promoter_filt_fg, Distal_Promoter_filt_fg, Promoter_Promoter_filt_fg, Promoter_Gene_filt_fg]).drop_duplicates()
        else:
            if network_distal_only=="true":
                Edges =  Factor_Distal_filt_fg.append([Distal_Promoter_filt_fg, Promoter_Gene_filt_fg]).drop_duplicates()
            else:
                Edges =  Factor_Distal_filt_fg.append([Factor_Promoter_filt_fg, Distal_Promoter_filt_fg, Promoter_Gene_filt_fg]).drop_duplicates()

    Edges.to_csv('Network_Edges_' + prefix + '_interactions.txt', index=False, sep='\t' )


    ### Creating node table for cytoscape
    if network_mode == "a":
    #Specifying node type for all nodes
        Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes.columns=['Node']
        if promoter_promoter =="true":
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Promoter_Gene['Source']) | Nodes['Node'].isin(Distal_Promoter['Target']) | Nodes['Node'].isin(Promoter_Promoter['Source']) | Nodes['Node'].isin(Promoter_Promoter['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))
        else:
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                     (np.where(Nodes['Node'].isin(Promoter_Gene['Source']) | Nodes['Node'].isin(Distal_Promoter['Target']) | Nodes['Node'].isin(Factor_Promoter['Target']), 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))

    elif network_mode == "f":
    # Specifying node type for all nodes that are associated with factor binding
        Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes.columns=['Node']
        if promoter_promoter =="true":
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Factor_Promoter['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))
        else:
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                    (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Factor_Promoter['Target']), 'Promoter',
                                      (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))

    elif network_mode == "g":
        # Specifying node type for all nodes that are associated with selected genes
        Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes.columns=['Node']
        if promoter_promoter =="true":
          Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                        (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                           (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_g['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_g['Target']), 'Promoter',
                                              (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))
        else:
          Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                      (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                         (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Source']), 'Promoter',
                                            (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))
    elif network_mode == "fg":
    # Specifying node type for all nodes that are associated with selected genes
        Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes.columns=['Node']
        if promoter_promoter == "true":
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_fg['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_fg['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter_filt_fg['Source']), 'Distal',
                                        (np.where((Nodes['Node'].isin(Promoter_Promoter_filt_fg['Target'])  | Nodes['Node'].isin(Promoter_Promoter_filt_fg['Source']) | Nodes['Node'].isin(Promoter_Gene_filt_fg['Source'])), 'Promoter',
                                            (np.where(Nodes['Node'].isin(Promoter_Gene_filt_fg['Target']), 'Gene', np.nan)))))))
        else:
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_fg['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_fg['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter_filt_fg['Source']), 'Distal',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene_filt_fg['Source']), 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene_filt_fg['Target']), 'Gene', np.nan)))))))

    Nodes.to_csv('Network_Nodes_' + prefix + '_interactions.txt', index=False, sep='\t' )

    ### Save files for UpSet PLOT
    if (upset_plot == 'true' or complete == 'true' or circos_plot == 'true'):
        Factor_Promoter.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv('UpSet_' + prefix + '_interactions_Promoter.txt', index=False, sep='\t' )
        Factor_Distal.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv('UpSet_' + prefix + '_interactions_Distal.txt', index=False, sep='\t' )

    if filter_genes == 'true':
        genes = pd.read_table(genes, header=None)
        if (circos_use_promoters == "true"):
            Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target'])]
            Promoter_Gene_filt_f_org = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Target'])]
            Promoter_Gene_filt_g_org = Promoter_Gene_filt_f_org[(Promoter_Gene_filt_f_org['Target'].isin(genes.iloc[:,0]))]
            Distal_Promoter_filt_g_org = Distal_Promoter_filt_f[Distal_Promoter_filt_f['Target'].isin(Promoter_Gene_filt_g_org['Source'])]
            Factor_Distal_filt_g_org = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g_org['Source'])]
            Factor_Promoter_filt_g_org = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g_org['Source'])]
            Factor_Distal_filt_g_org.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv('UpSet_' + prefix + '_interactions_Distal_genes.txt', index=False, sep='\t' )
            Factor_Promoter_filt_g_org.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv('UpSet_' + prefix + '_interactions_Promoter_genes.txt', index=False, sep='\t' )

        else:
            Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target'])]
            Promoter_Gene_filt_f_org = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Target'])]
            Promoter_Gene_filt_g_org = Promoter_Gene_filt_f_org[(Promoter_Gene_filt_f_org['Target'].isin(genes.iloc[:,0]))]
            Distal_Promoter_filt_g_org = Distal_Promoter_filt_f[Distal_Promoter_filt_f['Target'].isin(Promoter_Gene_filt_g_org['Source'])]
            Factor_Distal_filt_g_org = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g_org['Source'])]
            Factor_Promoter_filt_g_org = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g_org['Source'])]
            Factor_Distal_filt_g_org.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv('UpSet_' + prefix + '_interactions_Distal_genes.txt', index=False, sep='\t' )
            Factor_Promoter_filt_g_org.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv('UpSet_' + prefix + '_interactions_Promoter_genes.txt', index=False, sep='\t' )



# RUN FUNCTION
network_preprocessing_multiple(interactions_annotated=args.INTERACTIONS_ANNO_AGG, interactions_annotated_not_aggregated=args.INTERACTIONS_ANNO, genes=args.GENES, prefix=args.PREFIX, network_mode=args.NETWORK_MODE, promoter_promoter=args.PROMOTER_PROMOTER, upset_plot=args.UPSET_PLOT, circos_plot=args.CIRCOS_PLOT, filter_genes=args.FILTER_GENES, complete=args.COMPLETE, network_distal_only=args.NETWORK_DISTAL_ONLY, circos_use_promoters=args.CIRCOS_USE_PROMOTERS)
