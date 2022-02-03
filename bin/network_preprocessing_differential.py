#!/usr/bin/env python

### PROCESS 10 DIFFERENTIAL: CREATES NETWORK FILES FOR CYTOSCAPE ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'PREPROCESSING OF ANNOTATED INTERATION FOR CYTOSCAPE NETWORK VISUALIZATION IN DIFFERENTIAL MODE'
Epilog = """Example usage: network_preprocessing_differential.py <INTERACTIONS_ANNO_AGG> <INTERACTIONS_ANNO> --genes <GENES> --peak_differential <PEAK_DIFFERENTIAL> --expression <EXPRESSION> --sample <SAMPLE> --prefix <PREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# General arguments
argParser.add_argument('INTERACTIONS_ANNO_AGG', help="Annotated and aggregated interactions.")
argParser.add_argument('INTERACTIONS_ANNO', help="Annotated, not aggregated interactions")
argParser.add_argument('--genes', dest='GENES', help="Text file specifying genes for filtering.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--sample', dest="SAMPLE", help="Sample name.")
argParser.add_argument('--network_mode', dest="NETWORK_MODE", help="Defines mode network. Options are all (all interaction in the 2D-bed file), factor (all interaction with at least on peak overlap either anchor point) or genes (interactions associates with a gene list, provided by --genes)." , choices=['all', 'factor', 'genes', 'expression', 'differential', 'factorgenes', 'factorexpression'])
argParser.add_argument('--promoter_promoter', dest="PROMOTER_PROMOTER", help="If set to true, promoter-promoter interactions included in network (default: false).", choices=['true', 'false'])
argParser.add_argument('--complete', dest="COMPLETE", help="If set to true, all available processes for the selected mode and provided inputs are run.", choices=['true', 'false'])
argParser.add_argument('--network_distal_only', dest="NETWORK_DISTAL_ONLY", help="If true, only distal factor binding are shown in the netork.", choices=['true', 'false'])

# Differntial mode specific arguments
argParser.add_argument('--peak_differential', dest='PEAK_DIFFERENTIAL', help="Path to textfile that contain log2FC and adjusted p-value from differential analysis. The 1st column should contain peakID matching the peakID in the 4th column of the input bed file. Standard DESeq2 output is expected (with log2FC in the 3rd column and padj in the 9th column), but other formats are accepted as well is the column corresponding to log2FC and padj are specified with the arguments --log2FC_column and --padj column.")
argParser.add_argument('--expression', dest='EXPRESSION', help="Specifies path to file that contain information about differential expression between the two conditions. The first column must contain gene symbol. The column for log2FC/padj can be specified using --expression_log2FC_column and --expression_padj_column respectively, default: 3 & 9 (standard DESeq2 format). Note: make sure that you use the same direction for the comparison in --peak_differential and --expression.")
argParser.add_argument('--log2FC_column', dest="LOG2FC_COLUMN", help="Log2FC column for differential peaks.", type=int)
argParser.add_argument('--padj_column', dest="PADJ_COLUMN", help="Padj column for differential peaks.", type=int)
argParser.add_argument('--log2FC', dest="LOG2FC", help="Log2FC threshold for differential peaks.", type=float)
argParser.add_argument('--padj', dest="PADJ", help="Padj threshold for differential peaks.", type=float)
argParser.add_argument('--skip_expression', dest="SKIP_EXPRESSION", help="Use this argument if no --expression file is provided.", choices=['true', 'false'])
argParser.add_argument('--expression_log2FC_column', dest="EXPRESSION_LOG2FC_COLUMN", help="Log2FC column for differential expression.", type=int)
argParser.add_argument('--expression_padj_column', dest="EXPRESSION_PADJ_COLUMN", help="Padj column for differential expression.", type=int)

args = argParser.parse_args()

# DEFINE FUNCTION
def network_preprocessing_differential(interactions_annotated, interactions_annotated_not_aggregated, genes, prefix, sample, network_mode, promoter_promoter, peak_differential, expression, log2FC_column, padj_column, log2FC, padj, skip_expression, expression_log2FC_column, expression_padj_column, complete, network_distal_only):

    #Loading input file
    anchors_peaks_anno = pd.read_table(interactions_annotated_not_aggregated, index_col=0)
    interactions_anno = pd.read_table(interactions_annotated, index_col=0)

    if skip_expression == "false":
    expression = pd.read_table(expression, index_col=0).sort_index()
    expression = expression.iloc[:, [expression_log2FC_column-2, expression_padj_column-2]]
    expression.columns = ['log2FC', 'padj']
    expression_diff = expression.loc[(expression['padj'] <= padj) & (abs(expression['log2FC']) >= log2FC),:]

    # Aggregating interaction file to only incude one row per interaction
    interactions_anno = interactions_anno.iloc[:,np.r_[0:5,6,11:16,17, 22:len(anchors_peaks_anno.columns)]]
    interactions_anno['Anchor1'] = interactions_anno["chr1"].map(str) +':'+ (interactions_anno["s1"]).map(str) +'-'+ interactions_anno["e1"].map(str)
    interactions_anno['Anchor2'] = interactions_anno["chr2"].map(str) +':'+ (interactions_anno["s2"]).map(str) +'-'+ interactions_anno["e2"].map(str)
    interactions_anno = pd.concat([interactions_anno['Anchor1'], interactions_anno.iloc[:,3:6], interactions_anno['Anchor2'],interactions_anno.iloc[:,9:(len(interactions_anno.columns)-2)]], axis=1)

    ### Creating edge table for cytoscape
    #Factor-Interaction
    Factor_Interaction = anchors_peaks_anno.copy(deep=True)
    Factor_Interaction.loc[Factor_Interaction.Overlap_1 == 1, 'Peak1'] = sample
    Factor_Interaction.loc[Factor_Interaction.Overlap_2 == 1, 'Peak2'] = sample
    Factor_Interaction = Factor_Interaction[['chr1', 's1', 'e1','Gene_Name_1', 'Peak1','Peak1_ID','Peak1_score', 'log2FC_1', 'padj_1','chr2', 's2', 'e2',  'Gene_Name_2','Peak2','Peak2_ID','Peak2_score','log2FC_2', 'padj_2', 'Is_Promoter_1', 'Is_Promoter_2']]
    Factor_Interaction['Anchor1'] = Factor_Interaction['chr1'].map(str) +':'+ (Factor_Interaction['s1']).map(str) +'-'+ Factor_Interaction['e1'].map(str)
    Factor_Interaction['Anchor2'] = Factor_Interaction['chr2'].map(str) +':'+ (Factor_Interaction['s2']).map(str) +'-'+ Factor_Interaction['e2'].map(str)
    Factor_Interaction = Factor_Interaction.dropna(subset=['Peak1', 'Peak2'], thresh=1)

    #Factor-Distal
    Factor_Distal_1 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 0) & (Factor_Interaction['Is_Promoter_2'] == 1), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Distal_1.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_2 = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['Is_Promoter_2'] == 0), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Distal_2.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal = Factor_Distal_1.append(Factor_Distal_2)
    Factor_Distal['Edge_type'] = 'Factor-Distal'

    Factor_Distal_1_diff = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 0) & (Factor_Interaction['Is_Promoter_2'] == 1) & (Factor_Interaction['padj_1'] <= padj) & (abs(Factor_Interaction['log2FC_1']) >= log2FC), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Distal_1_diff.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_2_diff = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['Is_Promoter_2'] == 0) & (Factor_Interaction['padj_2'] <= padj) & (abs(Factor_Interaction['log2FC_2']) >= log2FC), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Distal_2_diff.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_diff = Factor_Distal_1_diff.append(Factor_Distal_2_diff)
    Factor_Distal_diff['Edge_type'] = 'Factor-Distal'

    Factor_Distal_1_up = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 0) & (Factor_Interaction['Is_Promoter_2'] == 1) & (Factor_Interaction['padj_1'] <= padj) & (Factor_Interaction['log2FC_1'] >= log2FC), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Distal_1_up.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_2_up = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['Is_Promoter_2'] == 0) & (Factor_Interaction['padj_2'] <= padj) & (Factor_Interaction['log2FC_2'] >= log2FC), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Distal_2_up.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_up = Factor_Distal_1_up.append(Factor_Distal_2_up)
    Factor_Distal_up['Edge_type'] = 'Factor-Distal'

    Factor_Distal_1_down = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 0) & (Factor_Interaction['Is_Promoter_2'] == 1) & (Factor_Interaction['padj_1'] <= padj) & (Factor_Interaction['log2FC_1'] <= -log2FC), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Distal_1_down.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_2_down = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['Is_Promoter_2'] == 0) & (Factor_Interaction['padj_2'] <= padj) & (Factor_Interaction['log2FC_2'] <= -log2FC), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Distal_2_down.columns = ['Source', 'Target', 'Edge_score']
    Factor_Distal_down = Factor_Distal_1_down.append(Factor_Distal_2_down)
    Factor_Distal_down['Edge_type'] = 'Factor-Distal'

    #Factor-Promoter
    Factor_Promoter_1 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_1'] == 1, ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Promoter_1.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_2 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_2'] == 1, ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Promoter_2.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter = Factor_Promoter_1.append(Factor_Promoter_2)
    Factor_Promoter['Edge_type'] = 'Factor-Promoter'

    Factor_Promoter_1_diff = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['padj_1'] <= padj) & ((abs(Factor_Interaction['log2FC_1']) >= log2FC)), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Promoter_1_diff.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_2_diff = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['padj_1'] <= padj) & ((abs(Factor_Interaction['log2FC_1']) >= log2FC)), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Promoter_2_diff.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_diff = Factor_Promoter_1_diff.append(Factor_Promoter_2_diff)
    Factor_Promoter_diff['Edge_type'] = 'Factor-Promoter'

    Factor_Promoter_1_up = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['padj_1'] <= padj) & (Factor_Interaction['log2FC_1'] >= log2FC), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Promoter_1_up.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_2_up = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['padj_2'] <= padj) & (Factor_Interaction['log2FC_2'] >= log2FC), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Promoter_2_up.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_up = Factor_Promoter_1_up.append(Factor_Promoter_2_up)
    Factor_Promoter_up['Edge_type'] = 'Factor-Promoter'

    Factor_Promoter_1_down = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['padj_1'] <= padj) & (Factor_Interaction['log2FC_1'] <= -log2FC), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
    Factor_Promoter_1_down.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_2_down = Factor_Interaction.loc[(Factor_Interaction['Is_Promoter_1'] == 1) & (Factor_Interaction['padj_2'] <= padj) & (Factor_Interaction['log2FC_2'] <= -log2FC), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
    Factor_Promoter_2_down.columns = ['Source', 'Target', 'Edge_score']
    Factor_Promoter_down = Factor_Promoter_1_down.append(Factor_Promoter_2_down)
    Factor_Promoter_down['Edge_type'] = 'Factor-Promoter'

    #Distal-Promoter
    DP_1 = interactions_anno.loc[(interactions_anno['Is_Promoter_1'] == 0) & (interactions_anno['Is_Promoter_2'] == 1), ['Anchor1','Anchor2', 'Interaction_score']]
    DP_1.columns = ['Source', 'Target', 'Edge_score']
    DP_2 = interactions_anno.loc[(interactions_anno['Is_Promoter_1'] == 1) & (interactions_anno['Is_Promoter_2'] == 0), ['Anchor2',  'Anchor1', 'Interaction_score']]
    DP_2.columns = ['Source', 'Target', 'Edge_score']
    Distal_Promoter = DP_1.append(DP_2)
    Distal_Promoter['Edge_type'] = 'Distal-Promoter'
    Distal_Promoter['Edge_score'] = - np.log10(Distal_Promoter['Edge_score'])

    #Promoter-Promoter
    Promoter_Promoter = interactions_anno.loc[(interactions_anno['Is_Promoter_1']==1) & (interactions_anno['Is_Promoter_2']==1),:][['Anchor1', 'Anchor2', 'Interaction_score']]
    Promoter_Promoter['Edge_type'] = 'Promoter-Promoter'
    Promoter_Promoter.columns = ['Source', 'Target', 'Edge_score', 'Edge_type']
    Promoter_Promoter['Edge_score'] = - np.log10(Promoter_Promoter['Edge_score'])

    #Promoter-Gene
    Promoter_Gene_1 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_1'] == 1, ['Anchor1', 'Gene_Name_1']].dropna(subset=['Gene_Name_1']).drop_duplicates()
    Promoter_Gene_1.columns = ['Source', 'Target']
    Promoter_Gene_2 = Factor_Interaction.loc[Factor_Interaction['Is_Promoter_2'] == 1, ['Anchor2',  'Gene_Name_2']].dropna(subset=['Gene_Name_2']).drop_duplicates()
    Promoter_Gene_2.columns = ['Source', 'Target']
    Promoter_Gene = Promoter_Gene_1.append(Promoter_Gene_2)
    Promoter_Gene['Edge_score'], Promoter_Gene['Edge_type'] = [1, 'Promoter-Gene']

    # Filtering of edges based on network mode
    if (network_mode == "factor" or network_mode == "factorgenes" or network_mode == "factorexpression"):
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

    elif (network_mode == "genes" or network_mode == "expression"):
        #Filter edges based on gene
        if network_mode == "genes":
            genes = pd.read_table(genes, header=None)

        elif network_mode == "expression":
            genes = pd.DataFrame(pd.unique(expression_diff.index.dropna().values.ravel('K')))

        Promoter_Gene_filt_g = Promoter_Gene[Promoter_Gene['Target'].isin(genes.iloc[:,0])]
        if promoter_promoter =="true":
            Promoter_Promoter_filt_g = Promoter_Promoter[Promoter_Promoter['Source'].isin(Promoter_Gene_filt_g['Source']) | Promoter_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
        Distal_Promoter_filt_g = Distal_Promoter[Distal_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
        Factor_Promoter_filt_g = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
        Factor_Distal_filt_g = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g['Source'])]


        if (network_mode == "factorgenes" or network_mode == "factorexpression"):
            #Filter edges based on gene
            if network_mode == "factorgenes":
                genes = pd.read_table(genes, header=None)

            elif network_mode == "factorexpression":
                genes = pd.DataFrame(pd.unique(expression_diff.index.dropna().values.ravel('K')))

            if promoter_promoter =="true":
                if network_distal_only=="true":
                    Promoter_Gene_filt_fg = Promoter_Gene_filt_f[(Promoter_Gene_filt_f['Target'].isin(genes.iloc[:,0])) & (Promoter_Gene_filt_f['Source'].isin(Distal_Promoter_filt_f['Target']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Source']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Target']))]
                else:
                    Promoter_Gene_filt_gf = Promoter_Gene_filt_f[(Promoter_Gene_filt_f['Target'].isin(genes.iloc[:,0])) & (Promoter_Gene_filt_f['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene_filt_f['Source'].isin(Distal_Promoter_filt_f['Target']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Source']) | Promoter_Gene_filt_f['Source'].isin(Promoter_Promoter_filt_f['Target']))]
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

    elif network_mode == "differential":
        #Filter edges based on differnetial peaks
        Distal_Promoter_filt_diff = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal_diff['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter_diff['Target'])]
        Distal_Promoter_filt_up = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal_up['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter_up['Target'])]
        Distal_Promoter_filt_down = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal_down['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter_down['Target'])]

        if promoter_promoter =="true":
            Promoter_Promoter_filt_diff = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter_diff['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter_diff['Target'])]
            if network_distal_only=="true":
                Promoter_Promoter_filt_diff = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_diff['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_diff['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_diff['Target'])]
                Promoter_Promoter_filt_up = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_up['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_up['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_up['Target'])]
                Promoter_Promoter_filt_down = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_down['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_down['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_down['Target'])]
            else:
                Promoter_Promoter_filt_diff = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_diff['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_diff['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_diff['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_diff['Target'])]
                Promoter_Promoter_filt_up = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_up['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_up['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_up['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_up['Target'])]
                Promoter_Promoter_filt_down = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_down['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_down['Target']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_down['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_down['Target'])]

        else:
            if network_distal_only=="true":
                Promoter_Promoter_filt_diff = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_diff['Target'])]
                Promoter_Promoter_filt_up = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_up['Target'])]
                Promoter_Promoter_filt_down = Promoter_Gene[Promoter_Gene['Source'].isin(Distal_Promoter_filt_down['Target'])]
            else:
                Promoter_Promoter_filt_diff = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_diff['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_diff['Target'])]
                Promoter_Promoter_filt_up = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_up['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_up['Target'])]
                Promoter_Promoter_filt_down = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_down['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_down['Target'])]

    elif (network_mode == "genes" or network_mode == "expression"):
          #Filter edges based on gene
          if network_mode == "genes":
              genes = pd.read_table(genes, header=None)

          elif network_mode == "expression":
              genes = pd.DataFrame(pd.unique(expression_diff.index.dropna().values.ravel('K')))

          Promoter_Gene_filt_g = Promoter_Gene[Promoter_Gene['Target'].isin(genes.iloc[:,0])]
          if promoter_promoter =="true":
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
            if network_distal_only=="true":
                Edges =  Factor_Distal.append([Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
            else:
                Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
        else:
            if network_distal_only=="true":
                Edges =  Factor_Distal.append([Distal_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
            else:
                Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()

    elif (network_mode == "factorgenes" or network_mode == "factorexpression"):
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

    elif network_mode == "differential":
        if promoter_promoter =="true":
            if network_distal_only=="true":
                Edges =  Factor_Distal_diff.append([Distal_Promoter_diff, Promoter_Promoter_filt_diff, Promoter_Gene_filt_diff]).drop_duplicates()
                Edges_up =  Factor_Distal_up.append([Distal_Promoter_up, Promoter_Promoter_filt_up, Promoter_Gene_filt_up]).drop_duplicates()
                Edges_down =  Factor_Distal_diff.append([Distal_Promoter_down, Promoter_Promoter_filt_down, Promoter_Gene_filt_down]).drop_duplicates()

            else:
                Edges =  Factor_Distal_diff.append([Factor_Promoter_diff, Distal_Promoter_filt_diff, Promoter_Promoter_filt_diff, Promoter_Gene_filt_diff]).drop_duplicates()
                Edges_up =  Factor_Distal_up.append([Factor_Promoter_up, Distal_Promoter_filt_up, Promoter_Promoter_filt_up, Promoter_Gene_filt_up]).drop_duplicates()
                Edges_down =  Factor_Distal_down.append([Factor_Promoter_down, Distal_Promoter_filt_down, Promoter_Promoter_filt_down, Promoter_Gene_filt_down]).drop_duplicates()
        else:
            if network_distal_only=="true":
                Edges =  Factor_Distal.append([Factor_Promoter_diff, Promoter_Gene_filt_diff]).drop_duplicates()
                Edges_up =  Factor_Distal.append([Factor_Promoter_up, Promoter_Gene_filt_up]).drop_duplicates()
                Edges_down =  Factor_Distal.append([Factor_Promoter_down, Promoter_Gene_filt_down]).drop_duplicates()
            else:
                Edges =  Factor_Distal.append([Factor_Promoter_diff, Distal_Promoter_filt_diff, Promoter_Gene_filt_diff]).drop_duplicates()
                Edges_up =  Factor_Distal.append([Factor_Promoter_up, Distal_Promoter_filt_up, Promoter_Gene_filt_up]).drop_duplicates()
                Edges_down =  Factor_Distal.append([Factor_Promoter_down, Distal_Promoter_filt_down, Promoter_Gene_filt_down]).drop_duplicates()

    elif (network_mode == "genes" or network_mode == "expression"):
        if promoter_promoter =="true":
            Egdes =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()
        else:
            Egdes =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()

        Egdes_up.to_csv('Network_Edges_' + prefix + '_interactions_up.txt', index=False, sep='\t' )
        Egdes_down.to_csv('Network_Edges_' + prefix + '_interactions_down.txt', index=False, sep='\t' )

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
                                     (np.where(Nodes['Node'].isin(Distal_Promoter['Target']) | nodes['Node'].isin(Factor_Promoter['Target']), 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))

    elif network_mode == "factor":
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

    elif (network_mode == "factorgenes" or network_mode == "factorexpression" ):
    # Specifying node type for all nodes that are associated with selected genes
        Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes.columns=['Node']
        if promoter_promoter == "true":
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                        (np.where((Nodes['Node'].isin(Promoter_Promoter_filt_g['Target'])  | Nodes['Node'].isin(Promoter_Promoter_filt_g['Source']) | Nodes['Node'].isin(Promoter_Gene_filt_g['Source'])), 'Promoter',
                                            (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))
        else:
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Source']), 'Promoter',
                                            (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))

    elif network_mode == "differential":
    # Specifying node type for all nodes that are associated with factor binding
        Nodes = pd.DataFrame(pd.unique(Edges[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes.columns=['Node']
        Nodes_up = pd.DataFrame(pd.unique(Edges_up[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes_up.columns=['Node']
        Nodes_down = pd.DataFrame(pd.unique(Edges_down[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes_down.columns=['Node']
        if promoter_promoter =="true":
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_diff['Source']) | Nodes['Node'].isin(Factor_Distal_diff['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter_filt_diff['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Promoter_Gene_filt_diff['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_diff['Target']) | Nodes['Node'].isin(Factor_Promoter_diff['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_diff['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_diff['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene_filt_diff['Target']), 'Gene', np.nan)))))))
            Nodes_up['Node_type'] = np.where(Nodes_up['Node'].isin(Factor_Distal_up['Source']) | Nodes['Node'].isin(Factor_Promoter_up['Source']), 'Factor',
                                  (np.where(Nodes_up['Node'].isin(Distal_Promoter_filt_up['Source']), 'Distal',
                                     (np.where(Nodes_up['Node'].isin(Promoter_Gene_filt_up['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_up['Target']) | Nodes['Node'].isin(Factor_Promoter_up'Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_up['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_up['Target']), 'Promoter',
                                        (np.where(Nodes_up['Node'].isin(Promoter_Gene_filt_up['Target']), 'Gene', np.nan)))))))
            Nodes_down['Node_type'] = np.where(Nodes_down['Node'].isin(Factor_Distal_down['Source']) | Nodes['Node'].isin(Factor_Promoter_down['Source']), 'Factor',
                                    (np.where(Nodes_down['Node'].isin(Distal_Promoter_filt_down['Source']), 'Distal',
                                       (np.where(Nodes_down['Node'].isin(Promoter_Gene_filt_down['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_down['Target']) | Nodes['Node'].isin(Factor_Promoter_down['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_down['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_down['Target']), 'Promoter',
                                          (np.where(Nodes_down['Node'].isin(Promoter_Gene_filt_down['Target']), 'Gene', np.nan)))))))

        else:
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_diff['Source']) | Nodes['Node'].isin(Factor_Distal_diff['Source']), 'Factor',
                                (np.where(Nodes['Node'].isin(Distal_Promoter_filt_diff['Source']), 'Distal',
                                    (np.where(Nodes['Node'].isin(Promoter_Gene_filt_diff['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_diff['Target']) | Nodes['Node'].isin(Factor_Promoter_diff['Target']), 'Promoter',
                                      (np.where(Nodes['Node'].isin(Promoter_Gene_filt_diff['Target']), 'Gene', np.nan)))))))
            Nodes_up['Node_type'] = np.where(Nodes_up['Node'].isin(Factor_Distal_up['Source']) | Nodes['Node'].isin(Factor_Distal_up['Source']), 'Factor',
                              (np.where(Nodes_up['Node'].isin(Distal_Promoter_filt_up['Source']), 'Distal',
                                  (np.where(Nodes_up['Node'].isin(Promoter_Gene_filt_up['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_up['Target']) | Nodes['Node'].isin(Factor_Promoter_up['Target']), 'Promoter',
                                    (np.where(Nodes_up['Node'].isin(Promoter_Gene_filt_up['Target']), 'Gene', np.nan)))))))
            Nodes_down['Node_type'] = np.where(Nodes_down['Node'].isin(Factor_Distal_down['Source']) | Nodes['Node'].isin(Factor_Distal_down['Source']), 'Factor',
                                (np.where(Nodes_down['Node'].isin(Distal_Promoter_filt_down['Source']), 'Distal',
                                    (np.where(Nodes_down['Node'].isin(Promoter_Gene_filt_down['Source']) | Nodes['Node'].isin(Distal_Promoter_filt_down['Target']) | Nodes['Node'].isin(Factor_Promoter_down['Target']), 'Promoter',
                                      (np.where(Nodes_down['Node'].isin(Promoter_Gene_filt_down['Target']), 'Gene', np.nan)))))))

        Factor_stat_up = Factor_stat.loc[(Factor_stat['log2FC'] >= log2FC) & (Factor_stat['padj'] <= padj),:]
        Factor_stat_up = Factor_stat_up.sort_values('padj').drop_duplicates(subset=['Anchor'],keep='first').sort_index()
        Nodes_peak_anno_up = Nodes_up[Nodes_up.Node_type.isin(('Distal', 'Promoter'))].merge(Factor_stat_up.loc[:,['Anchor', 'padj', 'log2FC']], how='left', right_on='Anchor', left_on='Node')
        Nodes_gene_anno_up = Nodes_up[Nodes_up.Node_type =='Gene'].merge(expression.loc[:,['padj', 'log2FC']], how='left', right_index=True, left_on='Node')
        Nodes_anno_up = Nodes_peak_anno_up.append(Nodes_gene_anno_up)
        Nodes_up=Nodes_up.merge(Nodes_anno_up[['Node','padj', 'log2FC']], on='Node', how='left')
        Nodes_up.to_csv('Network_Nodes_' + prefix + '_interactions_up.txt', index=False, sep='\t' )

        Factor_stat_down = Factor_stat.loc[(Factor_stat['log2FC'] <= -log2FC) & (Factor_stat['padj'] <= padj),:]
        Factor_stat_down = Factor_stat_down.sort_values('padj').drop_duplicates(subset=['Anchor'],keep='first').sort_index()
        Nodes_peak_anno_down = Nodes_down[Nodes_down.Node_type.isin(('Distal', 'Promoter'))].merge(Factor_stat_down.loc[:,['Anchor', 'padj', 'log2FC']], how='left', right_on='Anchor', left_on='Node')
        Nodes_gene_anno_down = Nodes_down[Nodes_down.Node_type =='Gene'].merge(expression.loc[:,['padj', 'log2FC']], how='left', right_index=True, left_on='Node')
        Nodes_anno_down = Nodes_peak_anno_down.append(Nodes_gene_anno_down)
        Nodes_down=Nodes_down.merge(Nodes_anno_down[['Node','padj', 'log2FC']], on='Node', how='left')
        Nodes_down.to_csv('Network_Nodes_' + prefix + '_interactions_down.txt', index=False, sep='\t' )

    elif (network_mode == "genes" or network_mode == "expression" ):
        # Specifying node type for all nodes that are associated with selected genes
        Nodes = pd.DataFrame(pd.unique(Egdes[['Source', 'Target']].dropna().values.ravel('K')))
        Nodes.columns=['Node']
        if promoter_promoter =="true":
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                        (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                           (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_g['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_g['Target']), 'Promoter',
                                              (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))
        else:
            Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_filt_g['Source']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                      (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                         (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Target']) | Nodes['Node'].isin(Factor_Promoter_filt_g['Target']), 'Promoter',
                                            (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))

    #Adding node stats
    Factor_stat = Factor_stat.sort_values('padj').drop_duplicates(subset=['Anchor'],keep='first').sort_index()
    Nodes_peak_anno = Nodes[Nodes.Node_type.isin(('Distal', 'Promoter'))].merge(Factor_stat.loc[:,['Anchor', 'padj', 'log2FC']], how='left', right_on='Anchor', left_on='Node')

    if skip_expression == "false":
        Nodes_gene_anno = Nodes[Nodes.Node_type =='Gene'].merge(expression.loc[:,['padj', 'log2FC']], how='left', right_index=True, left_on='Node')
        Nodes_anno = Nodes_peak_anno.append(Nodes_gene_anno)
    else:
        Nodes_anno = Nodes_peak_anno

    Nodes=Nodes.merge(Nodes_anno[['Node','padj', 'log2FC']], on='Node', how='left')
    Nodes.to_csv('Network_Nodes_' + prefix + '_interactions.txt', index=False, sep='\t' )

# RUN FUNCTION
network_preprocessing_differential(interactions_annotated=args.INTERACTIONS_ANNO_AGG, interactions_annotated_not_aggregated=args.INTERACTIONS_ANNO, genes=args.GENES, prefix=args.PREFIX, sample=args.SAMPLE, network_mode=args.NETWORK_MODE, promoter_promoter=args.PROMOTER_PROMOTER, peak_differential=args.PEAK_DIFFERENTIAL, expression=args.EXPRESSION, log2FC_column=args.LOG2FC_COLUMN, padj_column=args.PADJ_COLUMN, log2FC=args.LOG2FC, padj=args.PADJ, skip_expression=args.SKIP_EXPRESSION, expression_log2FC_column=args.EXPRESSION_LOG2FC_COLUMN, expression_padj_column=args.EXPRESSION_PADJ_COLUMN, complete=args.COMPLETE, network_distal_only=args.NETWORK_DISTAL_ONLY)
