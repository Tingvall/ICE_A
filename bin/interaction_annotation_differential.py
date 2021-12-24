#!/usr/bin/env python

### PROCESS 9 DIFFERENTIAL: INTERACTION CENTERED ANNOTATION - WITH PEAK OVERLAP ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'INTERACTION-CENTERED ANNOTATION WITH PEAK OVERLAP IN DIFFERENTIAL MODE'
Epilog = """Example usage: interaction_annotation_differential.py <ANCHOR1_PEAK> <ANCHOR2_PEAK> <BED2D> --prefix <PREFIX> --sample <SAMPLE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# General arguments
argParser.add_argument('ANCHOR1_PEAK', help="Annotated anchor1 regions with peak overlap.")
argParser.add_argument('ANCHOR2_PEAK', help="Annotated anchor2 regions with peak overlap.")
argParser.add_argument('BED2D', help="2D-bed interactions.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--sample', dest="SAMPLE", help="Name of sample.")
argParser.add_argument('--network', dest="NETWORK", help="Specifies if files for network visualization in Cytoscape should be created." , choices=['true', 'false'])
argParser.add_argument('--complete', dest="COMPLETE", help="If set to true, all available processes for the selected mode and provided inputs are run.", choices=['true', 'false'])
argParser.add_argument('--binsize', dest="BINSIZE", help="Specifies interaction binsize (default: 5000)", type=int)
argParser.add_argument('--promoter_start', dest="PROMOTER_START", help="Specifies the upstream of TSS considered as a promoter (default: 2500).", type=int)
argParser.add_argument('--promoter_end', dest="PROMOTER_END", help="Specifies the downstream of TSS considered as a promoter (default: 2500).", type=int)

# Differntial mode specific arguments
argParser.add_argument('--peak_differential', dest='PEAK_DIFFERENTIAL', help="Path to textfile that contain log2FC and adjusted p-value from differential analysis. The 1st column should contain peakID matching the peakID in the 4th column of the input bed file. Standard DESeq2 output is expected (with log2FC in the 3rd column and padj in the 9th column), but other formats are accepted as well is the column corresponding to log2FC and padj are specified with the arguments --log2FC_column and --padj column.")
argParser.add_argument('--log2FC_column', dest="LOG2FC_COLUMN", help="Log2FC column for differential peaks.", type=int)
argParser.add_argument('--padj_column', dest="PADJ_COLUMN", help="Padj column for differential peaks.", type=int)
argParser.add_argument('--log2FC', dest="LOG2FC", help="Log2FC threshold for differential peaks.", type=float)
argParser.add_argument('--padj', dest="PADJ", help="Padj threshold for differential peaks.", type=float)

args = argParser.parse_args()

# DEFINE FUNCTION
def interaction_annotation_differential(anchor_1_peak_collect, anchor_2_peak_collect, bed2D_index_anno, prefix, sample, network, complete, binsize, promoter_start, promoter_end, peak_differential, log2FC_column, padj_column, log2FC, padj):

    # Column names for loaded data
    anchor1_peak_name = ('Anchor1_Chr', 'Anchor1_Start', 'Anchor1_End', 'Peak1_Chr', 'Peak1_Start', 'Peak1_End', 'Peak1_name', 'Peak1_score')
    anchor2_peak_name = ('Anchor2_Chr', 'Anchor2_Start', 'Anchor2_End', 'Peak2_Chr', 'Peak2_Start', 'Peak2_End', 'Peak2_name', 'Peak2_score')

    # Load interaction centered peak overlaps 2and annotated 2D-bed
    anchor1_peaks = pd.read_table(anchor_1_peak_collect, index_col=3,names=anchor1_peak_name).sort_index()
    anchor2_peaks = pd.read_table(anchor_2_peak_collect, index_col=3,names=anchor2_peak_name).sort_index()
    bed2D_anno = pd.read_table(bed2D_index_anno, index_col=1).sort_index().iloc[:,1:]

    # Load differntial peaks and merge
    peak_differntial = pd.read_table(peak_differential, index_col=0).sort_index()
    peak_differntial = peak_differntial.iloc[:, [log2FC_column-2, padj_column-2]]
    peak_differntial.columns = ['log2FC', 'padj']
    anchor1_peaks=anchor1_peaks.merge(peak_differntial, left_on="Peak1_name", right_index=True, how = 'left')
    anchor2_peaks=anchor2_peaks.merge(peak_differntial, left_on="Peak2_name", right_index=True, how = 'left')

    # Create Peak columns (chr:start-end) for anchor 1 & 2
    anchor1_peaks["Peak1_ID"] = anchor1_peaks["Peak1_Chr"].map(str) +':'+ (anchor1_peaks["Peak1_Start"]-1).map(str) +'-'+ anchor1_peaks["Peak1_End"].map(str)
    anchor2_peaks["Peak2_ID"] = anchor2_peaks["Peak2_Chr"].map(str) +':'+ (anchor2_peaks["Peak2_Start"]-1).map(str) +'-'+ anchor2_peaks["Peak2_End"].map(str)

    # Merging anchor points
    anchor1_peaks_anno =bed2D_anno.loc[:,['chr1', 's1','e1', 'Entrez_ID_1', 'Gene_Name_1', 'Distance_to_TSS_1', 'TSS_1']].merge(anchor1_peaks.loc[:,['Peak1_ID', 'Peak1_score']], left_index=True, right_index=True, how = 'left')
    anchor2_peaks_anno =bed2D_anno.loc[:,['chr2', 's2','e2' ,'Entrez_ID_2', 'Gene_Name_2','Distance_to_TSS_2', "TSS_2"]].merge(anchor2_peaks.loc[:,['Peak2_ID', 'Peak2_score']], left_index=True, right_index=True, how = 'left').merge(bed2D_anno.loc[:,['Interaction_score']], left_index=True, right_index=True, how = 'left')
    anchors_peaks_anno = anchor1_peaks_anno.merge(anchor2_peaks_anno, left_index=True, right_index=True,how = 'outer').drop_duplicates()
    anchors_peaks_anno.rename(columns = {'TSS_1': 'Is_Promoter_1', 'TSS_2': 'Is_Promoter_2'}, inplace = True)
    anchors_peaks_anno['Is_Promoter_1'] = np.where(((anchors_peaks_anno['Distance_to_TSS_1'] <= 0) & (anchors_peaks_anno['Distance_to_TSS_1'] >= -(binsize/2+promoter_end))) | ((anchors_peaks_anno['Distance_to_TSS_1'] >= 0) & (anchors_peaks_anno['Distance_to_TSS_1'] <= (binsize/2+promoter_start))),1,0)
    anchors_peaks_anno['Is_Promoter_2'] = np.where(((anchors_peaks_anno['Distance_to_TSS_2'] <= 0) & (anchors_peaks_anno['Distance_to_TSS_2'] >= -(binsize/2+promoter_end))) | ((anchors_peaks_anno['Distance_to_TSS_2'] >= 0) & (anchors_peaks_anno['Distance_to_TSS_2'] <= (binsize/2+promoter_start))),1,0)

    # Creation and use of function for adding 2 columns for peak overlap in anchor points (overlap in anchor 1/2) with 1 if overlap
    def peak_in_anchor_1(row):
        if pd.isna(row['Peak1_ID']) :
            return ''
        else:
            return 1
    def peak_in_anchor_2(row):
        if pd.isna(row['Peak2_ID']) :
            return ''
        else:
            return 1

    anchors_peaks_anno.index.name = 'Interaction'
    anchors_peaks_anno['Overlap_1'] = anchors_peaks_anno.apply (lambda row: peak_in_anchor_1(row), axis=1)
    anchors_peaks_anno['Overlap_2'] = anchors_peaks_anno.apply (lambda row: peak_in_anchor_2(row), axis=1)
    anchors_peaks_anno_original = anchors_peaks_anno.copy(deep=True)
    anchors_peaks_anno_factor = anchors_peaks_anno[(anchors_peaks_anno['Overlap_1'] == 1) | (anchors_peaks_anno['Overlap_2'] == 1)]

    # Extracting interactions associates with differntial peaks
    anchors_peaks_anno_up=anchors_peaks_anno_factor[((anchors_peaks_anno_factor.padj_1 <= padj) & (anchors_peaks_anno_factor.log2FC_1 >= log2FC)) | ((anchors_peaks_anno_factor.padj_2 <= padj) & (anchors_peaks_anno_factor.log2FC_2 >= log2FC)) ]
    anchors_peaks_anno_down=anchors_peaks_anno_factor[((anchors_peaks_anno_factor.padj_1 <= padj) & (anchors_peaks_anno_factor.log2FC_1 <= log2FC)) | ((anchors_peaks_anno_factor.padj_2 <= padj) & (anchors_peaks_anno_factor.log2FC_2 <= -log2FC)) ]

    #Merging for one row per interaction
    anchors_peaks_anno_factor = anchors_peaks_anno_factor.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))
    anchors_peaks_anno_up = anchors_peaks_anno_up.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))
    anchors_peaks_anno_down = anchors_peaks_anno_down.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))
    anchors_peaks_anno = anchors_peaks_anno.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    # Saving annotated interactions files (all interaction and interactions with peak overlap)
    anchors_peaks_anno_factor.to_csv(sample + '_' + prefix + '_interactions.txt', index=False, sep='\t' )
    anchors_peaks_anno_up.to_csv(sample + '_' + prefix + '_interactions_up.txt', index=False, sep='\t' )
    anchors_peaks_anno_down.to_csv(sample + '_' + prefix + '_interactions_down.txt', index=False, sep='\t' )
    anchors_peaks_anno.to_csv(prefix + '_HOMER_annotated_interactions_with_peak_overlap.txt', index=True, sep='\t' )


    # Save files for Network
    if (network == 'true' or complete == 'true'):
      anchors_peaks_anno_original.to_csv(prefix + '_HOMER_annotated_interactions_with_peak_overlap_not_aggregated.txt', index=True, sep='\t' )


# RUN FUNCTION
interaction_annotation_differential(anchor_1_peak_collect=args.ANCHOR1_PEAK,anchor_2_peak_collect=args.ANCHOR2_PEAK,bed2D_index_anno=args.BED2D, prefix=args.PREFIX, sample=args.SAMPLE, network=args.NETWORK, complete=args.COMPLETE, binsize=args.BINSIZE, promoter_start=args.PROMOTER_START, promoter_end=args.PROMOTER_END, peak_differential=args.PEAK_DIFFERENTIAL, log2FC_column=args.LOG2FC_COLUMN, padj_column=args.PADJ_COLUMN, log2FC=args.LOG2FC, padj=args.PADJ)
