#!/usr/bin/env python

### PROCESS 9 MULTIPLE: INTERACTION CENTERED ANNOTATION - WITH PEAK OVERLAP ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'INTERACTION-CENTERED ANNOTATION WITH PEAK OVERLAP IN MULTIPLE MODE'
Epilog = """Example usage: interaction_annotation_multiple.py <ANCHOR1_PEAK> <ANCHOR2_PEAK> <BED2D> --prefix <PREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# General arguments
argParser.add_argument('ANCHOR1_PEAK', help="Annotated anchor1 regions with peak overlap.")
argParser.add_argument('ANCHOR2_PEAK', help="Annotated anchor2 regions with peak overlap.")
argParser.add_argument('BED2D', help="2D-bed interactions.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--network', dest="NETWORK", help="Specifies if files for network visualization in Cytoscape should be created." , choices=['true', 'false'])
argParser.add_argument('--complete', dest="COMPLETE", help="If set to true, all available processes for the selected mode and provided inputs are run.", choices=['true', 'false'])
argParser.add_argument('--binsize', dest="BINSIZE", help="Specifies interaction binsize (default: 5000)", type=int)
argParser.add_argument('--promoter_start', dest="PROMOTER_START", help="Specifies the upstream of TSS considered as a promoter (default: 2500).", type=int)
argParser.add_argument('--promoter_end', dest="PROMOTER_END", help="Specifies the downstream of TSS considered as a promoter (default: 2500).", type=int)
argParser.add_argument('--in_regions', dest="IN_REGIONS", help="Specify regions to filter input beds. Defult:Not_specified")

# Multiple mode specific arguments
argParser.add_argument('--upset_plot', dest="UPSET_PLOT", help="Specifies if Upset plot of peak overlap will be created.", choices=['true', 'false'])
argParser.add_argument('--circos_plot', dest="CIRCOS_PLOT", help="Specifies if Circos plot of peak overlap will be created.", choices=['true', 'false'])
argParser.add_argument('--circos_use_promoters', dest="CIRCOS_USE_PROMOTERS", help="Specifies if TF overlap in promoters (defined based on promoter_start/end) should be used in circos plot in multiple mode when regions are specified.", choices=['true', 'false'])

args = argParser.parse_args()

# DEFINE FUNCTION
def interaction_annotation_multiple(anchor_1_peak_collect, anchor_2_peak_collect, bed2D_index_anno, prefix, network, complete, binsize, promoter_start, promoter_end, upset_plot, circos_plot, circos_use_promoters, in_regions):

    # Column names for loaded data
    anchor1_peak_name = ('Anchor1_Chr', 'Anchor1_Start', 'Anchor1_End', 'Peak1', 'Peak1_Chr', 'Peak1_Start', 'Peak1_End', 'Peak1_ID', 'Peak1_score')
    anchor2_peak_name = ('Anchor2_Chr', 'Anchor2_Start', 'Anchor2_End', 'Peak2', 'Peak2_Chr', 'Peak2_Start', 'Peak2_End', 'Peak2_ID', 'Peak2_score')

    # Load interaction centered peak overlaps 2and annotated 2D-bed
    anchor1_peaks = pd.read_table(anchor_1_peak_collect, index_col=3,names=anchor1_peak_name).sort_index().drop_duplicates()
    anchor2_peaks = pd.read_table(anchor_2_peak_collect, index_col=3,names=anchor2_peak_name).sort_index().drop_duplicates()
    bed2D_anno = pd.read_table(bed2D_index_anno, index_col=1).sort_index().iloc[:,1:]

    # Create Peak columns (chr:start-end) for anchor 1 & 2
    anchor1_peaks["Peak1_ID"] = anchor1_peaks["Peak1_Chr"].map(str) +':'+ (anchor1_peaks["Peak1_Start"]-1).map(str) +'-'+ anchor1_peaks["Peak1_End"].map(str)
    anchor2_peaks["Peak2_ID"] = anchor2_peaks["Peak2_Chr"].map(str) +':'+ (anchor2_peaks["Peak2_Start"]-1).map(str) +'-'+ anchor2_peaks["Peak2_End"].map(str)

    # Merging anchor points

    if (in_regions !="Not_specified"):
        anchor1_peaks_anno =bed2D_anno.loc[:,['Entrez_ID_1', 'Gene_Name_1', 'Distance_to_TSS_1', 'TSS_1']].merge(anchor1_peaks.loc[:,['Anchor1_Chr', 'Anchor1_Start', 'Anchor1_End','Peak1','Peak1_ID', 'Peak1_score']], left_index=True, right_index=True, how = 'left')
        anchor2_peaks_anno =bed2D_anno.loc[:,['Entrez_ID_2', 'Gene_Name_2','Distance_to_TSS_2', "TSS_2"]].merge(anchor2_peaks.loc[:,['Anchor2_Chr', 'Anchor2_Start', 'Anchor2_End','Peak2','Peak2_ID', 'Peak2_score']], left_index=True, right_index=True, how = 'left').merge(bed2D_anno.loc[:,['Interaction_score']], left_index=True, right_index=True, how = 'left')
        anchors_peaks_anno = anchor1_peaks_anno.merge(anchor2_peaks_anno, left_index=True, right_index=True,how = 'outer')
        anchors_peaks_anno = anchors_peaks_anno.iloc[:,np.r_[4:7,0:4,7:10,14:17,10:14,17:21]]
        anchors_peaks_anno.rename(columns = {'Anchor1_Chr': 'chr1', 'Anchor1_Start': 's1','Anchor1_End': 'e1','Anchor2_Chr': 'chr2', 'Anchor2_Start': 's2','Anchor2_End': 'e2'}, inplace = True)
        anchors_peaks_anno = anchors_peaks_anno.dropna(subset=['chr1', 'chr2'])
        anchors_peaks_anno.index = anchors_peaks_anno["Peak1_ID"].map(str) +'_'+ (anchors_peaks_anno["Peak2_ID"]).map(str)
        anchors_peaks_anno.index.name = 'Interaction'
        anchors_peaks_anno= anchors_peaks_anno.reset_index().drop_duplicates().set_index('Interaction')

        anchors_peaks_anno.rename(columns = {'TSS_1': 'Is_Promoter_1', 'TSS_2': 'Is_Promoter_2'}, inplace = True)
        anchors_peaks_anno['Is_Promoter_1'] = np.where(anchors_peaks_anno['Peak1_score'] ==0, 1, 0)
        anchors_peaks_anno['Is_Promoter_2'] = np.where(anchors_peaks_anno['Peak2_score'] ==0, 1, 0)

    else:
        anchor1_peaks_anno =bed2D_anno.loc[:,['chr1', 's1', 'e1','Entrez_ID_1', 'Gene_Name_1', 'Distance_to_TSS_1', 'TSS_1']].merge(anchor1_peaks.loc[:,['Peak1','Peak1_ID', 'Peak1_score']], left_index=True, right_index=True, how = 'left')
        anchor2_peaks_anno =bed2D_anno.loc[:,['chr2', 's2', 'e2','Entrez_ID_2', 'Gene_Name_2','Distance_to_TSS_2', "TSS_2"]].merge(anchor2_peaks.loc[:,['Peak2','Peak2_ID', 'Peak2_score']], left_index=True, right_index=True, how = 'left').merge(bed2D_anno.loc[:,['Interaction_score']], left_index=True, right_index=True, how = 'left')
        anchors_peaks_anno = anchor1_peaks_anno.merge(anchor2_peaks_anno, left_index=True, right_index=True,how = 'outer')
        anchors_peaks_anno.index.name = 'Interaction'
        anchors_peaks_anno= anchors_peaks_anno.reset_index().drop_duplicates().set_index('Interaction')

        anchors_peaks_anno.rename(columns = {'TSS_1': 'Is_Promoter_1', 'TSS_2': 'Is_Promoter_2'}, inplace = True)
        anchors_peaks_anno['Is_Promoter_1'] = np.where(((anchors_peaks_anno['Distance_to_TSS_1'] <= 0) & (anchors_peaks_anno['Distance_to_TSS_1'] >= -(binsize/2+promoter_end))) | ((anchors_peaks_anno['Distance_to_TSS_1'] >= 0) & (anchors_peaks_anno['Distance_to_TSS_1'] <= (binsize/2+promoter_start))),1,0)
        anchors_peaks_anno['Is_Promoter_2'] = np.where(((anchors_peaks_anno['Distance_to_TSS_2'] <= 0) & (anchors_peaks_anno['Distance_to_TSS_2'] >= -(binsize/2+promoter_end))) | ((anchors_peaks_anno['Distance_to_TSS_2'] >= 0) & (anchors_peaks_anno['Distance_to_TSS_2'] <= (binsize/2+promoter_start))),1,0)

    # Creation and use of function for adding 2 columns for each factor (overlap in anchor 1/2) with 1 if overlap
    def peak_in_anchor_1(row):
        if row['Peak1'] == f :
            return 1
        else:
            return ''
    def peak_in_anchor_2(row):
        if row['Peak2'] == f :
            return 1
        else:
            return ''

    factor = pd.unique(anchors_peaks_anno[['Peak1', 'Peak2']].dropna().values.ravel('K'))
    factor = factor[ (factor != "REGIONS") & (factor != "ALL") ]

    for f in factor:
        anchors_peaks_anno[f + '_1'] = anchors_peaks_anno.apply (lambda row: peak_in_anchor_1(row), axis=1)
    for f in factor:
        anchors_peaks_anno[f + '_2'] = anchors_peaks_anno.apply (lambda row: peak_in_anchor_2(row), axis=1)

    anchors_peaks_anno.index.name = 'Interaction'
    anchors_peaks_anno_original = anchors_peaks_anno.copy(deep=True)

    # Creating dictionary with each factors as a key and associated df with interactions with factor overlap in at least one anchor point
    factor_dict={}
    for f in factor:
        factor_dict[f] = anchors_peaks_anno[(anchors_peaks_anno['Peak1'] == f) | (anchors_peaks_anno['Peak2'] == f)]
        factor_dict[f].loc[factor_dict[f].Peak1 !=f,['Peak1', 'Peak1_ID', 'Peak1_score']] = ''
        factor_dict[f].loc[factor_dict[f].Peak2 !=f,['Peak2', 'Peak2_ID', 'Peak2_score']] = ''
        if (circos_use_promoters == "true"):
            factor_dict[f] = factor_dict[f][(factor_dict[f]["Peak1_score"]==1) | (factor_dict[f]["Peak2_score"]==1)]
        factor_dict[f] = factor_dict[f].groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    # Saving annotated interactions files (all interaction and interactions with peak overlap)
    for f in factor:
        factor_dict[f].to_csv(str(f) + '_' + prefix + '_interactions.txt', index=False, sep='\t' )

    if (in_regions !="Not_specified" and circos_use_promoters == "true"):
        anchors_peaks_anno_promoters = anchors_peaks_anno.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))
        anchors_peaks_anno_promoters.to_csv(prefix + '_HOMER_annotated_interactions_with_peak_and_promoter.txt', index=True, sep='\t' )
        anchors_peaks_anno_filt = anchors_peaks_anno[(anchors_peaks_anno["Peak1_score"]==1) | (anchors_peaks_anno["Peak2_score"]==1)]
        anchors_peaks_anno_filt = anchors_peaks_anno_filt.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))
        anchors_peaks_anno_filt.to_csv(prefix + '_HOMER_annotated_interactions_with_peak_overlap.txt', index=True, sep='\t' )

    else:
        anchors_peaks_anno = anchors_peaks_anno.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))
        anchors_peaks_anno.to_csv(prefix + '_HOMER_annotated_interactions_with_peak_overlap.txt', index=True, sep='\t' )

    # Save files for Network
    if (network == 'true' or complete == 'true' or upset_plot =='true' or circos_plot =='true'):
        anchors_peaks_anno_original.to_csv(prefix + '_HOMER_annotated_interactions_with_peak_overlap_not_aggregated.txt', index=True, sep='\t' )


# RUN FUNCTION
interaction_annotation_multiple(anchor_1_peak_collect=args.ANCHOR1_PEAK,anchor_2_peak_collect=args.ANCHOR2_PEAK,bed2D_index_anno=args.BED2D, prefix=args.PREFIX, network=args.NETWORK, complete=args.COMPLETE, binsize=args.BINSIZE, promoter_start=args.PROMOTER_START, promoter_end=args.PROMOTER_END, upset_plot=args.UPSET_PLOT, circos_plot=args.CIRCOS_PLOT, circos_use_promoters=args.CIRCOS_USE_PROMOTERS, in_regions=args.IN_REGIONS)
