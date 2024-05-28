#!/usr/bin/env python

### PROCESS 7: PEAK_INTERACTION_BASED_ANNOTATION ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'ANNOTATION OF PEAKS BY PROXIMITY AND INTERACTION-BASED ANNOTATION'
Epilog = """Example usage: peak_annotaiton.py <PEAK_ANCHOR1> <PEAK_ANCHOR2> <PEAK_ANNO> <BED2D> --prefix <PREFIX>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# General arguments
argParser.add_argument('PEAK_ANCHOR1', help="Annotated anchor1 regions.")
argParser.add_argument('PEAK_ANCHOR2', help="Annotated anchor2 regions.")
argParser.add_argument('PEAK_ANNO', help="HOMER annotated peak file.")
argParser.add_argument('BED2D', help="2D-bed interactions.")
argParser.add_argument('--promoter_pos', dest="PROMOTER_POS", help="Text file containing promoter positions from HOMER.")
argParser.add_argument('--peak_name', dest="PEAK_NAME", help="Name of peak.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--proximity_unannotated', dest="PROXIMITY_UNANNOTATED", help="Specifies if unannotated distal peaks should be annotated by proximity annotation (default: false)")
argParser.add_argument('--mode', dest="MODE", help="Define which mode to run the pipeline in. The options are basic (default), multiple or differential.", choices=['basic', 'multiple', 'differential'])
argParser.add_argument('--multiple_anno', dest="MULTIPLE_ANNO", help="Defines how to handle peaks annotated to more than one promoter. Options are keep (all annotations are kept with one row for each annotation), concentrate (the annotated peak file is concentrated to only include one row per peak but information about all annotations are kept) and qvalue (only one annotation per peak is kept. The annotation is decided by the interaction with the lowest qvalue). Default is: concentrate.", choices=['keep', 'concentrate', 'qvalue'])
argParser.add_argument('--promoter_start', dest="PROMOTER_START", help="Specifies the upstream of TSS considered as a promoter (default: 2500).", type=int)
argParser.add_argument('--promoter_end', dest="PROMOTER_END", help="Specifies the downstream of TSS considered as a promoter (default: 2500).", type=int)
argParser.add_argument('--binsize', dest="BINSIZE", help="Specifies interaction binsize (default: 5000)", type=int)
argParser.add_argument('--close_peak_type', dest="CLOSE_PEAK_TYPE", help="Specifies how to handle interactions close to peaks. Options are bin (based on number of bins) or distance (distance from peaks start/end to bin). Default: bin.", choices=['overlap','bin', 'distance'])
argParser.add_argument('--close_peak_distance', dest="CLOSE_PEAK_DISTANCE", help="Specify distance for peak annoation with close interaction. If --close_peak_type is bin (default) the option specifies number of bins +/- overlapping bin and if close_peak_type is distance it specifies distance from peak start/end to bin. Default: 1.", type=int)
argParser.add_argument('--skip_promoter_promoter', dest="SKIP_PROMOTER_PROMOTER", help="SPecifies with interaction-based annotation of peaks located in promoter regions should be skiped (default: false).", choices=['true', 'false'])
argParser.add_argument('--interaction_threshold', dest="INTERACTION_THRESHOLD", help="Lower interaction distance threshold, regions with a distance to the closest TSS < interaction_threshold will be proximity annotated (default: 2*binsize).", type=int)
argParser.add_argument('--close_promoter_type', dest="CLOSE_PROMOTER_TYPE", help="Specifies how to handle interactions close to promoter. Options are bin (based on number of bins) or distance (distance from TSS to bin). Default: bin.", choices=['overlap', 'bin', 'distance'])
argParser.add_argument('--close_promoter_distance_start', dest="CLOSE_PROMOTER_DISTANCE_START", help="Specify distance for interaction close to but not overlapping TSS, upstreams of TSS. Default: 2500.", type=int)
argParser.add_argument('--close_promoter_distance_end', dest="CLOSE_PROMOTER_DISTANCE_END", help="Specify distance for interaction close to but not overlapping TSS, downstreams of TSS lf TSS. Default: 2500.", type=int)
argParser.add_argument('--close_promoter_bin', dest="CLOSE_PROMOTER_BIN", help="Specify distance for interaction close to but not overlapping TSS, upstream lf TSS if --close_promoter_type is bin. The option specifies number of bins +/- overlapping bin and if close_peaks_type is distance it specifies distance from TSS to bin. Default: 1.", type=int)
argParser.add_argument('--filter_close', dest="FILTER_CLOSE", help="Depending on the close peak/promoter options, the same peak can be annotated to the same gene from interactions in neighboring bins. This options specifies how to handle this: no_filter (no filtering), peak (filter on peak_bin_distance firs and than TSS_bin_distance), tss (filter on TSS_bin_distance firs and than peak_bin_distance) or sum (sum of absolite values of peak_bin_distance and TSS_bin_distance). Default: no_filter",choices=['peak','tss', 'sum', 'no_filter'])

#Multiple mode specific arguments
argParser.add_argument('--circos_use_promoters', dest="CIRCOS_USE_PROMOTERS", help="Specifies if TF overlap in promoters (defined based on promoter_start/end) should be used in circos plot in multiple mode when regions are specified.", choices=['true', 'false'])

# Differntial mode specific arguments
argParser.add_argument('--peak_differential', dest='PEAK_DIFFERENTIAL', help="Path to textfile that contain log2FC and adjusted p-value from differential analysis. The 1st column should contain peakID matching the peakID in the 4th column of the input bed file. Standard DESeq2 output is expected (with log2FC in the 3rd column and padj in the 9th column), but other formats are accepted as well is the column corresponding to log2FC and padj are specified with the arguments --log2FC_column and --padj column.")
argParser.add_argument('--log2FC_column', dest="LOG2FC_COLUMN", help="Log2FC column for differential peaks.", type=int)
argParser.add_argument('--padj_column', dest="PADJ_COLUMN", help="Padj column for differential peaks.", type=int)
argParser.add_argument('--log2FC', dest="LOG2FC", help="Log2FC threshold for differential peaks.", type=float)
argParser.add_argument('--padj', dest="PADJ", help="Padj threshold for differential peaks.", type=float)
argParser.add_argument('--skip_expression', dest="SKIP_EXPRESSION", help="Specify this argument if no --expression file is provided.", choices=['true', 'false'])

args = argParser.parse_args()

# DEFINE FUNCTION
def peak_annotation(peak_anno_anchor1,peak_anno_anchor2,peak_anno, bed2D_index_anno, promoter_pos, peak_name, prefix, proximity_unannotated, mode, multiple_anno, promoter_start, promoter_end, binsize, close_peak_type, close_peak_distance, skip_promoter_promoter, interaction_threshold, close_promoter_type, close_promoter_distance_start, close_promoter_distance_end, close_promoter_bin, filter_close, circos_use_promoters, peak_differential, log2FC_column, padj_column, log2FC, padj, skip_expression):

    # Column names for loaded data
    peak_anchor1_name = ('peak_chr', 'peak_start', 'peak_end','Peak_score', 'anchor1_chr', 'anchor1_start', 'anchor1_end', 'anchor1_id')
    peak_anchor2_name = ('peak_chr', 'peak_start', 'peak_end',  'Peak_score', 'anchor2_chr', 'anchor2_start', 'anchor2_end', 'anchor2_id')

    # Load peak overlap for anchor 1 & 2, as well as annotated peak & 2D-bed files
    peak_anchor1 = pd.read_table(peak_anno_anchor1, index_col=3, names=peak_anchor1_name).sort_index()
    peak_anchor2 = pd.read_table(peak_anno_anchor2, index_col=3, names=peak_anchor2_name).sort_index()

    #if (circos_use_promoters == "true" and peak_name!="ALL"):
    peak_anchor1 = peak_anchor1[peak_anchor1['Peak_score'] == 1]
    peak_anchor2 = peak_anchor2[peak_anchor2['Peak_score'] == 1]

    peak_anno = pd.read_table(peak_anno,index_col=0).sort_index()
    bed2D_anno = pd.read_table(bed2D_index_anno, index_col=1).sort_index().iloc[:,1:]
    bed2D_anno.rename(columns={bed2D_anno.columns[0]: 'chr1', bed2D_anno.columns[1]: 's1', bed2D_anno.columns[2]: 'e1', bed2D_anno.columns[3]: 'chr2', bed2D_anno.columns[4]: 's2', bed2D_anno.columns[5]: 'e2'}, inplace =True)
    promoter_pos = pd.read_table(promoter_pos, names=("TSS_chr", "TSS_start", "TSS_end", "TSS_strand"))

    # Match peaks with interactions annotations for overlap with anchor point 1 & 2 respectily - Then merge
    Peak_overlap_1 =peak_anno.loc[:,['Chr','Start','End', 'Peak Score', 'Distance to TSS','Entrez ID','Nearest Refseq','Nearest Ensembl','Gene Name']].merge(peak_anchor1.iloc[:,7:], left_index=True, right_index=True, how = 'outer')\
      .merge(bed2D_anno.loc[:,['chr1', 's1', 'e1', 'chr2', 's2', 'e2', 'Interaction_score','Distance_to_TSS_2', 'Entrez_ID_2', 'Nearest_Refseq_2', 'Nearest_Ensembl_2', 'Gene_Name_2','TSS_1', 'TSS_2']], left_on='anchor1_id', right_index=True, how = 'left').drop_duplicates()
    Peak_overlap_1['overlap'] = 1
    Peak_overlap_1.columns = ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS', 'EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'InteractionID',  'Anchor_Overlap_Chr', 'Anchor_Overlap_Start', 'Anchor_Overlap_End', 'Anchor_Interaction_Chr', 'Anchor_Interaction_Start', 'Anchor_Interaction_End', 'Interaction_score', 'TSS_bin_distance', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Anchor_Overlap_TSS', 'Promoter_bin', 'Anchor_Overlap']
    Peak_overlap_2 =peak_anno.loc[:,['Chr','Start','End', 'Peak Score', 'Distance to TSS','Entrez ID','Nearest Refseq','Nearest Ensembl','Gene Name']].merge(peak_anchor2.iloc[:,7:], left_index=True, right_index=True, how = 'outer')\
      .merge(bed2D_anno.loc[:,['chr2', 's2', 'e2', 'chr1', 's1', 'e1', 'Interaction_score','Distance_to_TSS_1', 'Entrez_ID_1', 'Nearest_Refseq_1', 'Nearest_Ensembl_1', 'Gene_Name_1','TSS_2', 'TSS_1']], left_on='anchor2_id', right_index=True, how = 'left').drop_duplicates()
    Peak_overlap_2['overlap'] = 2
    Peak_overlap_2.columns = ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS', 'EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'InteractionID',  'Anchor_Overlap_Chr', 'Anchor_Overlap_Start', 'Anchor_Overlap_End','Anchor_Interaction_Chr', 'Anchor_Interaction_Start', 'Anchor_Interaction_End', 'Interaction_score', 'TSS_bin_distance', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Anchor_Overlap_TSS', 'Promoter_bin', 'Anchor_Overlap']
    Peak_overlap_merge = pd.concat([Peak_overlap_1, Peak_overlap_2], axis=0).sort_index()

    # Create a new column that specify type of annotation for each peak: Promoter, proximal annotation (Homer) or interactuion-based annotation
    Peak_overlap_merge['Peak_type'] = np.where((Peak_overlap_merge['Distance_to_TSS'] >= -promoter_start) & (Peak_overlap_merge['Distance_to_TSS'] <= promoter_end), 'Promoter', (np.where(abs(Peak_overlap_merge['Distance_to_TSS']) <= interaction_threshold, 'Proximal', 'Distal')))

    # Extrating promoter and proximity annotated peak, adding Q_value column (for filtering) and renaming columns
    Proximal = Peak_overlap_merge.loc[Peak_overlap_merge['Peak_type'].isin(['Promoter','Proximal']),['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'Peak_type']].drop_duplicates()
    Proximal['Peak_bin_distance'], Proximal['Q-Peak_bin'], Proximal['Annotation_method'], Proximal['Interaction_score'], Proximal['Promoter_bin'], Proximal['TSS_bin_distance'] = [np.nan, np.nan, 'Proximal_anno', np.nan, np.nan, np.nan]
    Proximal.columns = ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID', 'Refseq','Ensembl', 'Gene', 'Peak_type', 'Peak_bin_distance', 'Peak_bin', 'Annotation_method', 'Interaction_score', 'Promoter_bin', 'TSS_bin_distance']

    # Extracting interaction annotated peaks
    Distal = Peak_overlap_merge.dropna(subset=['InteractionID'])
    # Create two new columns for close interaction anntation: peak_bin_distance (distance between peak center and bin (negative value means peak is upstream of interaction)) and bin (+/- number of bins from peak to bin with interaction)
    Distal['Peak_bin_distance'] = (Distal.Start+((Distal.End-Distal.Start)/2)) - (Distal.Anchor_Overlap_Start+((Distal.Anchor_Overlap_End-Distal.Anchor_Overlap_Start)/2))
    Distal['Peak_bin'] = round(Distal.Peak_bin_distance/binsize)
    Distal['Peak_bin_distance'] = np.where(Distal['Peak_bin_distance'] < -binsize/2, Distal['Peak_bin_distance'] + binsize/2, (np.where(Distal['Peak_bin_distance'] > binsize/2, Distal['Peak_bin_distance'] - binsize/2, 0)))

    #Filtering based on close peak/promoter distance
    if close_peak_type == 'overlap':
        Distal = Distal.loc[(abs(Distal['Peak_bin_distance'])-((Distal.End-Distal.Start)/2)) < 0,:]
    elif close_peak_type == 'bin':
        Distal = Distal.loc[~(abs(Distal['Peak_bin']) > close_peak_distance),:]
    elif close_peak_type == 'distance':
        Distal = Distal.loc[~(abs(Distal['Peak_bin_distance']) > close_peak_distance),:]

    if close_promoter_type == 'overlap':
        Distal = Distal.loc[(((Distal['TSS_bin_distance'] <= 0) & (Distal['TSS_bin_distance'] >= -(binsize/2+promoter_end))) | ((Distal['TSS_bin_distance'] >= 0) & (Distal['TSS_bin_distance'] <= (binsize/2+promoter_start)))),['Chr', 'Start', 'End', 'Peak_score','Distance_to_TSS', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Peak_type', 'Peak_bin_distance', 'Peak_bin','Interaction_score', 'Promoter_bin','TSS_bin_distance']].drop_duplicates()
    elif close_promoter_type == 'bin':
        Distal = Distal.loc[~(abs(Distal['Promoter_bin']) > close_promoter_bin),['Chr', 'Start', 'End', 'Peak_score','Distance_to_TSS', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Peak_type', 'Peak_bin_distance', 'Peak_bin','Interaction_score', 'Promoter_bin','TSS_bin_distance']].drop_duplicates()
    elif close_promoter_type == 'distance':
        Distal = Distal.loc[(((Distal['TSS_bin_distance'] <= 0) & (Distal['TSS_bin_distance'] >= -(binsize/2+close_promoter_distance_start))) | ((Distal['TSS_bin_distance'] >= 0) & (Distal['TSS_bin_distance'] <= (binsize/2+close_promoter_distance_end)))),['Chr', 'Start', 'End', 'Peak_score','Distance_to_TSS', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Peak_type', 'Peak_bin_distance', 'Peak_bin','Interaction_score', 'Promoter_bin','TSS_bin_distance']].drop_duplicates()

    Distal['Annotation_method'] = 'Interaction_anno'
    Distal = Distal.loc[:,['Chr', 'Start', 'End', 'Peak_score','Distance_to_TSS', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Peak_type', 'Peak_bin_distance', 'Peak_bin','Annotation_method', 'Interaction_score', 'Promoter_bin','TSS_bin_distance']]
    Distal.columns= ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID', 'Refseq','Ensembl', 'Gene', 'Peak_type', 'Peak_bin_distance', 'Peak_bin','Annotation_method', 'Interaction_score', 'Promoter_bin','TSS_bin_distance']

    # Merge proximity and PLAC-seq annotated peaks
    Proximal_Distal = pd.concat([Proximal, Distal]).sort_index().rename_axis('Peak')

    # Annotate unannotated distal peaks by proximity annotation
    if proximity_unannotated == 'true':
        # Extracting unannotated distal peaks (not overlapping 2D-bed)
        Unannotated = Peak_overlap_merge.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'Peak_type']][~Peak_overlap_merge.index.isin(Proximal_Distal.index)].drop_duplicates()
        Unannotated['Peak_bin_distance'], Unannotated['Peak_bin'], Unannotated['Annotation_method'], Unannotated['Interaction_score'], Unannotated['Promoter_bin'], Unannotated['TSS_bin_distance'] = [np.NaN, np.NaN,'Proximal_anno', np.NaN, np.NaN, np.NaN]
        Unannotated.columns=['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID', 'Refseq','Ensembl', 'Gene', 'Peak_type','Peak_bin_distance', 'Peak_bin', 'Annotation_method', 'Interaction_score', 'Promoter_bin','TSS_bin_distance']
        Proximal_Distal = pd.concat([Proximal_Distal, Unannotated]).sort_index().rename_axis('Peak')

    Proximal_Distal['Start'] = Proximal_Distal['Start']-1
    if (skip_promoter_promoter=='true'):
        #Proximal_Distal = Proximal_Distal.drop(Proximal_Distal[(Proximal_Distal.Peak_type =='Promoter') & (Proximal_Distal.Annotation_method =='Interaction_anno')].index)
        Proximal_Distal = Proximal_Distal.loc[~((Proximal_Distal.Peak_type =='Promoter') & (Proximal_Distal.Annotation_method =='Interaction_anno')),:]

    # Assigning distance to TSS for the annotation (not closest gene) using HOMER TSS position
    promoter_pos['TSS'] = np.where(promoter_pos['TSS_strand'] == 1, promoter_pos.TSS_end-promoter_end, promoter_pos.TSS_end-promoter_start)
    Proximal_Distal=Proximal_Distal.merge(promoter_pos.loc[:,["TSS", "TSS_strand"]], left_on='Refseq', right_index=True, how = 'left')
    Proximal_Distal["Distance_to_TSS_int"] = np.where(Proximal_Distal['TSS_strand'] ==0, ((Proximal_Distal.Start+Proximal_Distal.End)/2)-Proximal_Distal.TSS, np.where(Proximal_Distal['TSS_strand'] ==1,-(((Proximal_Distal.Start+Proximal_Distal.End)/2)-Proximal_Distal.TSS),np.nan))
    Proximal_Distal.loc[Proximal_Distal['Annotation_method'] =='Interaction_anno', 'Distance_to_TSS'] = Proximal_Distal.loc[Proximal_Distal['Annotation_method'] =='Interaction_anno', 'Distance_to_TSS_int']
    Proximal_Distal=Proximal_Distal.drop(columns=['TSS',"TSS_strand", "Distance_to_TSS_int"])

    # Filer for multiple annoation of the same peak to the same gene (result of close peak/bin options)
    if filter_close == "peak":
        Proximal_Distal =Proximal_Distal.sort_values(['Peak_bin_distance', 'TSS_bin_distance'], key=abs)
        Proximal_Distal = Proximal_Distal.drop_duplicates(subset=['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS', 'EntrezID', 'Refseq', 'Ensembl', 'Gene', 'Peak_type', 'Annotation_method', 'Interaction_score'], keep="first")
    elif filter_close == "tss":
        Proximal_Distal =Proximal_Distal.sort_values(['TSS_bin_distance', 'Peak_bin_distance'], key=abs)
        Proximal_Distal = Proximal_Distal.drop_duplicates(subset=['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS', 'EntrezID', 'Refseq', 'Ensembl', 'Gene', 'Peak_type', 'Annotation_method', 'Interaction_score'], keep="first")
    elif filter_close == "sum":
        Proximal_Distal["sum"]= abs(Proximal_Distal["Peak_bin_distance"])+abs(Proximal_Distal["TSS_bin_distance"])
        Proximal_Distal =Proximal_Distal.sort_values('sum')
        Proximal_Distal = Proximal_Distal.drop_duplicates(subset=['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS', 'EntrezID', 'Refseq', 'Ensembl', 'Gene', 'Peak_type', 'Annotation_method', 'Interaction_score'], keep="first")
        Proximal_Distal=Proximal_Distal.drop(columns=['sum'])

    if peak_name == 'ALL':
        peak_anno_2 = peak_anno.iloc[:,0:5].rename_axis('factor').reset_index().replace({'factor': r'-[0-9]*$'}, {'factor': ''}, regex=True).drop_duplicates()
        peak_anno_2["peakID"] = peak_anno_2["Chr"] + ":" + (peak_anno_2["Start"]-1).astype(str) + "-" + peak_anno_2["End"].astype(str)
        peak_anno_2 = peak_anno_2.loc[:,['peakID', 'factor', 'Peak Score']].pivot( index='peakID', columns='factor', values='Peak Score')
        peak_anno_2['Peak_Score'] = peak_anno_2.sum(axis=1)

        Proximal_Distal.index = Proximal_Distal["Chr"] + ":" + Proximal_Distal["Start"].astype(str) + "-" + Proximal_Distal["End"].astype(str)
        Proximal_Distal = Proximal_Distal.merge(peak_anno_2, left_index=True, right_index=True, how = 'left')
        Proximal_Distal = Proximal_Distal.iloc[:, np.r_[0,1,2,len(Proximal_Distal.columns)-1,16:len(Proximal_Distal.columns)-1, 4:16]]

    # Organizing and saving annotated files/genelsits
    # Basic/Multiple mode: Handling of peaks annotating to several genes
    if (mode=='basic' or mode=='multiple'):
        if multiple_anno == 'keep':
            Proximal_Distal = Proximal_Distal
            Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()
        elif multiple_anno == 'concentrate':
            Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()
            Proximal_Distal_merge = Proximal_Distal.groupby('Peak')[['Distance_to_TSS','EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type','Peak_bin_distance', 'Peak_bin', 'Annotation_method', 'Interaction_score', 'Promoter_bin','TSS_bin_distance']].agg(lambda x: ', '.join(list(x.astype(str))))
            Proximal_Distal = Proximal_Distal.loc[:,['Chr', 'Start', 'End', 'Peak_score']].merge(Proximal_Distal_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
        elif multiple_anno == 'one_annotation':
            Proximal_Distal_promoter = Proximal_Distal.loc[(Proximal_Distal["Peak_type"] == 'Promoter') & (Proximal_Distal["Annotation_method"] == 'Proximal_anno')]
            Proximal_Distal['abs_bin'] = abs(Proximal_Distal['Peak_bin'])
            Proximal_Distal_remaining =  Proximal_Distal[~Proximal_Distal.index.isin(Proximal_Distal_promoter.index)].sort_values(["abs_bin", "Q-Value"], ascending = (True, True)).reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
            Proximal_Distal = pd.concat([Proximal_Distal_promoter, Proximal_Distal_remaining]).sort_index().drop_duplicates()
            Proximal_Distal=Proximal_Distal.drop(columns=['abs_bin'])
            Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()

        Proximal_Distal =Proximal_Distal.dropna(subset=['EntrezID','Refseq','Ensembl', 'Gene'], how='all')
        Proximal_Distal.to_csv(peak_name + '_' + prefix + '_annotated.txt', index=False, sep='\t' )
        pd.DataFrame(Genelist).to_csv(peak_name + '_' + prefix + '_annotated_genelist.txt', index=False, header=False,sep='\t' )

    # Differntial mode: Creting peak annotation/genelists for UP/DOWN peak_start
    if mode == 'differential':
      peak_differntial = pd.read_table(peak_differential, index_col=0).sort_index()
      peak_differntial = peak_differntial.iloc[:, [log2FC_column-2, padj_column-2]]
      peak_differntial.columns = ['log2FC', 'padj']
      Proximal_Distal_differential=Proximal_Distal.merge(peak_differntial, left_index=True, right_index=True, how = 'left')
      Proximal_Distal_differential.index.name = 'Peak'
      Proximal_Distal_up=Proximal_Distal_differential[(Proximal_Distal_differential.padj <= padj) & (Proximal_Distal_differential.log2FC >=log2FC)]
      Proximal_Distal_down=Proximal_Distal_differential[(Proximal_Distal_differential.padj <= padj) & (Proximal_Distal_differential.log2FC <= -log2FC)]

      if skip_expression == 'false':
        Proximal_Distal_differential_for_differntial_expression = Proximal_Distal_differential.copy(deep=True)
        Proximal_Distal_differential_for_differntial_expression.to_csv(peak_name + '_' + prefix + '_annotated_for_differential_expression.txt', index=True, sep='\t' )

      # Differnetial mode: Handling of peaks annotating to several genes
      if multiple_anno == 'keep':
        Genelist = Proximal_Distal_differential.loc[:,'Gene'].unique().tolist()
        Genelist_up = Proximal_Distal_up.loc[:,'Gene'].unique().tolist()
        Genelist_down = Proximal_Distal_down.loc[:,'Gene'].unique().tolist()
      elif multiple_anno == 'concentrate':
        Genelist = Proximal_Distal_differential.loc[:,'Gene'].unique().tolist()
        Genelist_up = Proximal_Distal_up.loc[:,'Gene'].unique().tolist()
        Genelist_down = Proximal_Distal_down.loc[:,'Gene'].unique().tolist()
        Proximal_Distal_differential_merge = Proximal_Distal_differential.groupby('Peak')[['EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type','Peak_bin_distance', 'Peak_bin','Annotation_method', 'Interaction_score', 'Promoter_bin','TSS_bin_distance']].agg(lambda x: ', '.join(list(x.astype(str))))
        Proximal_Distal_differential = Proximal_Distal_differential.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS']].merge(Proximal_Distal_differential_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
        Proximal_Distal_up_merge = Proximal_Distal_up.groupby('Peak')[['EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type', 'Peak_bin_distance', 'Peak_bin', 'Annotation_method', 'Interaction_score', 'Promoter_bin','TSS_bin_distance']].agg(lambda x: ', '.join(list(x.astype(str))))
        Proximal_Distal_up = Proximal_Distal_up.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS']].merge(Proximal_Distal_up_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
        Proximal_Distal_down_merge = Proximal_Distal_down.groupby('Peak')[['EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type','Peak_bin_distance', 'Peak_bin', 'Annotation_method', 'Interaction_score','Promoter_bin','TSS_bin_distance']].agg(lambda x: ', '.join(list(x.astype(str))))
        Proximal_Distal_down = Proximal_Distal_down.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS']].merge(Proximal_Distal_down_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
      elif multiple_anno == 'one_annotation':
        Proximal_Distal_differential_promoter = Proximal_Distal_differential.loc[(Proximal_Distal_differential["Peak_type"] == 'Promoter') & (Proximal_Distal_differential["Annotation_method"] == 'Proximal_anno')]
        Proximal_Distal_differential['abs_bin'] = abs(Proximal_Distal_differential['Peak_bin'])
        Proximal_Distal_differential_remaining =  Proximal_Distal_differential[~Proximal_Distal_differential.index.isin(Proximal_Distal_differential_promoter.index)].sort_values(["abs_bin", "Q-Value"], ascending = (True, True)).reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
        Proximal_Distal_differential = pd.concat([Proximal_Distal_differential_promoter, Proximal_Distal_differential_remaining]).sort_index().drop_duplicates()
        Proximal_Distal_differential=Proximal_Distal_differential.drop(columns=['abs_bin'])
        Genelist = Proximal_Distal_differential.loc[:,'Gene'].unique().tolist()
        Proximal_Distal_up_promoter = Proximal_Distal_up.loc[(Proximal_Distal_up["Peak_type"] == 'Promoter') & (Proximal_Distal_up["Annotation_method"] == 'Proximal_anno')]
        Proximal_Distal_up['abs_bin'] = abs(Proximal_Distal_up['Peak_bin'])
        Proximal_Distal_up_remaining =  Proximal_Distal_up[~Proximal_Distal_up.index.isin(Proximal_Distal_up_promoter.index)].sort_values(["abs_bin", "Q-Value"], ascending = (True, True)).reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
        Proximal_Distal_up = pd.concat([Proximal_Distal_up_promoter, Proximal_Distal_up_remaining]).sort_index().drop_duplicates()
        Proximal_Distal_up=Proximal_Distal_up.drop(columns=['abs_bin'])
        Genelist_up = Proximal_Distal_up.loc[:,'Gene'].unique().tolist()
        Proximal_Distal_down_promoter = Proximal_Distal_down.loc[(Proximal_Distal_down["Peak_type"] == 'Promoter') & (Proximal_Distal_down["Annotation_method"] == 'Proximal_anno')]
        Proximal_Distal_down['abs_bin'] = abs(Proximal_Distal_down['Peak_bin'])
        Proximal_Distal_down_remaining =  Proximal_Distal_down[~Proximal_Distal_down.index.isin(Proximal_Distal_down_promoter.index)].sort_values(["abs_bin", "Q-Value"], ascending = (True, True)).reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
        Proximal_Distal_down = pd.concat([Proximal_Distal_down_promoter, Proximal_Distal_down_remaining]).sort_index().drop_duplicates()
        Proximal_Distal_down=Proximal_Distal_down.drop(columns=['abs_bin'])
        Genelist_down = Proximal_Distal_down.loc[:,'Gene'].unique().tolist()

      Proximal_Distal_differential =Proximal_Distal_differential.dropna(subset=['EntrezID','Refseq','Ensembl', 'Gene'], how='all')
      Proximal_Distal_differential.to_csv(peak_name + '_'+prefix + '_annotated.txt', index=False, sep='\t' )
      Proximal_Distal_up =Proximal_Distal_up.dropna(subset=['EntrezID','Refseq','Ensembl', 'Gene'], how='all')
      Proximal_Distal_up.to_csv(peak_name + '_' + prefix + '_annotated_up.txt', index=False, sep='\t' )
      Proximal_Distal_down =Proximal_Distal_down.dropna(subset=['EntrezID','Refseq','Ensembl', 'Gene'], how='all')
      Proximal_Distal_down.to_csv(peak_name + '_' + prefix + '_annotated_down.txt', index=False, sep='\t' )
      pd.DataFrame(Genelist).to_csv(peak_name + '_' + prefix + '_annotated_genelist.txt', index=False, header=False,sep='\t' )
      pd.DataFrame(Genelist_up).to_csv(peak_name + '_' + prefix + '_annotated_genelist_up.txt', index=False, header=False,sep='\t' )
      pd.DataFrame(Genelist_down).to_csv(peak_name + '_' + prefix + '_annotated_genelist_down.txt', index=False, header=False,sep='\t' )


# RUN FUNCTION
peak_annotation(peak_anno_anchor1=args.PEAK_ANCHOR1,peak_anno_anchor2=args.PEAK_ANCHOR2,peak_anno=args.PEAK_ANNO,bed2D_index_anno=args.BED2D, promoter_pos=args.PROMOTER_POS, peak_name=args.PEAK_NAME, prefix=args.PREFIX, proximity_unannotated=args.PROXIMITY_UNANNOTATED, mode=args.MODE, multiple_anno=args.MULTIPLE_ANNO, promoter_start=args.PROMOTER_START, promoter_end=args.PROMOTER_END, binsize=args.BINSIZE, close_peak_type=args.CLOSE_PEAK_TYPE, close_peak_distance=args.CLOSE_PEAK_DISTANCE, skip_promoter_promoter= args.SKIP_PROMOTER_PROMOTER, interaction_threshold=args.INTERACTION_THRESHOLD,  close_promoter_type=args.CLOSE_PROMOTER_TYPE, close_promoter_distance_start=args.CLOSE_PROMOTER_DISTANCE_START, close_promoter_distance_end=args.CLOSE_PROMOTER_DISTANCE_END, close_promoter_bin=args.CLOSE_PROMOTER_BIN, filter_close=args.FILTER_CLOSE, circos_use_promoters=args.CIRCOS_USE_PROMOTERS, peak_differential=args.PEAK_DIFFERENTIAL, log2FC_column=args.LOG2FC_COLUMN, padj_column=args.PADJ_COLUMN, log2FC=args.LOG2FC, padj=args.PADJ, skip_expression=args.SKIP_EXPRESSION)
