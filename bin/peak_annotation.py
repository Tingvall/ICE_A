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
argParser.add_argument('--peak_name', dest="PEAK_NAME", help="Name of peak.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--proximity_unannotated', dest="PROXIMITY_UNANNOTATED", help="Specifies if unannotated distal peaks should be annotated by proximity annotation (default: false)")
argParser.add_argument('--mode', dest="MODE", help="Define which mode to run the pipeline in. The options are basic (default), multiple or differential.", choices=['basic', 'multiple', 'differential'])
argParser.add_argument('--multiple_anno', dest="MULTIPLE_ANNO", help="Defines how to handle peaks annotated to more than one promoter. Options are keep (all annotations are kept with one row for each annotation), concentrate (the annotated peak file is concentrated to only include one row per peak but information about all annotations are kept) and qvalue (only one annotation per peak is kept. The annotation is decided by the interaction with the lowest qvalue). Default is: concentrate.", choices=['keep', 'concentrate', 'qvalue'])
argParser.add_argument('--promoter_start', dest="PROMOTER_START", help="Specifies the upstream of TSS considered as a promoter (default: 2500).", type=int)
argParser.add_argument('--promoter_end', dest="PROMOTER_END", help="Specifies the downstream of TSS considered as a promoter (default: 2500).", type=int)
argParser.add_argument('--skip_promoter_promoter', dest="SKIP_PROMOTER_PROMOTER", help="SPecifies with interaction-based annotation of peaks located in promoter regions should be skiped (default: false).", choices=['true', 'false'])
argParser.add_argument('--interaction_threshold', dest="INTERACTION_THRESHOLD", help="Lower interaction distance threshold, regions with a distance to the closest TSS < interaction_threshold will be proximity annotated (default: 2*binsize).", type=int)

# Differntial mode specific arguments
argParser.add_argument('--peak_differential', dest='PEAK_DIFFERENTIAL', help="Path to textfile that contain log2FC and adjusted p-value from differential analysis. The 1st column should contain peakID matching the peakID in the 4th column of the input bed file. Standard DESeq2 output is expected (with log2FC in the 3rd column and padj in the 9th column), but other formats are accepted as well is the column corresponding to log2FC and padj are specified with the arguments --log2FC_column and --padj column.")
argParser.add_argument('--log2FC_column', dest="LOG2FC_COLUMN", help="Log2FC column for differential peaks.", type=int)
argParser.add_argument('--padj_column', dest="PADJ_COLUMN", help="Padj column for differential peaks.", type=int)
argParser.add_argument('--log2FC', dest="LOG2FC", help="Log2FC threshold for differential peaks.", type=float)
argParser.add_argument('--padj', dest="PADJ", help="Padj threshold for differential peaks.", type=float)
argParser.add_argument('--skip_expression', dest="SKIP_EXPRESSION", help="Specify this argument if no --expression file is provided.", choices=['true', 'false'])

args = argParser.parse_args()

# DEFINE FUNCTION
def peak_annotation(peak_anno_anchor1,peak_anno_anchor2,peak_anno, bed2D_index_anno, peak_name, prefix, proximity_unannotated, mode, multiple_anno, promoter_start, promoter_end, skip_promoter_promoter, interaction_threshold, peak_differential, log2FC_column, padj_column, log2FC, padj, skip_expression):
    # Column names for loaded data
    peak_anchor1_name = ('peak_chr', 'peak_start', 'peak_end','Peak_score', 'anchor1_chr', 'anchor1_start', 'anchor1_end', 'anchor1_id')
    peak_anchor2_name = ('peak_chr', 'peak_start', 'peak_end',  'Peak_score', 'anchor2_chr', 'anchor2_start', 'anchor2_end', 'anchor2_id')

    # Load peak overlap for anchor 1 & 2, as well as annotated peak & 2D-bed files
    peak_anchor1 = pd.read_table(peak_anno_anchor1, index_col=3, names=peak_anchor1_name).sort_index()
    peak_anchor2 = pd.read_table(peak_anno_anchor2, index_col=3, names=peak_anchor2_name).sort_index()
    peak_anno = pd.read_table(peak_anno,index_col=0).sort_index()
    bed2D_anno = pd.read_table(bed2D_index_anno, index_col=1).sort_index().iloc[:,1:]

    # Match peaks with interactions annotations for overlap with anchor point 1 & 2 respectily - Then merge
    Peak_overlap_1 =peak_anno.loc[:,['Chr','Start','End', 'Peak Score', 'Distance to TSS','Entrez ID','Nearest Refseq','Nearest Ensembl','Gene Name']].merge(peak_anchor1.iloc[:,7:], left_index=True, right_index=True, how = 'outer')\
      .merge(bed2D_anno.loc[:,['chr2', 's2', 'e2', 'cc', 'P-Value_Bias', 'Q-Value_Bias','Entrez_ID_2', 'Nearest_Refseq_2', 'Nearest_Ensembl_2', 'Gene_Name_2','TSS_1', 'TSS_2']], left_on='anchor1_id', right_index=True, how = 'left').drop_duplicates()
    Peak_overlap_1['overlap'] = 1
    Peak_overlap_1.columns = ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS', 'EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'InteractionID', 'Anchor_Interaction_Chr', 'Anchor_Interaction_Start', 'Anchor_Interaction_End', 'cc', 'P-Value','Q-Value', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Anchor_Overlap_TSS', 'Anchor_Interaction_TSS', 'Anchor_Overlap']
    Peak_overlap_2 =peak_anno.loc[:,['Chr','Start','End', 'Peak Score', 'Distance to TSS','Entrez ID','Nearest Refseq','Nearest Ensembl','Gene Name']].merge(peak_anchor2.iloc[:,7:], left_index=True, right_index=True, how = 'outer')\
      .merge(bed2D_anno.loc[:,['chr1', 's1', 'e1', 'cc', 'P-Value_Bias', 'Q-Value_Bias','Entrez_ID_1', 'Nearest_Refseq_1', 'Nearest_Ensembl_1', 'Gene_Name_1','TSS_2', 'TSS_1']], left_on='anchor2_id', right_index=True, how = 'left').drop_duplicates()
    Peak_overlap_2['overlap'] = 2
    Peak_overlap_2.columns = ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS', 'EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'InteractionID', 'Anchor_Interaction_Chr', 'Anchor_Interaction_Start', 'Anchor_Interaction_End', 'cc', 'P-Value','Q-Value', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Anchor_Overlap_TSS', 'Anchor_Interaction_TSS', 'Anchor_Overlap']
    Peak_overlap_merge = pd.concat([Peak_overlap_1, Peak_overlap_2], axis=0).sort_index()

    # Create a new column that specify type of annotation for each peak: Promoter, proximal annotation (Homer) or PLAC-seq based annotation
    Peak_overlap_merge['Peak_type'] = np.where((Peak_overlap_merge['Distance_to_TSS'] >= -promoter_start) & (Peak_overlap_merge['Distance_to_TSS'] <= promoter_end), 'Promoter', (np.where(abs(Peak_overlap_merge['Distance_to_TSS']) <= interaction_threshold, 'Proximal', 'Distal')))

    # Extrating promoter and proximity annotated peak, adding Q_value column (for filtering) and renaming columns
    Proximal = Peak_overlap_merge.loc[Peak_overlap_merge['Peak_type'].isin(['Promoter','Proximal']),['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'Peak_type']].drop_duplicates()
    Proximal['Q-value'] = np.nan
    Proximal['Annotation_method'] = 'Proximal_anno'
    Proximal.columns = ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID', 'Refseq','Ensembl', 'Gene', 'Peak_type', 'Q-value', 'Annotation_method']

    # Extracting interaction annotated peaks
    Distal = Peak_overlap_merge.dropna(subset=['InteractionID'])
    Distal = Distal.loc[Distal['Anchor_Interaction_TSS'] == 1,['Chr', 'Start', 'End', 'Peak_score','Distance_to_TSS', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Peak_type', 'Q-Value']].drop_duplicates()
    Distal['Annotation_method'] = 'Interaction_anno'
    Distal.columns = ['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID', 'Refseq','Ensembl', 'Gene', 'Peak_type', 'Q-value', 'Annotation_method']

    # Merge proximity and PLAC-seq annotated peaks
    Proximal_Distal = pd.concat([Proximal, Distal]).sort_index().rename_axis('Peak')

    # Annotate unannotated distal peaks by proximity annotation
    if proximity_unannotated == 'true':
        # Extracting unannotated distal peaks (not overlapping 2D-bed)
        Unannotated = Peak_overlap_merge.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'Peak_type']][~Peak_overlap_merge.index.isin(Proximal_Distal.index)].drop_duplicates()
        Unannotated['Q-value'], Unannotated['Annotation_method'] = [np.NaN, 'Proximal_anno']
        Unannotated.columns=['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS','EntrezID', 'Refseq','Ensembl', 'Gene', 'Peak_type', 'Q-value', 'Annotation_method']
        Proximal_Distal = pd.concat([Proximal_Distal, Unannotated]).sort_index().rename_axis('Peak')

    Proximal_Distal['Start'] = Proximal_Distal['Start']-1
    if (skip_promoter_promoter=='true'):
        Proximal_Distal = Proximal_Distal.drop(Proximal_Distal[(Proximal_Distal.Peak_type =='Promoter') & (Proximal_Distal.Annotation_method =='Interaction_anno')].index)

    # Organizing and saving annotated files/genelsits
    # Basic/Multiple mode: Handling of peaks annotating to several genes
    if (mode=='basic' or mode=='multiple'):
        if multiple_anno == 'keep':
            Proximal_Distal = Proximal_Distal
            Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()
        elif multiple_anno == 'concentrate':
            Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()
            Proximal_Distal_merge = Proximal_Distal.groupby('Peak')[['EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type', 'Q-value', 'Annotation_method']].agg(lambda x: ', '.join(list(x.astype(str))))
            Proximal_Distal = Proximal_Distal.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS']].merge(Proximal_Distal_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
        elif multiple_anno == 'one_annotation':
            Proximal_Distal_promoter = Proximal_Distal.loc[(Proximal_Distal["Peak_type"] == 'Promoter') & (Proximal_Distal["Annotation_method"] == 'Proximal_anno')]
            Proximal_Distal_remaining =  Proximal_Distal[~Proximal_Distal.index.isin(Proximal_Distal_promoter.index)].sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
            Proximal_Distal = pd.concat([Proximal_Distal_promoter, Proximal_Distal_remaining]).sort_index().drop_duplicates()
            Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()

      Proximal_Distal.to_csv(peak_name + '_' + prefix + '_annotated.txt', index=False, sep='\t' )
      pd.DataFrame(Genelist).to_csv(peak_name + '_' + prefix + '_annotated_genelist.txt', index=False, header=False,sep='\t' )

    # Differntial mode: Creting peak annotation/genelists for UP/DOWN peak_start
    if mode == 'differential':
      peak_differntial = pd.read_table("peak_differential", index_col=0).sort_index()
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
        Proximal_Distal_differential_merge = Proximal_Distal_differential.groupby('Peak')[['EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type', 'Q-value', 'Annotation_method']].agg(lambda x: ', '.join(list(x.astype(str))))
        Proximal_Distal_differential = Proximal_Distal_differential.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS']].merge(Proximal_Distal_differential_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
        Proximal_Distal_up_merge = Proximal_Distal_up.groupby('Peak')[['EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type', 'Q-value', 'Annotation_method']].agg(lambda x: ', '.join(list(x.astype(str))))
        Proximal_Distal_up = Proximal_Distal_up.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS']].merge(Proximal_Distal_up_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
        Proximal_Distal_down_merge = Proximal_Distal_down.groupby('Peak')[['EntrezID','Refseq','Ensembl', 'Gene', 'Peak_type', 'Q-value', 'Annotation_method']].agg(lambda x: ', '.join(list(x.astype(str))))
        Proximal_Distal_down = Proximal_Distal_down.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'Distance_to_TSS']].merge(Proximal_Distal_down_merge, left_index=True, right_index=True, how='outer').drop_duplicates()
    elif multiple_anno == 'one_annotation':
        Proximal_Distal_differential_promoter = Proximal_Distal_differential.loc[(Proximal_Distal_differential["Peak_type"] == 'Promoter') & (Proximal_Distal_differential["Annotation_method"] == 'Proximal_anno')]
        Proximal_Distal_differential_remaining =  Proximal_Distal_differential[~Proximal_Distal_differential.index.isin(Proximal_Distal_differential_promoter.index)].sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
        Proximal_Distal_differential = pd.concat([Proximal_Distal_differential_promoter, Proximal_Distal_differential_remaining]).sort_index().drop_duplicates()
        Proximal_Distal_up_promoter = Proximal_Distal_up.loc[(Proximal_Distal_up["Peak_type"] == 'Promoter') & (Proximal_Distal_up["Annotation_method"] == 'Proximal_anno')]
        Proximal_Distal_up_remaining =  Proximal_Distal_up[~Proximal_Distal_up.index.isin(Proximal_Distal_up_promoter.index)].sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
        Proximal_Distal_up = pd.concat([Proximal_Distal_up_promoter, Proximal_Distal_up_remaining]).sort_index().drop_duplicates()
        Proximal_Distal_down_promoter = Proximal_Distal_down.loc[(Proximal_Distal_down["Peak_type"] == 'Promoter') & (Proximal_Distal_down["Annotation_method"] == 'Proximal_anno')]
        Proximal_Distal_down_remaining =  Proximal_Distal_down[~Proximal_Distal_down.index.isin(Proximal_Distal_down_promoter.index)].sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak')
        Proximal_Distal_down = pd.concat([Proximal_Distal_down_promoter, Proximal_Distal_down_remaining]).sort_index().drop_duplicates()
        Genelist = Proximal_Distal_differential.loc[:,'Gene'].unique().tolist()
        Genelist_up = Proximal_Distal_up.loc[:,'Gene'].unique().tolist()
        Genelist_down = Proximal_Distal_down.loc[:,'Gene'].unique().tolist()

      Proximal_Distal_differential.to_csv(peak_name + '_'+prefix + '_annotated.txt', index=False, sep='\t' )
      Proximal_Distal_up.to_csv(peak_name + '_' + prefix + '_annotated_up.txt', index=False, sep='\t' )
      Proximal_Distal_down.to_csv(peak_name + '_' + prefix + '_annotated_down.txt', index=False, sep='\t' )
      pd.DataFrame(Genelist).to_csv(peak_name + '_' + prefix + '_annotated_genelist.txt', index=False, header=False,sep='\t' )
      pd.DataFrame(Genelist_up).to_csv(peak_name + '_' + prefix + '_annotated_genelist_up.txt', index=False, header=False,sep='\t' )
      pd.DataFrame(Genelist_down).to_csv(peak_name + '_' + prefix + '_annotated_genelist_down.txt', index=False, header=False,sep='\t' )


# RUN FUNCTION
peak_annotation(peak_anno_anchor1=args.PEAK_ANCHOR1,peak_anno_anchor2=args.PEAK_ANCHOR2,peak_anno=args.PEAK_ANNO,bed2D_index_anno=args.BED2D, peak_name=args.PEAK_NAME, prefix=args.PREFIX, proximity_unannotated=args.PROXIMITY_UNANNOTATED, mode=args.MODE, multiple_anno=args.MULTIPLE_ANNO, promoter_start=args.PROMOTER_START, promoter_end=args.PROMOTER_END, skip_promoter_promoter= args.SKIP_PROMOTER_PROMOTER, interaction_threshold=args.INTERACTION_THRESHOLD, peak_differential=args.PEAK_DIFFERENTIAL, log2FC_column=args.LOG2FC_COLUMN, padj_column=args.PADJ_COLUMN, log2FC=args.LOG2FC, padj=args.PADJ, skip_expression=args.SKIP_EXPRESSION)
