#!/usr/bin/env python

### PROCESS 3: JOIN_ANNOTATED_INTERACTIONS ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'JOIN INTERACTION ANCHOR POINTS AFTER HOMER ANNOTATION'
Epilog = """Example usage: join_annotated_interactions.py <ANCHOR1_ANNO> <ANCHOR2_ANNO> <BED2D> --prefix <PREFIX> --binsize <BINSIZE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

argParser.add_argument('ANCHOR1_ANNO', help="Annotated anchor1 regions.")
argParser.add_argument('ANCHOR2_ANNO', help="Annotated anchor2 regions.")
argParser.add_argument('BED2D', help="2D-bed interactions.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--binsize', dest="BINSIZE", help="Interaction bin size.", type=int)

args = argParser.parse_args()

# DEFINE FUNCTION
def join_annotated_interactions(anchor1_anno,anchor2_anno,bed2D_index, prefix, binsize):
    anchor1_anno = pd.read_table(anchor1_anno, index_col=0).sort_index()
    anchor2_anno = pd.read_table(anchor2_anno, index_col=0).sort_index()
    bed2D = pd.read_table(bed2D_index,index_col=0).sort_index()

    anchor_merge = pd.concat([bed2D[bed2D.columns[0:9]], anchor1_anno[anchor1_anno.columns[6: len(anchor1_anno.columns)]], anchor2_anno[anchor2_anno.columns[6: len(anchor2_anno.columns)]]], axis=1)
    anchor_merge.columns = ['chr1', 's1', 'e1', 'chr2', 's2', 'e2', 'cc', 'P-Value_Bias',
       'Q-Value_Bias', 'Annotation_1', 'Detailed_Annotation_1', 'Distance_to_TSS_1',
       'Nearest_PromoterID_1', 'Entrez_ID_1', 'Nearest_Unigene_1', 'Nearest_Refseq_1',
       'Nearest_Ensembl_1', 'Gene_Name_1', 'Gene_Alias_1', 'Gene_Description_1',
       'Gene_Type_1', 'Annotation_2', 'Detailed_Annotation_2', 'Distance_to_TSS_2',
       'Nearest_PromoterID_2', 'Entrez_ID_2', 'Nearest_Unigene_2', 'Nearest_Refseq_2',
       'Nearest_Ensembl_2', 'Gene_Name_2', 'Gene_Alias_2', 'Gene_Description_2',
       'Gene_Type_2']

    anchor_merge['TSS_1'] = np.where(abs(anchor_merge['Distance_to_TSS_1']) <= binsize/2, 1, 0) #Defining anchor point 1 as promoter bins (TSS_1=1) if the bin contain a TSS (distance to closest TSS < half the binsize)
    anchor_merge['TSS_2'] = np.where(abs(anchor_merge['Distance_to_TSS_2']) <= binsize/2, 1, 0) #Defining anchor point 2 as promoter bins (TSS_2=1) if the bin contain a TSS (distance to closest TSS < half the binsize)
    anchor_merge.to_csv(prefix + '_HOMER_annotated_interactions.txt', index=True, sep='\t' )

# RUN FUNCTION
join_annotated_interactions(anchor1_anno=args.ANCHOR1_ANNO,anchor2_anno=args.ANCHOR2_ANNO,bed2D_index=args.BED2D, prefix=args.PREFIX, binsize=args.BINSIZE)
