#!/usr/bin/env python

### PROCESS 14 DIFFERENTIAL: DIFFERTIAL PEAKS ASSOCIATES WITH DIFFERENTIAL GENE EXPRESSION ###
import pandas as pd
import numpy as np
import argparse

# PARSE ARGUMENTS
Description = 'IDENTIFIES DIFFERTIAL PEAKS ASSOCIATES WITH DIFFERENTIAL GENE EXPRESSION'
Epilog = """Example usage: differential_expression_associated_peaks.py <ANNOTATED_PEAKS> <EXPRESSION> --prefix <PREFIX> --sample <SAMPLE> """

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

# Arguments
argParser.add_argument('ANNOTATED_PEAKS', help="Annotated and aggregated interactions.")
argParser.add_argument('EXPRESSION', help="Specifies path to file that contain information about differential expression between the two conditions. The first column must contain gene symbol. The column for log2FC/padj can be specified using --expression_log2FC_column and --expression_padj_column respectively, default: 3 & 9 (standard DESeq2 format). Note: make sure that you use the same direction for the comparison in --peak_differential and --expression.")
argParser.add_argument('--prefix', dest="PREFIX", help="Prefix for output file.")
argParser.add_argument('--sample', dest="SAMPLE", help="Name of sample.")
argParser.add_argument('--log2FC_column', dest="LOG2FC_COLUMN", help="Log2FC column for differential peaks.", type=int)
argParser.add_argument('--padj_column', dest="PADJ_COLUMN", help="Padj column for differential peaks.", type=int)
argParser.add_argument('--log2FC', dest="LOG2FC", help="Log2FC threshold for differential peaks.", type=float)
argParser.add_argument('--padj', dest="PADJ", help="Padj threshold for differential peaks.", type=float)
argParser.add_argument('--expression_log2FC_column', dest="EXPRESSION_LOG2FC_COLUMN", help="Log2FC column for differential expression.", type=int)
argParser.add_argument('--expression_padj_column', dest="EXPRESSION_PADJ_COLUMN", help="Padj column for differential expression.", type=int)
argParser.add_argument('--expression_log2FC', dest="EXPRESSION_LOG2FC", help="Log2FC threshold for differential expression.", type=float)
argParser.add_argument('--expression_padj', dest="EXPRESSION_PADJ", help="Padj threshold for differential expression.", type=float)

args = argParser.parse_args()

# DEFINE FUNCTION
def differential_expression_associated_peaks(annotated_peaks, expression, prefix, sample, log2FC_column, padj_column, log2FC, padj, expression_log2FC_column, expression_padj_column, expression_log2FC, expression_padj):

  #Loading data
  annotated_peaks = pd.read_table(annotated_peaks, index_col=0).sort_index()
  expression = pd.read_table(expression, index_col=0).sort_index()
  expression = expression.iloc[:, [expression_log2FC_column-2, expression_padj_column-2]]
  expression.columns = ['Gene_log2FC', 'Gene_padj']

  #Adding gene stats to annotated peaks
  annotated_peaks_expression = annotated_peaks.reset_index().merge(expression, how='left', left_on='Gene', right_on='symbol').set_index('Peak')

  #Finding differntial peaks associated with changes in expression (sepeartion based on.annotation proximal/distal and activating/prepressive)
  annotated_peaks_expression_proximal_activating = annotated_peaks_expression[(annotated_peaks_expression.Annotation_method.isin(['Proximal_anno'])) & (annotated_peaks_expression.padj < padj) & (annotated_peaks_expression.Gene_padj < expression_padj) & (((annotated_peaks_expression.log2FC > log2FC) & (annotated_peaks_expression.Gene_log2FC > expression_log2FC)) | ((annotated_peaks_expression.log2FC < -log2FC) & (annotated_peaks_expression.Gene_log2FC < -expression_log2FC)))]
  annotated_peaks_expression_proximal_repressive = annotated_peaks_expression[(annotated_peaks_expression.Annotation_method.isin(['Proximal_anno'])) & (annotated_peaks_expression.padj < padj) & (annotated_peaks_expression.Gene_padj < expression_padj) & (((annotated_peaks_expression.log2FC < -log2FC) & (annotated_peaks_expression.Gene_log2FC > expression_log2FC)) | ((annotated_peaks_expression.log2FC > log2FC) & (annotated_peaks_expression.Gene_log2FC < -expression_log2FC)))]
  annotated_peaks_expression_distal_activating = annotated_peaks_expression[(annotated_peaks_expression.Annotation_method == 'Interaction_anno') & (annotated_peaks_expression.padj < padj) & (annotated_peaks_expression.Gene_padj < expression_padj) & (((annotated_peaks_expression.log2FC >log2FC) & (annotated_peaks_expression.Gene_log2FC > expression_log2FC)) | ((annotated_peaks_expression.log2FC < -log2FC) & (annotated_peaks_expression.Gene_log2FC < -expression_log2FC)))]
  annotated_peaks_expression_distal_repressive = annotated_peaks_expression[(annotated_peaks_expression.Annotation_method =='Interaction_anno') & (annotated_peaks_expression.padj < padj) & (annotated_peaks_expression.Gene_padj < expression_padj) & (((annotated_peaks_expression.log2FC < -log2FC) & (annotated_peaks_expression.Gene_log2FC > expression_log2FC)) | ((annotated_peaks_expression.log2FC > log2FC) & (annotated_peaks_expression.Gene_log2FC < -expression_log2FC)))]

  #Adding activatin/repressive funktion to all annotated peaks
  annotated_peaks_expression['Activating_or_repressive'] = np.where((annotated_peaks_expression.index.isin(annotated_peaks_expression_proximal_activating.index)) | (annotated_peaks_expression.index.isin(annotated_peaks_expression_distal_activating.index)), 'Activating', np.where((annotated_peaks_expression.index.isin(annotated_peaks_expression_proximal_repressive.index)) | (annotated_peaks_expression.index.isin(annotated_peaks_expression_distal_repressive.index)), 'Repressive', ''))

  #Save files
  annotated_peaks_expression.to_csv(sample + '_' + prefix + '_annotated_differential_expression.txt', index=True, sep='\t' )
  annotated_peaks_expression_proximal_activating.to_csv(sample + '_' + prefix + '_annotated_differential_expression_proximal_activating.txt', index=True, sep='\t' )
  annotated_peaks_expression_proximal_repressive.to_csv(sample + '_' + prefix + '_annotated_differential_expression_proximal_repressive.txt', index=True, sep='\t' )
  annotated_peaks_expression_distal_activating.to_csv(sample + '_' + prefix + '_annotated_differential_expression_distal_activating.txt', index=True, sep='\t' )
  annotated_peaks_expression_distal_repressive.to_csv(sample + '_' + prefix + '_annotated_differential_expression_distal_repressive.txt', index=True, sep='\t' )

# RUN FUNCTION
differential_expression_associated_peaks(annotated_peaks=args.ANNOTATED_PEAKS, expression=args.EXPRESSION, prefix=args.PREFIX, sample=args.SAMPLE, log2FC_column=args.LOG2FC_COLUMN, padj_column=args.PADJ_COLUMN, log2FC=args.LOG2FC, padj=args.PADJ, expression_log2FC_column=args.EXPRESSION_LOG2FC_COLUMN, expression_padj_column=args.EXPRESSION_PADJ_COLUMN, expression_log2FC=args.EXPRESSION_LOG2FC, expression_padj=args.EXPRESSION_PADJ)
