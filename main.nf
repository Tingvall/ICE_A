#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  Usage:
  The typical command for running the pipeline is as follows:
    nextflow run PLAC_anno_2.nf --bed2D interactions.bed  --genome mm10 --peaks peaks.txt

  Required inputs:
    --peaks [file]                  Patj to text file specifying the name of the peak files(s) in the first column and the path to the file(s) in the second column (For example see: [peaks.txt](example_files/peaks.txt)). The recomended input format for the peak files are 6 column bed files (more columns are allowed but will be ignored): chr, start, end, peakID, peak score, strand. It is also possible to use a 3 column bed but then peakIDs are automaically generated and peak score and strand information are not available.
    --bed2D [file]                  Path 2D-BED file containing genomic interactions. Currently the pipeline is designed for 2D-bed files created by FitHIC (https://github.com/ay-lab/fithic). This input can be replcaed by bed2D_anno if an alredy annotated 2D-bed file is available and the arguement skip_anno is used.
    --genome                        Specification of genome for annotation (e.g. mm10).

  Optional input:
    --genes [file]                  Path to textfile with gene names (specify???), that is used for filtering of interactions associated with the specified genes. nly used when the option --filtering_genes is specified or if --network_mode is set to genes.
    --bed2D_anno [file]             Specifies path to annotated 2D-bed file if --skip_anno is used.

  Arguments - General:
    --mode                          Define which mode to run the pipeline in. The options are basic (default), multiple or differntial.
    --outdir                        Specifies the output driectory (default: ./results).
    --prefix                        Prefix used for interactions (default: PLACseq).
    --proximity_unannotated         Specifies if unannotated distal peaks should be annotated by proximity annotation (default: false).
    --multiple_anno                 Defines how to handle peaks annotated to more than one promoter. Options are keep (all anotations are kept with one row for each annotation), concetrate (the annotated peak file is concetrated to only incude one row per peak but information about all annotations are kept) and qvalue (only one annotation per peak is kept. The annotation is decided by the interaction with the lowset qvalue). Default is: concentrate.
    --skip_anno                     Skip the HOMER annotation of the interactions (requires specification of path to annotated 2D-bed by using the argumnet --bed2D_anno).
    --annotate_interactions         Specifes if interaction-centered annotation with peak overlap should be performed. Only valid if --complete is set to false.
    --network                       Specifes if files for network visualization in Cytoskape should be created. Only valid if --complete is set to false.
    --network_mode                  Defines mode network. Options are all (all interaction in the 2D-bed file), factor (all interaction with at least on peak overlap either anchor point) or genes (interactions associates with a genelist, provided by --genes).
    --complete                      If specified, all available processes for the selected mode and provided inputs are run.
    --save_tmp                      If used, all intermediate files are saved in the directory ./tmp. Can be useful for investigating promblems. Default: false.
    --help                          Help message is shown.

  Arguments - Multiple mode specific:
    --upset_plot                    If specified, Upset plot of peak overlap will be created. Only valid if --complete is set to false.
    --circos_plot                   If specified, Circos plot of peak overlap will be created. Only valid if --complete is set to false.
    --filter_genes                  Specifies if additional plot (Upset and/or Circos plots) should be created based on interactions filtered by provided genelist (default: false). This option requires that a genelist is provided with the argument --genes.

    Arguments - Multiple mode specific:
    --peak_differntial [file]       Path to textfile that contain log2FC and adjusted p-value from differntial analyis. The 1st column should contain peakID matching the peakID in the 4th column of the input bed file. Standard DESeq2 output is expected (with log2FC in the 3rd column and padj in the 9th column), but other formats are accepted as well is the column corresponding to log2FC and padj are specified with the aruguments --log2FC_column and --padj column.
    --log2FC_column                 Specifies which column in the differntial peak file that contain log2FC values (default:3, standard DESeq2 output)
    --padj_column                   Specifies which column in the differntial peak file that contain adjusted p-values (default:9, standard DESeq2 output)
    --log2FC                        Log2FC treshold for differtial peak analysis
    --padj                          Significance level for differtial peak analysis
    --expression [file]

  """.stripIndent()
  }


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * VALIDATE INPUTS
 */
if (!params.skip_anno) {
  if (params.bed2D)     { ch_bed2D = Channel.fromPath(params.bed2D, checkIfExists: true) } else { exit 1, '2D-bed file not found' }
}
else{
  ch_bed2D = Channel.empty()
}

if (params.peaks)     { ch_peaks = Channel.fromPath(params.peaks, checkIfExists: true) } else { exit 1, 'Peaks not specified' }
    ch_peaks
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.peak_name, [ file(row.peak_file) ] ] }
        .set { ch_peaks_split}

if (!params.genome)      { exit 1, 'Refence genome not specified' }


if (params.network && params.network_mode == 'genes') {
  if (params.genes)     { ch_genes = Channel.fromPath(params.genes, checkIfExists: true) } else { exit 1, 'Genes not specified' }
}
else {
  ch_genes = file(params.genes)
}

if (params.mode == 'differential') {
  if (params.peak_differential)     { ch_peak_differential = Channel.fromPath(params.peak_differential, checkIfExists: true) } else { exit 1, 'Peak file log2FC and adjusted p-value not specified' }
}
else {
  ch_peak_differential = file(params.genes)
}
ch_peak_differential.into { ch_peak_differential_1; ch_peak_differential_2 }

if (params.mode == 'basic') {
println ("""
        ===========================================================================================
                                            PLAC-SEQ ANNOTATION
                                           ---------------------
                                           Basic Annotation Mode
        ===========================================================================================
        General parameters:
        Prefix: ${params.prefix}
        Skip 2D-bed annotation: ${params.skip_anno}
        2D-bed: ${params.bed2D}
        Reference genome: ${params.genome}
        Outdir: ${params.outdir}
        Peak file:  ${params.peaks}
        Proximity annotate unannotated distal regions: ${params.proximity_unannotated}
        Mode for multiple annotation of Peaks: ${params.multiple_anno}

        Basic Annotation mode:
        Perform interaction-based annotation: ${params.annotate_interactions}
        Create network of peak annotations: ${params.network}
        Mode for network: ${params.network_mode}
        Genelist for netowrk: ${params.genes}
        Run in complete mode (all available processes): ${params.complete}
        ===========================================================================================
        """)
}

else if (params.mode == 'multiple') {
println ("""
        ===========================================================================================
                                            PLAC-SEQ ANNOTATION
                                          ------------------------
                                          Multiple Annotation Mode
        ===========================================================================================
        Prefix: ${params.prefix}
        Skip 2D-bed annotation: ${params.skip_anno}
        2D-bed: ${params.bed2D}
        Reference genome: ${params.genome}
        Outdir: ${params.outdir}
        Peak file:  ${params.peaks}
        Proximity annotate unannotated distal regions: ${params.proximity_unannotated}
        Mode for multiple annotation of Peaks: ${params.multiple_anno}

        Multiple Annotation mode:
        Perform interaction-based annotation: ${params.annotate_interactions}
        Create network of peak annotations: ${params.network}
        Mode for network: ${params.network_mode}
        Genelist for netowrk: ${params.genes}
        Upset plot: ${params.upset_plot}
        Circos plot: ${params.circos_plot}
        Run in complete mode (all available processes): ${params.complete}
        ===========================================================================================
        """)
}

else if (params.mode == 'differntial') {
println ("""
        ===========================================================================================
                                            PLAC-SEQ ANNOTATION
                                        ---------------------------
                                        Differntial Annotation Mode
        ===========================================================================================
        Prefix: ${params.prefix}
        Skip 2D-bed annotation: ${params.skip_anno}
        2D-bed: ${params.bed2D}
        Reference genome: ${params.genome}
        Outdir: ${params.outdir}
        Peak file:  ${params.peaks}
        Proximity annotate unannotated distal regions: ${params.proximity_unannotated}
        Mode for multiple annotation of Peaks: ${params.multiple_anno}

        Multiple Annotation mode:
        ===========================================================================================
        """)
}


/*
 * 1. 2D-BED SPLIT: SPLIT 2D-BED FILE INTO 2 BED FILES FOR ANNOTATION
 */
process BED2D_SPLIT {
    publishDir "${params.outdir}/tmp/process1", mode: 'copy', enabled: params.save_tmp


    when:
    !params.skip_anno

    input:
    path bed2D from ch_bed2D
    val prefix from Channel.value(params.prefix)

    output:
    path "${prefix}_anchor1.bed" into ch_anchor1
    path "${prefix}_anchor2.bed" into ch_anchor2
    path "${prefix}_index.bed" into ch_bed2D_index      // 2D-bed file with index column first

    script:
    """
    awk -v OFS='\t' 'FNR==1{a="index"} FNR>1{a=NR-1} {print a,\$0}' $bed2D > ${prefix}_index.bed
    awk -v OFS='\t' '{if (NR!=1) {print \$2,\$3,\$4,\$1 }}' ${prefix}_index.bed >  ${prefix}_anchor1.bed
    awk -v OFS='\t' '{if (NR!=1) {print \$5,\$6,\$7,\$1}}' ${prefix}_index.bed >  ${prefix}_anchor2.bed
    """
}


/*
 * 2. HOMER ANNOTATION PLAC-seq: ANNOTATION OF EACH ANCHOR REGION USING HOMER
 */
process ANNOTATE_INTERACTION {
    publishDir "${params.outdir}/tmp/process2", mode: 'copy', enabled: params.save_tmp

    when:
    !params.skip_anno

    input:
    path anchor1 from ch_anchor1
    path anchor2 from ch_anchor2
    val genome from Channel.value(params.genome)

    output:
    path "${anchor1.baseName}_anno.txt" into ch_anchor1_anno       // Annotated anchor bed files
    path "${anchor2.baseName}_anno.txt" into ch_anchor2_anno

    script:
    """
    annotatePeaks.pl $anchor1 $genome > ${anchor1.baseName}_anno.txt
    annotatePeaks.pl $anchor2 $genome > ${anchor2.baseName}_anno.txt
    """
}

/*
 * 3. MERGE ANNOTATED ANCHOR REGIONS
 */
process JOIN_ANNOTATED_INTERACTIONS {
    publishDir "${params.outdir}/tmp/process3", mode: 'copy', enabled: params.save_tmp
    publishDir "${params.outdir}/Interaction_Annotation/", mode: 'copy', enabled: params.annotate_interactions | params.complete

    when:
    !params.skip_anno

    input:
    path anchor1_anno from ch_anchor1_anno
    path anchor2_anno from ch_anchor2_anno
    path bed2D_index from ch_bed2D_index
    val prefix from Channel.value(params.prefix)

    output:
    path "${prefix}_HOMER_annotated_interactions.txt" into ch_bed2D_anno

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np

    anchor1_anno = pd.read_table("${anchor1_anno}", index_col=0).sort_index()
    anchor2_anno = pd.read_table("${anchor2_anno}", index_col=0).sort_index()
    bed2D = pd.read_table("${bed2D_index}",index_col=0).sort_index()

    anchor_merge = pd.concat([bed2D[bed2D.columns[0:9]], anchor1_anno[anchor1_anno.columns[6: len(anchor1_anno.columns)]], anchor2_anno[anchor2_anno.columns[6: len(anchor2_anno.columns)]]], axis=1)
    anchor_merge.columns = ['chr1', 's1', 'e1', 'chr2', 's2', 'e2', 'cc', 'P-Value_Bias',
       'Q-Value_Bias', 'Annotation_1', 'Detailed_Annotation_1', 'Distance_to_TSS_1',
       'Nearest_PromoterID_1', 'Entrez_ID_1', 'Nearest_Unigene_1', 'Nearest_Refseq_1',
       'Nearest_Ensembl_1', 'Gene_Name_1', 'Gene_Alias_1', 'Gene_Description_1',
       'Gene_Type_1', 'Annotation_2', 'Detailed_Annotation_2', 'Distance_to_TSS_2',
       'Nearest_PromoterID_2', 'Entrez_ID_2', 'Nearest_Unigene_2', 'Nearest_Refseq_2',
       'Nearest_Ensembl_2', 'Gene_Name_2', 'Gene_Alias_2', 'Gene_Description_2',
       'Gene_Type_2']

    anchor_merge['TSS_1'] = np.where(abs(anchor_merge['Distance_to_TSS_1']) <= 2500, 1, 0)
    anchor_merge['TSS_2'] = np.where(abs(anchor_merge['Distance_to_TSS_2']) <= 2500, 1, 0)
    anchor_merge.to_csv("${prefix}_HOMER_annotated_interactions.txt", index=True, sep='\t' )
    """
}

if (params.skip_anno) {
  if (params.bed2D_anno)     { ch_bed2D_anno = Channel.fromPath(params.bed2D_anno, checkIfExists: true) } else { exit 1, 'Annotated 2D-bed file not found' }
}

  /*
   * 4. HOMER ANNOTATION PEAKS: ANNOTATION OF PEAK files USING HOMER
   */
  process ANNOTATE_PEAKS {
      publishDir "${params.outdir}/tmp/process4", mode: 'copy', enabled: params.save_tmp

      input:
      tuple val(peak_name), path(peak_file) from ch_peaks_split
      val genome from Channel.value(params.genome)


      output:
      tuple val(peak_name), file("${peak_name}_anno.txt") into ch_peak_anno
      tuple val(peak_name), file("${peak_name}.bed") into ch_peak_bed_1, ch_peak_bed_2,ch_peak_bed_3,ch_peak_bed_4

      script:
      """
      annotatePeaks.pl $peak_file $genome > ${peak_name}_anno.txt
      awk -v OFS='\t' '{if (NR!=1) {print \$2,\$3,\$4,\$1,\$6 }}' ${peak_name}_anno.txt >  ${peak_name}.bed
      """
  }

  /*
   * 5. SPLIT ANNOTATED 2D-BED: ANNOTATED 2D-BED SPLIT FOR PEAK OVERLAP
   */
  process SPLIT_ANNOTATED_INTERACTIONS {
    publishDir "${params.outdir}/tmp/process5", mode: 'copy', enabled: params.save_tmp

    input:
    path bed2D_anno from ch_bed2D_anno
    val prefix from Channel.value(params.prefix)

    output:
    path "${prefix}_anchor1_anno.bed" into ch_bed2D_anno_split_anchor1_1, ch_bed2D_anno_split_anchor1_2
    path "${prefix}_anchor2_anno.bed" into ch_bed2D_anno_split_anchor2_1, ch_bed2D_anno_split_anchor2_2
    path "${prefix}_index_anno.bed" into ch_bed2D_index_anno_1, ch_bed2D_index_anno_2

    script:
    """
    awk -v OFS='\t' 'FNR==1{a="index"} FNR>1{a=NR-1} {print a,\$0}' $bed2D_anno > ${prefix}_index_anno.bed
    awk -v OFS='\t' '{if (NR!=1) {print \$3,\$4,\$5,\$2 }}' ${prefix}_index_anno.bed >  ${prefix}_anchor1_anno.bed
    awk -v OFS='\t' '{if (NR!=1) {print \$6,\$7,\$8,\$2}}' ${prefix}_index_anno.bed >  ${prefix}_anchor2_anno.bed
    """
  }

  /*
   * 6. BEDTOOLS INTERSECT PEAK CENTERED: OVERLAPPING PEAKS WITH 2D-BED ANCHOR POINTS
   */
  process PEAK_INTERACTION_INTERSECT {
    publishDir "${params.outdir}/tmp/process6", mode: 'copy', enabled: params.save_tmp

    input:
    set val(peak_name), file(peak_bed), file(bed2D_anno_split_anchor1), file(bed2D_anno_split_anchor2) from ch_peak_bed_1.combine(ch_bed2D_anno_split_anchor1_1).combine(ch_bed2D_anno_split_anchor2_1).groupTuple()


    output:
    tuple val(peak_name), path("${peak_name}_anchor_1.bed") into ch_peak_anno_anchor1
    tuple val(peak_name), path("${peak_name}_anchor_2.bed") into ch_peak_anno_anchor2

    script:
    """
    bedtools intersect -wa -wb -a $peak_bed -b $bed2D_anno_split_anchor1 > ${peak_name}_anchor_1.bed
    bedtools intersect -wa -wb -a $peak_bed -b $bed2D_anno_split_anchor2 > ${peak_name}_anchor_2.bed

    """
  }

  /*
   * 7. PEAK ANNOTATION: PEAKS ANNOTATED BY PROXIMITY OR PLAC-SEQ BASED ANNOTATION
   */
  process PEAK_INTERACTION_BASED_ANNOTATION {
    publishDir "${params.outdir}/tmp/process7", mode: 'copy', enabled: params.save_tmp
    publishDir "${params.outdir}/Peak_Annotation/${peak_name}", mode: 'copy'

    input:
    set val(peak_name), file(peak_anno_anchor1), file(peak_anno_anchor2), file(peak_anno), file(bed2D_index_anno) from ch_peak_anno_anchor1.join(ch_peak_anno_anchor2).join(ch_peak_anno).combine(ch_bed2D_index_anno_1)
    val proximity_unannotated from Channel.value(params.proximity_unannotated)
    val multiple_anno from Channel.value(params.multiple_anno)
    val prefix from Channel.value(params.prefix)
    val mode from Channel.value(params.mode)

    //Differntial mode specific
    path peak_differential from ch_peak_differential_1
    val log2FC_column from Channel.value(params.log2FC_column)
    val padj_column from Channel.value(params.padj_column)
    val log2FC from Channel.value(params.log2FC)
    val padj from Channel.value(params.padj)

    output:
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated.txt") into ch_peak_PLACseq_annotated
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_up.txt") optional true into ch_peak_PLACseq_annotated_up
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_down.txt") optional true into ch_peak_PLACseq_annotated_down

    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_genelist.txt") into ch_peak_PLACseq_annotated_genelist
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_genelist_up.txt") optional true into ch_peak_PLACseq_annotated_genelist_up
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_genelist_down.txt") optional true into ch_peak_PLACseq_annotated_genelist_down


    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np

    # Column names for loaded data
    peak_anchor1_name = ('peak_chr', 'peak_start', 'peak_end','Peak_score', 'anchor1_chr', 'anchor1_start', 'anchor1_end', 'anchor1_id')
    peak_anchor2_name = ('peak_chr', 'peak_start', 'peak_end',  'Peak_score', 'anchor2_chr', 'anchor2_start', 'anchor2_end', 'anchor2_id')

    # Load peak overlap for anchor 1 & 2, as well as annotated peak & 2D-bed files
    peak_anchor1 = pd.read_table("${peak_anno_anchor1}", index_col=3, names=peak_anchor1_name).sort_index()
    peak_anchor2 = pd.read_table("${peak_anno_anchor2}", index_col=3, names=peak_anchor2_name).sort_index()
    peak_anno = pd.read_table("${peak_anno}",index_col=0).sort_index()
    bed2D_anno = pd.read_table("${bed2D_index_anno}", index_col=1).sort_index().iloc[:,1:]

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
    Peak_overlap_merge['Annotation'] = np.where(abs(Peak_overlap_merge['Distance_to_TSS']) <= 2500, 'Promoter', (np.where(abs(Peak_overlap_merge['Distance_to_TSS']) <= 10000, 'Proximal_anno', 'Plac_anno')))

    # Extrating promoter and proximity annotated peak, adding Q_value column (for filtering) and renaming columns
    Proximal = Peak_overlap_merge.loc[Peak_overlap_merge['Annotation'].isin(['Promoter','Proximal_anno']),['Chr', 'Start', 'End', 'Peak_score', 'EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal', 'Annotation']].drop_duplicates()
    Proximal['Q-value'] = np.nan
    Proximal.columns = ['Chr', 'Start', 'End', 'Peak_score', 'EntrezID', 'Refseq','Ensembl', 'Gene', 'Annotation', 'Q-value']

    # Extracting PLAC-seq annotated peaks
    Distal = Peak_overlap_merge.loc[Peak_overlap_merge['Annotation'].isin(['Plac_anno']),:].dropna(subset=['InteractionID'])
    Distal = Distal.loc[(Distal['Anchor_Overlap_TSS'] == 0) & (Distal['Anchor_Interaction_TSS'] == 1),['Chr', 'Start', 'End', 'Peak_score', 'EntrezID_Interaction', 'Refseq_Interaction','Ensembl_Interaction', 'Gene_Interaction', 'Annotation', 'Q-Value']].drop_duplicates()
    Distal.columns = ['Chr', 'Start', 'End', 'Peak_score', 'EntrezID', 'Refseq','Ensembl', 'Gene', 'Annotation', 'Q-value']

    # Merge proximity and PLAC-seq annotated peaks
    Proximal_Distal = pd.concat([Proximal, Distal]).sort_index().rename_axis('Peak')

    # Extracting unannotated distal peaks (not overlapping 2D-bed)
    Unannotated = Peak_overlap_merge.loc[:,['Chr', 'Start', 'End', 'Peak_score', 'EntrezID_Proximal', 'Refseq_Proximal','Ensembl_Proximal', 'Gene_Proximal']][~Peak_overlap_merge.index.isin(Proximal_Distal.index)].drop_duplicates()
    Unannotated['Annotation'], Unannotated['Q-value'] = ['Distal_no_Interaction', np.NaN]
    Unannotated.columns=['Chr', 'Start', 'End', 'Peak_score', 'EntrezID', 'Refseq','Ensembl', 'Gene', 'Annotation', 'Q-value']

    # Annotate unannotated distal peaks by proximity annotation
    proximity_unannotated = "${proximity_unannotated}"
    if proximity_unannotated == 'true':
        Proximal_Distal = pd.concat([Proximal_Distal, Unannotated]).sort_index().rename_axis('Peak')


    Proximal_Distal['Start'] = Proximal_Distal['Start']-1

    # Organizing and saving annotated files/genelsits
    mode = "${mode}"

    # Basic/Multiple mode: Handling of peaks annotating to several genes
    if (mode=='basic' or mode=='multiple'):
      multiple_anno = "${multiple_anno}"
      if multiple_anno == 'keep':
        Proximal_Distal = Proximal_Distal
        Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()
      elif multiple_anno == 'qvalue':
        Proximal_Distal = Proximal_Distal.sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak').sort_index()
        Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()
      elif multiple_anno == 'concentrate':
        Genelist = Proximal_Distal.loc[:,'Gene'].unique().tolist()
        Proximal_Distal = Proximal_Distal.groupby('Peak').agg(lambda x: ', '.join(list(x.unique().astype(str))))

      Proximal_Distal.to_csv("${peak_name}_${prefix}_annotated.txt", index=False, sep='\t' )
      pd.DataFrame(Genelist).to_csv("${peak_name}_${prefix}_annotated_genelist.txt", index=False, header=False,sep='\t' )

    # Differntial mode: Creting peak annotation/genelists for UP/DOWN peak_start
    if mode == 'differential':
      peak_differntial = pd.read_table("${peak_differential}", index_col=0).sort_index()
      peak_differntial = peak_differntial.iloc[:, [${log2FC_column}-2, ${padj_column}-2]]
      peak_differntial.columns = ['log2FC', 'padj']
      Proximal_Distal_differential=Proximal_Distal.merge(peak_differntial, left_index=True, right_index=True, how = 'left')
      Proximal_Distal_differential.index.name = 'Peak'
      Proximal_Distal_up=Proximal_Distal_differential[(Proximal_Distal_differential.padj <= ${padj}) & (Proximal_Distal_differential.log2FC >=${log2FC})]
      Proximal_Distal_down=Proximal_Distal_differential[(Proximal_Distal_differential.padj <= ${padj}) & (Proximal_Distal_differential.log2FC <= -${log2FC})]

      # Differnetial mode: Handling of peaks annotating to several genes
      multiple_anno = "${multiple_anno}"
      if multiple_anno == 'keep':
        Genelist = Proximal_Distal_differential.loc[:,'Gene'].unique().tolist()
        Genelist_up = Proximal_Distal_up.loc[:,'Gene'].unique().tolist()
        Genelist_down = Proximal_Distal_down.loc[:,'Gene'].unique().tolist()
      elif multiple_anno == 'qvalue':
        Proximal_Distal_differential = Proximal_Distal_differential.sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak').sort_index()
        Proximal_Distal_up = Proximal_Distal_up.sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak').sort_index()
        Proximal_Distal_down = Proximal_Distal_down.sort_values('Q-value').reset_index().drop_duplicates(subset=['Peak'],keep='first').set_index('Peak').sort_index()
        Genelist = Proximal_Distal_differential.loc[:,'Gene'].unique().tolist()
        Genelist_up = Proximal_Distal_up.loc[:,'Gene'].unique().tolist()
        Genelist_down = Proximal_Distal_down.loc[:,'Gene'].unique().tolist()
      elif multiple_anno == 'concentrate':
        Genelist = Proximal_Distal_differential.loc[:,'Gene'].unique().tolist()
        Genelist_up = Proximal_Distal_up.loc[:,'Gene'].unique().tolist()
        Genelist_down = Proximal_Distal_down.loc[:,'Gene'].unique().tolist()
        Proximal_Distal_differential = Proximal_Distal_differential.groupby('Peak').agg(lambda x: ', '.join(list(x.unique().astype(str))))
        Proximal_Distal_up = Proximal_Distal_up.groupby('Peak').agg(lambda x: ', '.join(list(x.unique().astype(str))))
        Proximal_Distal_down = Proximal_Distal_down.groupby('Peak').agg(lambda x: ', '.join(list(x.unique().astype(str))))

      Proximal_Distal_differential.to_csv("${peak_name}_${prefix}_annotated.txt", index=False, sep='\t' )
      Proximal_Distal_up.to_csv("${peak_name}_${prefix}_annotated_up.txt", index=False, sep='\t' )
      Proximal_Distal_down.to_csv("${peak_name}_${prefix}_annotated_down.txt", index=False, sep='\t' )
      pd.DataFrame(Genelist).to_csv("${peak_name}_${prefix}_annotated_genelist.txt", index=False, header=False,sep='\t' )
      pd.DataFrame(Genelist_up).to_csv("${peak_name}_${prefix}_annotated_genelist_up.txt", index=False, header=False,sep='\t' )
      pd.DataFrame(Genelist_down).to_csv("${peak_name}_${prefix}_annotated_genelist_down.txt", index=False, header=False,sep='\t' )
    """
}

def criteria = multiMapCriteria {
                     peak_names: it[0]
                     peaks_beds: it[1]
                   }

ch_peak_bed_2.multiMap(criteria).set {ch_t_1}
ch_peak_bed_3.multiMap(criteria).set {ch_t_2}
ch_peak_bed_4.multiMap(criteria).set {ch_t_3}


/*
 * 8. BEDTOOLS INTERSECT INTERACTION CENTERED: OVERLAPPING PEAKS WITH 2D-BED ANCHOR POINTS
 */
process INTERACTION_PEAK_INTERSECT {
  publishDir "${params.outdir}/tmp/process8", mode: 'copy', enabled: params.save_tmp

  when:
  params.annotate_interactions | params.network | params.upset_plot | params.circos_plot | params.complete

  input:
  val peak_names from ch_t_1.peak_names.collect().map{ it2 -> it2.join(' ')}
  val peak_beds from ch_t_1.peaks_beds.collect().map{ it2 -> it2.join(' ')}
  path bed2D_anno_split_anchor1 from ch_bed2D_anno_split_anchor1_2
  path bed2D_anno_split_anchor2 from ch_bed2D_anno_split_anchor2_2


  output:
  path "Anchor_1_peak_collect.bed" into ch_anchor_1_peak_collect
  path "Anchor_2_peak_collect.bed" into ch_anchor_2_peak_collect


  script:
  if (params.mode == 'basic' | params.mode == 'differential')
    """
    bedtools intersect -wa -wb -a $bed2D_anno_split_anchor1 -b $peak_beds > Anchor_1_peak_collect.bed
    bedtools intersect -wa -wb -a $bed2D_anno_split_anchor2 -b $peak_beds > Anchor_2_peak_collect.bed
    """

  else if (params.mode == 'multiple')
    """
    bedtools intersect -wa -wb -a $bed2D_anno_split_anchor1 -b $peak_beds -names $peak_names > Anchor_1_peak_collect.bed
    bedtools intersect -wa -wb -a $bed2D_anno_split_anchor2 -b $peak_beds -names $peak_names > Anchor_2_peak_collect.bed
    """
}


/*
 * 9. INTERACTION CENTERED ANNOTATION - WITH PEAK OVERLAP
 */
process ANNOTATE_INTERACTION_WITH_PEAKS {
  publishDir "${params.outdir}/tmp/process9", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Interaction_Annotation", mode: 'copy', pattern: '*.txt', enabled: params.annotate_interactions | params.complete

  when:
  params.annotate_interactions | params.network | params.upset_plot | params.circos_plot | params.complete

  input:
  path anchor_1_peak_collect from ch_anchor_1_peak_collect
  path anchor_2_peak_collect from ch_anchor_2_peak_collect
  path bed2D_index_anno from ch_bed2D_index_anno_2
  val prefix from Channel.value(params.prefix)
  val peak_names from ch_t_2.peak_names.collect().map{ it2 -> it2.join(' ')}

  //Differntial mode specific
  path peak_differential from ch_peak_differential_2
  val log2FC_column from Channel.value(params.log2FC_column)
  val padj_column from Channel.value(params.padj_column)
  val log2FC from Channel.value(params.log2FC)
  val padj from Channel.value(params.padj)

  output:
  path "*_${prefix}_interactions.txt" into ch_interactions_by_factor
  path "${prefix}_HOMER_annotated_interactions_with_peak_overlap.txt" into ch_interactions_all

  script:
  if (params.mode == 'basic')
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np

    # Column names for loaded data
    anchor1_peak_name = ('Anchor1_Chr', 'Anchor1_Start', 'Anchor1_End', 'Peak1_Chr', 'Peak1_Start', 'Peak1_End', 'Peak1_ID', 'Peak1_score')
    anchor2_peak_name = ('Anchor2_Chr', 'Anchor2_Start', 'Anchor2_End', 'Peak2_Chr', 'Peak2_Start', 'Peak2_End', 'Peak2_ID', 'Peak2_score')

    # Load interaction centered peak overlaps 2and annotated 2D-bed
    anchor1_peaks = pd.read_table("${anchor_1_peak_collect}", index_col=3,names=anchor1_peak_name).sort_index()
    anchor2_peaks = pd.read_table("${anchor_2_peak_collect}", index_col=3,names=anchor2_peak_name).sort_index()
    bed2D_anno = pd.read_table("${bed2D_index_anno}", index_col=1).sort_index().iloc[:,1:]

    # Create Peak columns (chr:start-end) for anchor 1 & 2
    anchor1_peaks["Peak1_ID"] = anchor1_peaks["Peak1_Chr"].map(str) +':'+ (anchor1_peaks["Peak1_Start"]-1).map(str) +'-'+ anchor1_peaks["Peak1_End"].map(str)
    anchor2_peaks["Peak2_ID"] = anchor2_peaks["Peak2_Chr"].map(str) +':'+ (anchor2_peaks["Peak2_Start"]-1).map(str) +'-'+ anchor2_peaks["Peak2_End"].map(str)

    # Merging anchor points
    anchor1_peaks_anno =bed2D_anno.loc[:,['chr1', 's1','e1' ,'Entrez_ID_1', 'Gene_Name_1']].merge(anchor1_peaks.loc[:,['Peak1_ID', 'Peak1_score']], left_index=True, right_index=True, how = 'left')
    anchor2_peaks_anno =bed2D_anno.loc[:,['chr2', 's2','e2' ,'Entrez_ID_2', 'Gene_Name_2','Q-Value_Bias','TSS_1', 'TSS_2']].merge(anchor2_peaks.loc[:,['Peak2_ID', 'Peak2_score']], left_index=True, right_index=True, how = 'left')
    anchors_peaks_anno = anchor1_peaks_anno.merge(anchor2_peaks_anno, left_index=True, right_index=True,how = 'outer').drop_duplicates()

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
    anchors_peaks_anno_factor = anchors_peaks_anno[(anchors_peaks_anno['Overlap_1'] == 1) | (anchors_peaks_anno['Overlap_2'] == 1)]
    anchors_peaks_anno_factor = anchors_peaks_anno_factor.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    anchors_peaks_anno = anchors_peaks_anno.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    # Saving annotated interactions files (all interaction and interactions with peak overlap)
    anchors_peaks_anno_factor.to_csv("${peak_names}_${prefix}_interactions.txt", index=False, sep='\t' )
    anchors_peaks_anno.to_csv("${prefix}_HOMER_annotated_interactions_with_peak_overlap.txt", index=True, sep='\t' )
    """

  else if (params.mode == 'multiple')
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np

    # Column names for loaded data
    anchor1_peak_name = ('Anchor1_Chr', 'Anchor1_Start', 'Anchor1_End', 'Peak1', 'Peak1_Chr', 'Peak1_Start', 'Peak1_End', 'Peak1_ID', 'Peak1_score')
    anchor2_peak_name = ('Anchor2_Chr', 'Anchor2_Start', 'Anchor2_End', 'Peak2', 'Peak2_Chr', 'Peak2_Start', 'Peak2_End', 'Peak2_ID', 'Peak2_score')

    # Load interaction centered peak overlaps 2and annotated 2D-bed
    anchor1_peaks = pd.read_table("${anchor_1_peak_collect}", index_col=3,names=anchor1_peak_name).sort_index()
    anchor2_peaks = pd.read_table("${anchor_2_peak_collect}", index_col=3,names=anchor2_peak_name).sort_index()
    bed2D_anno = pd.read_table("${bed2D_index_anno}", index_col=1).sort_index().iloc[:,1:]

    # Create Peak columns (chr:start-end) for anchor 1 & 2
    anchor1_peaks["Peak1_ID"] = anchor1_peaks["Peak1_Chr"].map(str) +':'+ (anchor1_peaks["Peak1_Start"]-1).map(str) +'-'+ anchor1_peaks["Peak1_End"].map(str)
    anchor2_peaks["Peak2_ID"] = anchor2_peaks["Peak2_Chr"].map(str) +':'+ (anchor2_peaks["Peak2_Start"]-1).map(str) +'-'+ anchor2_peaks["Peak2_End"].map(str)

    # Merging anchor points
    anchor1_peaks_anno =bed2D_anno.loc[:,['chr1', 's1','e1' ,'Entrez_ID_1', 'Gene_Name_1']].merge(anchor1_peaks.loc[:,['Peak1','Peak1_ID', 'Peak1_score']], left_index=True, right_index=True, how = 'left')
    anchor2_peaks_anno =bed2D_anno.loc[:,['chr2', 's2','e2' ,'Entrez_ID_2', 'Gene_Name_2','Q-Value_Bias','TSS_1', 'TSS_2']].merge(anchor2_peaks.loc[:,['Peak2','Peak2_ID', 'Peak2_score']], left_index=True, right_index=True, how = 'left')
    anchors_peaks_anno = anchor1_peaks_anno.merge(anchor2_peaks_anno, left_index=True, right_index=True,how = 'outer').drop_duplicates()

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

    for f in factor:
        anchors_peaks_anno[f + '_1'] = anchors_peaks_anno.apply (lambda row: peak_in_anchor_1(row), axis=1)
    for f in factor:
        anchors_peaks_anno[f + '_2'] = anchors_peaks_anno.apply (lambda row: peak_in_anchor_2(row), axis=1)

    anchors_peaks_anno.index.name = 'Interaction'

    # Creating dictionary with each factors as a key and associated df with interactions with factor overlap in at least one anchor point
    factor_dict={}
    for f in factor:
        factor_dict[f] = anchors_peaks_anno[(anchors_peaks_anno['Peak1'] == f) | (anchors_peaks_anno['Peak2'] == f)]
        factor_dict[f].loc[factor_dict[f].Peak1 !=f,['Peak1', 'Peak1_ID', 'Peak1_score']] = ''
        factor_dict[f].loc[factor_dict[f].Peak2 !=f,['Peak2', 'Peak2_ID', 'Peak2_score']] = ''
        factor_dict[f] = factor_dict[f].groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    anchors_peaks_anno = anchors_peaks_anno.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    # Saving annotated interactions files (all interaction and interactions with peak overlap)
    for f in factor:
        factor_dict[f].to_csv(str(f) + "_${prefix}_interactions.txt", index=False, sep='\t' )
    anchors_peaks_anno.to_csv("${prefix}_HOMER_annotated_interactions_with_peak_overlap.txt", index=True, sep='\t' )
    """

  else if (params.mode == 'differential')
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np

    # Column names for loaded data
    anchor1_peak_name = ('Anchor1_Chr', 'Anchor1_Start', 'Anchor1_End', 'Peak1_Chr', 'Peak1_Start', 'Peak1_End', 'Peak1_name', 'Peak1_score')
    anchor2_peak_name = ('Anchor2_Chr', 'Anchor2_Start', 'Anchor2_End', 'Peak2_Chr', 'Peak2_Start', 'Peak2_End', 'Peak2_name', 'Peak2_score')

    # Load interaction centered peak overlaps 2and annotated 2D-bed
    anchor1_peaks = pd.read_table("${anchor_1_peak_collect}", index_col=3,names=anchor1_peak_name).sort_index()
    anchor2_peaks = pd.read_table("${anchor_2_peak_collect}", index_col=3,names=anchor2_peak_name).sort_index()
    bed2D_anno = pd.read_table("${bed2D_index_anno}", index_col=1).sort_index().iloc[:,1:]

    # Load differntial peaks and merge
    peak_differntial = pd.read_table("${peak_differential}", index_col=0).sort_index()
    peak_differntial = peak_differntial.iloc[:, [${log2FC_column}-2, ${padj_column}-2]]
    peak_differntial.columns = ['log2FC', 'padj']
    anchor1_peaks=anchor1_peaks.merge(peak_differntial, left_on="Peak1_name", right_index=True, how = 'left')
    anchor2_peaks=anchor2_peaks.merge(peak_differntial, left_on="Peak2_name", right_index=True, how = 'left')

    # Create Peak columns (chr:start-end) for anchor 1 & 2
    anchor1_peaks["Peak1_ID"] = anchor1_peaks["Peak1_Chr"].map(str) +':'+ (anchor1_peaks["Peak1_Start"]-1).map(str) +'-'+ anchor1_peaks["Peak1_End"].map(str)
    anchor2_peaks["Peak2_ID"] = anchor2_peaks["Peak2_Chr"].map(str) +':'+ (anchor2_peaks["Peak2_Start"]-1).map(str) +'-'+ anchor2_peaks["Peak2_End"].map(str)

    # Merging anchor points
    anchor1_peaks_anno =bed2D_anno.loc[:,['chr1', 's1','e1' ,'Entrez_ID_1', 'Gene_Name_1']].merge(anchor1_peaks.loc[:,['Peak1_ID', 'Peak1_score', 'log2FC', 'padj']], left_index=True, right_index=True, how = 'left')
    anchor2_peaks_anno =bed2D_anno.loc[:,['chr2', 's2','e2' ,'Entrez_ID_2', 'Gene_Name_2','Q-Value_Bias','TSS_1', 'TSS_2']].merge(anchor2_peaks.loc[:,['Peak2_ID', 'Peak2_score', 'log2FC', 'padj']], left_index=True, right_index=True, how = 'left')
    anchors_peaks_anno = anchor1_peaks_anno.merge(anchor2_peaks_anno, left_index=True, right_index=True,how = 'outer', suffixes = ('_1', '_2')).drop_duplicates()

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
    anchors_peaks_anno_factor = anchors_peaks_anno[(anchors_peaks_anno['Overlap_1'] == 1) | (anchors_peaks_anno['Overlap_2'] == 1)]
    anchors_peaks_anno_factor = anchors_peaks_anno_factor.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    anchors_peaks_anno = anchors_peaks_anno.groupby('Interaction').agg(lambda x: ', '.join(filter(None, list(x.unique().astype(str)))))

    # Saving annotated interactions files (all interaction and interactions with peak overlap)
    anchors_peaks_anno_factor.to_csv("${peak_names}_${prefix}_interactions.txt", index=False, sep='\t' )
    anchors_peaks_anno.to_csv("${prefix}_HOMER_annotated_interactions_with_peak_overlap.txt", index=True, sep='\t' )
    """
}

/*
* 10. CREATES NETWORK FILES FOR CYTOSCAPE
*/
process NETWORK {
publishDir "${params.outdir}/tmp/process10", mode: 'copy', enabled: params.save_tmp
publishDir "${params.outdir}/Network", mode: 'copy', pattern: 'Network_*_interactions_*.txt', enabled: params.network | params.complete

when:
params.network | params.upset_plot | params.circos_plot | params.complete

input:
path interactions_annotated from ch_interactions_all
path genes from ch_genes
val prefix from Channel.value(params.prefix)
val peak_names from ch_t_3.peak_names.collect().map{ it2 -> it2.join(' ')}
val network_mode from Channel.value(params.network_mode)
val upset_plot from Channel.value(params.upset_plot)
val circos_plot from Channel.value(params.circos_plot)
val filter_genes from Channel.value(params.filter_genes)
val complete from Channel.value(params.complete)


output:
path "Network_Edges_${prefix}_interactions_*.txt" into ch_edges
path "Network_Nodes_${prefix}_interactions_*.txt" into ch_nodes

//For Upset PLOT (only created in multiple mode)
path "UpSet_${prefix}_interactions_Promoter.txt" optional true into ch_upset_promoter
path "UpSet_${prefix}_interactions_Distal.txt" optional true into ch_upset_distal
path "UpSet_${prefix}_interactions_Promoter_genes.txt" optional true into ch_upset_promoter_genes
path "UpSet_${prefix}_interactions_Distal_genes.txt" optional true into ch_upset_distal_genes

script:
if (params.mode == 'basic')
  """
  #!/usr/bin/env python

  import pandas as pd
  import numpy as np

  #Loading input file
  anchors_peaks_anno = pd.read_table("${interactions_annotated}", index_col=0)

  # Aggregating interaction file to only incude one row per interaction
  interactions_anno = anchors_peaks_anno.iloc[:,np.r_[0:5,7:15, 17:18]].groupby(by=['chr1', 's1', 'e1', 'Entrez_ID_1', 'Gene_Name_1', 'chr2', 's2', 'e2','Entrez_ID_2', 'Gene_Name_2', 'Q-Value_Bias', 'TSS_1', 'TSS_2'], axis=0, as_index=False).max()
  interactions_anno['Anchor1'] = interactions_anno["chr1"].map(str) +':'+ (interactions_anno["s1"]).map(str) +'-'+ interactions_anno["e1"].map(str)
  interactions_anno['Anchor2'] = interactions_anno["chr2"].map(str) +':'+ (interactions_anno["s2"]).map(str) +'-'+ interactions_anno["e2"].map(str)
  interactions_anno = pd.concat([interactions_anno['Anchor1'], interactions_anno.iloc[:,3:5], interactions_anno['Anchor2'],interactions_anno.iloc[:,8:(len(interactions_anno.columns)-2)]], axis=1)

  ### Creating edge table for cytoscape
  #Factor-Interaction
  Factor_Interaction = anchors_peaks_anno.copy(deep=True)
  Factor_Interaction.loc[Factor_Interaction.Overlap_1 == 1, 'Peak1'] = "${peak_names}"
  Factor_Interaction.loc[Factor_Interaction.Overlap_2 == 1, 'Peak2'] = "${peak_names}"
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

  network_mode = "${network_mode}"

  if network_mode == "factor":
    #Filter edges based on factor
    Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter['Target'])]
    Promoter_Promoter_filt_f = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter['Target'])]
    Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target'])]

  elif network_mode == "genes":
    #Filter edges based on gene
    genes = pd.read_table("${genes}", header=None)
    Promoter_Gene_filt_g = Promoter_Gene[Promoter_Gene['Target'].isin(genes.iloc[:,0])]
    Promoter_Promoter_filt_g = Promoter_Promoter[Promoter_Promoter['Source'].isin(Promoter_Gene_filt_g['Source']) | Promoter_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Distal_Promoter_filt_g = Distal_Promoter[Distal_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Factor_Promoter_filt_g = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Factor_Distal_filt_g = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g['Source'])]

  ### Creating edge table for cytoscape
  if network_mode == "all":
    Egdes_all = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Promoter, Promoter_Gene]).drop_duplicates()
    Egdes_all.to_csv("Network_Edges_${prefix}_interactions_all.txt", index=False, sep='\t' )

  elif network_mode == "factor":
    Egdes_factor =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
    Egdes_factor.to_csv("Network_Edges_${prefix}_interactions_factors.txt", index=False, sep='\t' )

  elif network_mode == "genes":
    Egdes_genes =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()
    Egdes_genes.to_csv("Network_Edges_${prefix}_interactions_genes.txt", index=False, sep='\t' )


  ### Creating node table for cytoscape

  if network_mode == "all":
    #Specifying node type for all nodes
    nodes_all = pd.DataFrame(pd.unique(Egdes_all[['Source', 'Target']].dropna().values.ravel('K')))
    nodes_all.columns=['Node']
    nodes_all['Node_type'] = np.where(nodes_all['Node'].isin(Factor_Distal['Source']) | nodes_all['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(nodes_all['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                       (np.where(nodes_all['Node'].isin(Distal_Promoter['Target']) | nodes_all['Node'].isin(Promoter_Promoter['Source']) | nodes_all['Node'].isin(Promoter_Promoter['Target']), 'Promoter',
                                          (np.where(nodes_all['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))
    nodes_all.to_csv("Network_Nodes_${prefix}_interactions_all.txt", index=False, sep='\t' )

  elif network_mode == "factor":
    # Specifying node type for all nodes that are associated with factor binding
    nodes_factor = pd.DataFrame(pd.unique(Egdes_factor[['Source', 'Target']].dropna().values.ravel('K')))
    nodes_factor.columns=['Node']
    nodes_factor['Node_type'] = np.where(nodes_factor['Node'].isin(Factor_Distal['Source']) | nodes_factor['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(nodes_factor['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                       (np.where(nodes_factor['Node'].isin(Distal_Promoter_filt_f['Target']) | nodes_factor['Node'].isin(Promoter_Promoter_filt_f['Source']) | nodes_factor['Node'].isin(Promoter_Promoter_filt_f['Target']), 'Promoter',
                                          (np.where(nodes_factor['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))
    nodes_factor.to_csv("Network_Nodes_${prefix}_interactions_factors.txt", index=False, sep='\t' )

  elif network_mode == "genes":
    # Specifying node type for all nodes that are associated with selected genes
    nodes_genes = pd.DataFrame(pd.unique(Egdes_genes[['Source', 'Target']].dropna().values.ravel('K')))
    nodes_genes.columns=['Node']
    nodes_genes['Node_type'] = np.where(nodes_genes['Node'].isin(Factor_Distal_filt_g['Source']) | nodes_genes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                    (np.where(nodes_genes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                       (np.where(nodes_genes['Node'].isin(Distal_Promoter_filt_g['Target']) | nodes_genes['Node'].isin(Promoter_Promoter_filt_g['Source']) | nodes_genes['Node'].isin(Promoter_Promoter_filt_g['Target']), 'Promoter',
                                          (np.where(nodes_genes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))
    nodes_genes.to_csv("Network_Nodes_${prefix}_interactions_genes.txt", index=False, sep='\t' )
  """

else if (params.mode == 'multiple')
  """
  #!/usr/bin/env python

  import pandas as pd
  import numpy as np

  #Loading input file
  anchors_peaks_anno = pd.read_table("${interactions_annotated}",index_col=0)

  # Aggregating interaction file to only incude one row per interaction
  interactions_anno = anchors_peaks_anno.iloc[:,np.r_[0:5,8:16, 19:len(anchors_peaks_anno.columns)]].groupby(by=['chr1', 's1', 'e1', 'Entrez_ID_1', 'Gene_Name_1', 'chr2', 's2', 'e2','Entrez_ID_2', 'Gene_Name_2', 'Q-Value_Bias', 'TSS_1', 'TSS_2'], axis=0, as_index=False).sum()
  interactions_anno['Anchor1'] = interactions_anno["chr1"].map(str) +':'+ (interactions_anno["s1"]).map(str) +'-'+ interactions_anno["e1"].map(str)
  interactions_anno['Anchor2'] = interactions_anno["chr2"].map(str) +':'+ (interactions_anno["s2"]).map(str) +'-'+ interactions_anno["e2"].map(str)
  interactions_anno = pd.concat([interactions_anno['Anchor1'], interactions_anno.iloc[:,3:5], interactions_anno['Anchor2'],interactions_anno.iloc[:,8:(len(interactions_anno.columns)-2)]], axis=1)

  # Factor-Interaction
  Factor_Interaction_all = anchors_peaks_anno[['chr1', 's1', 'e1','Gene_Name_1', 'Peak1','Peak1_ID', 'Peak1_score', 'chr2', 's2', 'e2',  'Gene_Name_2','Peak2','Peak2_ID','Peak2_score', 'TSS_1', 'TSS_2']]
  Factor_Interaction_all['Anchor1'] = Factor_Interaction_all['chr1'].map(str) +':'+ (Factor_Interaction_all['s1']).map(str) +'-'+ Factor_Interaction_all['e1'].map(str)
  Factor_Interaction_all['Anchor2'] = Factor_Interaction_all['chr2'].map(str) +':'+ (Factor_Interaction_all['s2']).map(str) +'-'+ Factor_Interaction_all['e2'].map(str)
  Factor_Interaction = Factor_Interaction_all.dropna(subset=['Peak1', 'Peak2'], thresh=1)

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

  network_mode = "${network_mode}"

  if network_mode == "factor":
    #Filter edges based on factor
    Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter['Target'])]
    Promoter_Promoter_filt_f = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter['Target'])]
    Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target'])]

  elif network_mode == "genes":
    #Filter edges based on gene
    genes = pd.read_table("${genes}", header=None)
    Promoter_Gene_filt_g = Promoter_Gene[Promoter_Gene['Target'].isin(genes.iloc[:,0])]
    Promoter_Promoter_filt_g = Promoter_Promoter[Promoter_Promoter['Source'].isin(Promoter_Gene_filt_g['Source']) | Promoter_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Distal_Promoter_filt_g = Distal_Promoter[Distal_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Factor_Promoter_filt_g = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
    Factor_Distal_filt_g = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g['Source'])]

  ### Creating edge table for cytoscape
  if network_mode == "all":
    Egdes_all = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Promoter, Promoter_Gene]).drop_duplicates()
    Egdes_all.to_csv("Network_Edges_${prefix}_interactions_all.txt", index=False, sep='\t' )

  elif network_mode == "factor":
    Egdes_factor =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
    Egdes_factor.to_csv("Network_Edges_${prefix}_interactions_factors.txt", index=False, sep='\t' )

  elif network_mode == "genes":
    Egdes_genes =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()
    Egdes_genes.to_csv("Network_Edges_${prefix}_interactions_genes.txt", index=False, sep='\t' )


  ### Creating node table for cytoscape

  if network_mode == "all":
    #Specifying node type for all nodes
    nodes_all = pd.DataFrame(pd.unique(Egdes_all[['Source', 'Target']].dropna().values.ravel('K')))
    nodes_all.columns=['Node']
    nodes_all['Node_type'] = np.where(nodes_all['Node'].isin(Factor_Distal['Source']) | nodes_all['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(nodes_all['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                       (np.where(nodes_all['Node'].isin(Distal_Promoter['Target']) | nodes_all['Node'].isin(Promoter_Promoter['Source']) | nodes_all['Node'].isin(Promoter_Promoter['Target']), 'Promoter',
                                          (np.where(nodes_all['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))
    nodes_all.to_csv("Network_Nodes_${prefix}_interactions_all.txt", index=False, sep='\t' )

  elif network_mode == "factor":
    # Specifying node type for all nodes that are associated with factor binding
    nodes_factor = pd.DataFrame(pd.unique(Egdes_factor[['Source', 'Target']].dropna().values.ravel('K')))
    nodes_factor.columns=['Node']
    nodes_factor['Node_type'] = np.where(nodes_factor['Node'].isin(Factor_Distal['Source']) | nodes_factor['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(nodes_factor['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                       (np.where(nodes_factor['Node'].isin(Distal_Promoter_filt_f['Target']) | nodes_factor['Node'].isin(Promoter_Promoter_filt_f['Source']) | nodes_factor['Node'].isin(Promoter_Promoter_filt_f['Target']), 'Promoter',
                                          (np.where(nodes_factor['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))
    nodes_factor.to_csv("Network_Nodes_${prefix}_interactions_factors.txt", index=False, sep='\t' )

  elif network_mode == "genes":
    # Specifying node type for all nodes that are associated with selected genes
    nodes_genes = pd.DataFrame(pd.unique(Egdes_genes[['Source', 'Target']].dropna().values.ravel('K')))
    nodes_genes.columns=['Node']
    nodes_genes['Node_type'] = np.where(nodes_genes['Node'].isin(Factor_Distal_filt_g['Source']) | nodes_genes['Node'].isin(Factor_Promoter_filt_g['Source']), 'Factor',
                                    (np.where(nodes_genes['Node'].isin(Distal_Promoter_filt_g['Source']), 'Distal',
                                       (np.where(nodes_genes['Node'].isin(Distal_Promoter_filt_g['Target']) | nodes_genes['Node'].isin(Promoter_Promoter_filt_g['Source']) | nodes_genes['Node'].isin(Promoter_Promoter_filt_g['Target']), 'Promoter',
                                          (np.where(nodes_genes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))
    nodes_genes.to_csv("Network_Nodes_${prefix}_interactions_genes.txt", index=False, sep='\t' )

  ### Save files for UpSet PLOT
  plot_upset = "${upset_plot}"
  circos_plot = "${circos_plot}"
  complete = "${complete}"
  if (plot_upset == 'true' or complete == 'true' or circos_plot == 'true'):
    Factor_Promoter.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv("UpSet_${prefix}_interactions_Promoter.txt", index=False, sep='\t' )
    Factor_Distal.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv("UpSet_${prefix}_interactions_Distal.txt", index=False, sep='\t' )

    filter_genes = "${filter_genes}"
    if filter_genes == 'true':
      Factor_Promoter_filt_g.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv("UpSet_${prefix}_interactions_Promoter_genes.txt", index=False, sep='\t' )
      Factor_Distal_filt_g.loc[:,['Source', 'Target', 'Edge_type']].sort_index().reset_index().drop_duplicates().to_csv("UpSet_${prefix}_interactions_Distal_genes.txt", index=False, sep='\t' )
  """
}

if ({params.upset_plot | params.circos_plot | params.complete} && !params.filter_genes ){
  ch_upset_promoter_genes = file("No_input_upset_promoter_genes")
  ch_upset_distal_genes = file("No_input_upset_distal_genes")
}

/*
 * 11. UPSET PLOT FOR FACTOR BINDING IN ANCHOR POINTS - ALL INTERACTIONS WITH FACTOR AND INTERACTIONS FILTERES FOR GENES
 */
process UPSET_PLOT {
  publishDir "${params.outdir}/tmp/process11", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Figures/Upset", mode: 'copy', pattern: 'Upset_plot_*.pdf', enabled: params.upset_plot | params.complete

  when:
  params.upset_plot | params.circos_plot | params.complete

  input:
  path upset_promoter from ch_upset_promoter
  path upset_distal from ch_upset_distal
  path upset_promoter_g from ch_upset_promoter_genes
  path upset_distal_g from ch_upset_distal_genes
  val prefix from Channel.value(params.prefix)
  val circos_plot from Channel.value(params.circos_plot)
  val filter_genes from Channel.value(params.filter_genes)
  val complete from Channel.value(params.complete)

  output:
  //Upset plots
  path "Upset_plot_Promoter_all.pdf" into ch_upset_plot_promoter_all
  path "Upset_plot_Distal_all.pdf" into ch_upset_plot_distal_all
  path "Upset_plot_Promoter_genelist.pdf" optional true into ch_upset_plot_promoter_genes
  path "Upset_plot_Distal_genelist.pdf" optional true into ch_upset_plot_distal_genes

  //For circos plot
  path "Circos_factors_${prefix}_interactions.txt" optional true into ch_circos_f
  path "Circos_genes_${prefix}_interactions.txt" optional true into ch_circos_g


  script:
  """
  #!/usr/bin/env python

  import pandas as pd
  import numpy as np
  from upsetplot import plot
  import matplotlib.pyplot as plt

  ### Loading and organizing data
  upset_promoter = pd.read_table("${upset_promoter}")
  upset_distal = pd.read_table("${upset_distal}")

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

  filter_genes = "${filter_genes}"
  if filter_genes == 'true':
    upset_promoter_g = pd.read_table("${upset_promoter_g}")
    upset_distal_g = pd.read_table("${upset_distal_g}")

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

  ### Preperations for circis PLOTS
  circos_plot = "${circos_plot}"
  complete = "${complete}"
  if (complete == 'true' or circos_plot == 'true'):
    upset_promoter = upset_promoter.iloc[:,np.r_[0,4:len(upset_promoter.columns)]]
    upset_promoter = upset_promoter.groupby(upset_promoter.columns[0]).max()
    upset_promoter['promoter_cat'] = 'Promoter'+upset_promoter.eq(True).dot('_'+upset_promoter.columns)
    upset_distal = upset_distal.iloc[:,np.r_[0,4:len(upset_distal.columns)]]
    upset_distal = upset_distal.groupby(upset_distal.columns[0]).max()
    upset_distal['distal_cat'] = 'distal'+upset_distal.eq(True).dot('_'+upset_distal.columns)

    circos_f = upset_promoter.merge(upset_distal, left_index=True, right_index=True, how = 'outer')
    circos_f.fillna(value={'promoter_cat': 'Promoter_NoBinding', 'distal_cat': 'Distal_NoBinding'}, inplace=True)
    circos_f.fillna(False,inplace=True)
    circos_f = circos_f.groupby(list(circos_f.columns)).size().to_frame('size').reset_index()
    circos_f.to_csv("Circos_fenes_${prefix}_interactions.txt", index=False, sep='\t' )

    filter_genes = "${filter_genes}"
    if filter_genes == 'true':
      upset_promoter_g = upset_promoter_g.iloc[:,np.r_[0,4:len(upset_promoter_g.columns)]]
      upset_promoter_g = upset_promoter_g.groupby(upset_promoter_g.columns[0]).max()
      upset_promoter_g['promoter_cat'] = 'Promoter'+upset_promoter_g.eq(True).dot('_'+upset_promoter_g.columns)
      upset_distal_g = upset_distal_g.iloc[:,np.r_[0,4:len(upset_distal_g.columns)]]
      upset_distal_g = upset_distal_g.groupby(upset_distal_g.columns[0]).max()
      upset_distal_g['distal_cat'] = 'distal'+upset_distal_g.eq(True).dot('_'+upset_distal_g.columns)

      circos_g = upset_promoter_g.merge(upset_distal_g, left_index=True, right_index=True, how = 'outer')
      circos_g.fillna(value={'promoter_cat': 'Promoter_NoBinding', 'distal_cat': 'Distal_NoBinding'}, inplace=True)
      circos_g.fillna(False,inplace=True)
      circos_g = circos_g.groupby(list(circos_g.columns)).size().to_frame('size').reset_index()
      circos_g.to_csv("Circos_genes_${prefix}_interactions.txt", index=False, sep='\t' )
  """
}

if ({params.circos_plot | params.complete} && !params.filter_genes ){
  ch_circos_g = file("No_input_ch_circos_g")
}

/*
 * 12. CIRCOS PLOTS - ALL INTERACTIONS WITH FACTOR AND INTERACTIONS FILTERED FOR GENES
 */
process CIRCOS_PLOT {
  publishDir "${params.outdir}/tmp/process12", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Figures/Circos", mode: 'copy', enabled: params.circos_plot | params.complete

  when:
  params.circos_plot | params.complete

  input:
  path circos_f from ch_circos_f
  path circos_g from ch_circos_g
  val filter_genes from Channel.value(params.filter_genes)


  output:
  path  "Circos_plot_*.pdf" into ch_circos_plot


  script:
  """
  #!/usr/local/bin/Rscript --vanilla

  require(circlize)
  require(viridis)
  require(inlmisc)
  require(stringr)
  require(mgsub)

  ## Regions with factor in at least one anchor point
  circos_data_all <- read.table("${circos_f}", header=TRUE, sep="\t")
  circos_data_all_2 <- circos_data_all[,c("promoter_cat", "distal_cat", "size")]
  nf_all <- (ncol(circos_data_all)-3)/2
  circos_data_all_p <- unique(circos_data_all[,1:nf_all])
  colnames(circos_data_all_p) <- sub("_x", "", colnames(circos_data_all_p))
  circos_data_all_d <- unique(circos_data_all[,(nf_all+2):(nf_all*2+1)])
  colnames(circos_data_all_d) <- sub("_y", "", colnames(circos_data_all_d))
  circos_data_all_pd <- rbind(circos_data_all_p,circos_data_all_d)
  circos_data_all_pd <- mgsub(circos_data_all_pd, c("True", "False"), c("black", "white"))

  np_all <- nrow(circos_data_all_p)
  nd_all <- nrow(circos_data_all_d)

  cols_all <- c(GetColors(n = np_all, start = 0.2, end = 0.9), GetColors(n = nd_all,start = 0.2, end = 0.9))

  factor_anno <- list()
  for (f in 1:nf_all){
    factor_anno[[f]] <- list(track.height = 0.05, bg.border = "black", bg.col=circos_data_all_pd[,f])
  }
  x=1.82
  pdf("Circos_plot_factors.pdf")
  circos.par(start.degree = 0)
  chordDiagram(circos_data_all_2, big.gap = 25, directional = 1, grid.col = cols_all, transparency = 0.5,annotationTrack = "grid", grid.border="black", annotationTrackHeight=0.05,
      preAllocateTracks = factor_anno, xmax=0.1)
  for (n in rev(colnames(circos_data_all_pd))){
    circos.text(6, x, n, facing="bending.inside", cex=0.75)
    x=x+1.42
  }
  dev.off()

  ## Regions associated with genelist
  filter_genes = "${filter_genes}"
  if (filter_genes == 'true'){
    circos_data_genes <- read.table("${circos_g}", header=TRUE, sep="\t")

    circos_data_genes_2 <- circos_data_genes[,c("promoter_cat", "distal_cat", "size")]
    nf_genes<- (ncol(circos_data_genes)-3)/2
    circos_data_genes_p <- unique(circos_data_genes[,1:nf_genes])
    colnames(circos_data_genes_p) <- sub("_x", "", colnames(circos_data_genes_p))
    circos_data_genes_d <- unique(circos_data_genes[,(nf_genes+2):(nf_genes*2+1)])
    colnames(circos_data_genes_d) <- sub("_y", "", colnames(circos_data_genes_d))
    circos_data_genes_pd <- rbind(circos_data_genes_p,circos_data_genes_d)
    circos_data_genes_pd <- mgsub(circos_data_genes_pd, c("True", "False"), c("black", "white"))

    np_genes <- nrow(circos_data_genes_p)
    nd_genes <- nrow(circos_data_genes_d)

    cols_genes <- c(GetColors(n = np_genes, start = 0.2, end = 0.9), GetColors(n = nd_genes,start = 0.2, end = 0.9))

    factor_anno <- list()
    for (f in 1:nf_genes){
      factor_anno[[f]] <- list(track.height = 0.05, bg.border = "black", bg.col=circos_data_genes_pd[,f])
    }
    x=1.82
    pdf("Circos_plot_genelist.pdf")
    circos.par(start.degree = 0)
    chordDiagram(circos_data_genes_2, big.gap = 25, directional = 1, grid.col = cols_genes, transparency = 0.5,annotationTrack = "grid", grid.border="black", annotationTrackHeight=0.05,
        preAllocateTracks = factor_anno)
    for (n in rev(colnames(circos_data_genes_pd))){
      circos.text(6, x, n, facing="bending.inside", cex=0.75)
      x=x+1.42
    }
    dev.off()
  }
  """
}
