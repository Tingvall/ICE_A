#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  Usage:
  The typical command for running the pipeline is as follows:
    nextflow run interaction_anno_2.nf --bed2D interactions.bed  --genome mm10 --peaks peaks.txt

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
    --promoter_distance             Distance +/- TSS considered a promoter (default: 2500).
    --binsize                       Binsize used for interactions (default: 5000).
    --interaction_threshold          Lower interaction distance threshold, regions with a distance to the closest TSS < interaction_threshold will be proximity annotated (default: 10000).
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

    Arguments - Differntial mode specific:
    --peak_differntial [file]       Path to textfile that contain log2FC and adjusted p-value from differntial analyis. The 1st column should contain peakID matching the peakID in the 4th column of the input bed file. Standard DESeq2 output is expected (with log2FC in the 3rd column and padj in the 9th column), but other formats are accepted as well is the column corresponding to log2FC and padj are specified with the aruguments --log2FC_column and --padj column.
    --log2FC_column 	              Specifies which column in --peak_differential that contain the log2FC values. Deafult: 3 (standard DESEq2 output).
    --padj_column 	                Specifies which column in --peak_differential that contain the adjusted p-value values. Deafult: 9 (standard DESEq2 output).
    --log2FC 	                      Log2FC threshold for differntial peaks. Default: 1.5
    --padj 	                        Adjusted p-value threshold for differntial peaks. Default: 0.05
    --expression [file]
    --skip_expression 	            Use this argumnet if no --expression file is provided.
    --expression_log2FC_column 	    Specifies which column in --expression that contain the log2FC values. Deafult: 3 (standard DESEq2 output).
    --expression_padj_column 	      Specifies which column in --expression that contain the adjusted p-value values. Deafult: 9 (standard DESEq2 output).
    --expression_log2FC 	          Set the log2FC threshold for differntial genes. Default: 1.5
    --expression_padj 	            Set the adjusted p-value threshold for differntial genes. Default: 0.05
    --expression [file]             Only used in differntial mode when `--skip_expression` is false (default). Specifies path to file that contain information about differntial expression between the two conditions. The first column must contain gene symbol. The column for log2FC/padj can be specified using `--expression_log2FC_column` and `--expression_padj_column`respectivly, default: 3 & 9 (stadanrds DESeq2 format). Note: make sure that you use the same direction for the comparison in `--peak_differential` and `--expression`.

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
        .map { row -> [ row.sample, [ file(row.path) ] ] }
        .set { ch_peaks_split}

if (!params.genome)      { exit 1, 'Refence genome not specified' }


if (params.network_mode == 'genes') {
  if (params.genes)     { ch_genes = Channel.fromPath(params.genes, checkIfExists: true) } else { exit 1, 'Genes not specified' }
}
else {
  ch_genes = file(params.genes)
}

if (params.mode == 'differential') {
  if (params.peak_differential)     { ch_peak_differential = Channel.fromPath(params.peak_differential, checkIfExists: true) } else { exit 1, 'Peak file log2FC and adjusted p-value not provided' }
  ch_peak_differential.into { ch_peak_differential_1; ch_peak_differential_2; ch_peak_differential_3 }
  if (!params.skip_expression) {
    if (params.expression)     { ch_expression = Channel.fromPath(params.expression, checkIfExists: true) } else { exit 1, 'Expression file not provided' }
    ch_expression.into {ch_expression_1; ch_expression_2}
  }
  else {ch_expression = file(params.expression)
    ch_expression_1 = file(params.expression)
    ch_expression_2 = file(params.expression)
  }
}
else {
  ch_peak_differential_1 = file(params.peak_differential)
  ch_peak_differential_2 = file(params.peak_differential)
  ch_peak_differential_3 = file(params.peak_differential)

  ch_expression_1 = file(params.expression)
  ch_expression_2 = file(params.expression)

}

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
    publishDir "${params.outdir}/Interaction_annotation/All_interactions", mode: 'copy', enabled: params.annotate_interactions | params.complete

    when:
    !params.skip_anno

    input:
    path anchor1_anno from ch_anchor1_anno
    path anchor2_anno from ch_anchor2_anno
    path bed2D_index from ch_bed2D_index
    val prefix from Channel.value(params.prefix)
    val binsize from Channel.value(params.binsize)

    output:
    path "${prefix}_HOMER_annotated_interactions.txt" into ch_bed2D_anno

    script:
    """
    join_annotated_interactions.py ${anchor1_anno} ${anchor2_anno} ${bed2D_index} --prefix=${prefix} --binsize=${binsize}
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

if (!params.interaction_threshold){
  ch_interaction_threshold = Channel.value(2*params.binsize)
}
else{
  ch_interaction_threshold = Channel.value(params.interaction_threshold)
}
  /*
   * 7. PEAK ANNOTATION: PEAKS ANNOTATED BY PROXIMITY OR INTERACION-BASED ANNOTATION
   */
  process PEAK_INTERACTION_BASED_ANNOTATION {
    publishDir "${params.outdir}/tmp/process7", mode: 'copy', enabled: params.save_tmp
    publishDir "${params.outdir}/Peak_annotation/${peak_name}", mode: 'copy'

    input:
    set val(peak_name), file(peak_anno_anchor1), file(peak_anno_anchor2), file(peak_anno), file(bed2D_index_anno) from ch_peak_anno_anchor1.join(ch_peak_anno_anchor2).join(ch_peak_anno).combine(ch_bed2D_index_anno_1)
    val proximity_unannotated from Channel.value(params.proximity_unannotated)
    val multiple_anno from Channel.value(params.multiple_anno)
    val prefix from Channel.value(params.prefix)
    val mode from Channel.value(params.mode)
    val promoter_distance from Channel.value(params.promoter_distance)
    val interaction_threshold from ch_interaction_threshold

    //Differntial mode specific
    path peak_differential from ch_peak_differential_1
    val log2FC_column from Channel.value(params.log2FC_column)
    val padj_column from Channel.value(params.padj_column)
    val log2FC from Channel.value(params.log2FC)
    val padj from Channel.value(params.padj)
    val skip_expression from Channel.value(params.skip_expression)

    output:
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated.txt") into ch_peak_PLACseq_annotated
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_up.txt") optional true into ch_peak_PLACseq_annotated_up
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_down.txt") optional true into ch_peak_PLACseq_annotated_down
    tuple val(peak_name), path( "${peak_name}_${prefix}_annotated_for_differential_expression.txt") optional true into ch_peak_PLACseq_annotated_for_differntial_expression

    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_genelist.txt") into ch_peak_PLACseq_annotated_genelist
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_genelist_up.txt") optional true into ch_peak_PLACseq_annotated_genelist_up
    tuple val(peak_name), path("${peak_name}_${prefix}_annotated_genelist_down.txt") optional true into ch_peak_PLACseq_annotated_genelist_down


    script:
    """
    peak_annotation.py ${peak_anno_anchor1} ${peak_anno_anchor2} ${peak_anno} ${bed2D_index_anno} --peak_name ${peak_name} --prefix ${prefix} --proximity_unannotated ${proximity_unannotated} --mode ${mode} --multiple_anno ${multiple_anno} --promoter_distance ${promoter_distance} --interaction_threshold ${interaction_threshold} --peak_differential ${peak_differential} --log2FC_column ${log2FC_column} --padj_column ${padj_column} --log2FC ${log2FC} --padj ${padj} --skip_expression ${skip_expression}
    """
}

def criteria = multiMapCriteria {
                     sample: it[0]
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
  val sample from ch_t_1.sample.collect().map{ it2 -> it2.join(' ')}
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
    bedtools intersect -wa -wb -a $bed2D_anno_split_anchor1 -b $peak_beds -names $sample > Anchor_1_peak_collect.bed
    bedtools intersect -wa -wb -a $bed2D_anno_split_anchor2 -b $peak_beds -names $sample > Anchor_2_peak_collect.bed
    """
}


/*
 * 9. INTERACTION CENTERED ANNOTATION - WITH PEAK OVERLAP
 */
process ANNOTATE_INTERACTION_WITH_PEAKS {
  publishDir "${params.outdir}/tmp/process9", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Interaction_annotation/Peak_sepcific_interactions", mode: 'copy', pattern: '*{interactions, up, down}.txt', enabled: params.annotate_interactions | params.complete
  publishDir "${params.outdir}/Interaction_annotation/All_interactions", mode: 'copy', pattern: '*overlap.txt', enabled: params.annotate_interactions | params.complete

  when:
  params.annotate_interactions | params.network | params.upset_plot | params.circos_plot | params.complete

  input:
  path anchor_1_peak_collect from ch_anchor_1_peak_collect
  path anchor_2_peak_collect from ch_anchor_2_peak_collect
  path bed2D_index_anno from ch_bed2D_index_anno_2
  val prefix from Channel.value(params.prefix)
  val sample from ch_t_2.sample.collect().map{ it2 -> it2.join(' ')}
  val network from Channel.value(params.network)
  val complete from Channel.value(params.complete)

  //Multiple mode specific
  val upset_plot from Channel.value(params.upset_plot)
  val circos_plot from Channel.value(params.circos_plot)

  //Differntial mode specific
  path peak_differential from ch_peak_differential_2
  val log2FC_column from Channel.value(params.log2FC_column)
  val padj_column from Channel.value(params.padj_column)
  val log2FC from Channel.value(params.log2FC)
  val padj from Channel.value(params.padj)

  output:
  path "*_${prefix}_interactions.txt" into ch_interactions_by_factor
  path "${prefix}_HOMER_annotated_interactions_with_peak_overlap.txt" into ch_interactions_all
  path "${prefix}_HOMER_annotated_interactions_with_peak_overlap_not_aggregated.txt"  optional true into ch_interactions_all_not_aggregated
  path "*_${prefix}_interactions_up.txt" optional true into ch_interactions_up
  path "*_${prefix}_interactions_down.txt" optional true into ch_interactions_down


  script:
  if (params.mode == 'basic')
    """
    interaction_annotation_basic.py ${anchor_1_peak_collect} ${anchor_2_peak_collect} ${bed2D_index_anno} --prefix ${prefix} --sample ${sample} --network ${network} --complete ${complete}
    """

  else if (params.mode == 'multiple')
    """
    interaction_annotation_multiple.py ${anchor_1_peak_collect} ${anchor_2_peak_collect} ${bed2D_index_anno} --prefix ${prefix} --network ${network} --complete ${complete} --upset_plot ${upset_plot} --circos_plot ${circos_plot}
    """

  else if (params.mode == 'differential')
    """
    interaction_annotation_multiple.py ${anchor_1_peak_collect} ${anchor_2_peak_collect} ${bed2D_index_anno} --prefix ${prefix} --sample ${sample} --network ${network} --complete ${complete} --peak_differential ${peak_differential} --log2FC_column ${log2FC_column} --padj_column ${padj_column} --log2FC ${log2FC} --padj ${padj}
    """
}

/*
* 10. CREATES NETWORK FILES FOR CYTOSCAPE
*/
process NETWORK {
publishDir "${params.outdir}/tmp/process10", mode: 'copy', enabled: params.save_tmp
publishDir "${params.outdir}/Network/Files", mode: 'copy', pattern: 'Network_*_interaction*.txt', enabled: params.network | params.complete

when:
params.network | params.upset_plot | params.circos_plot | params.complete

input:
path interactions_annotated from ch_interactions_all
path interactions_annotated_not_aggregated from ch_interactions_all_not_aggregated
path genes from ch_genes
val prefix from Channel.value(params.prefix)
val sample from ch_t_3.sample.collect().map{ it2 -> it2.join(' ')}
val network_mode from Channel.value(params.network_mode)
val complete from Channel.value(params.complete)
val promoter_promoter from Channel.value(params.promoter_promoter)

//Multiple mode specific
val upset_plot from Channel.value(params.upset_plot)
val circos_plot from Channel.value(params.circos_plot)
val filter_genes from Channel.value(params.filter_genes)

//Differntial mode specific
path peak_differential from ch_peak_differential_3
path expression from ch_expression_1
val log2FC_column from Channel.value(params.log2FC_column)
val padj_column from Channel.value(params.padj_column)
val log2FC from Channel.value(params.log2FC)
val padj from Channel.value(params.padj)
val skip_expression from Channel.value(params.skip_expression)
val expression_padj from Channel.value(params.expression_padj)
val expression_log2FC from Channel.value(params.expression_log2FC)
val expression_padj_column from Channel.value(params.expression_padj_column)
val expression_log2FC_column from Channel.value(params.expression_log2FC_column)


output:
path "Network_Edges_${prefix}_interactions.txt" into ch_edges
path "Network_Edges_${prefix}_interactions_up.txt" optional true into ch_edges_up
path "Network_Edges_${prefix}_interactions_down.txt" optional true into ch_edges_down

path "Network_Nodes_${prefix}_interactions.txt" into ch_nodes
path "Network_Nodes_${prefix}_interactions_up.txt" optional true into ch_nodes_up
path "Network_Nodes_${prefix}_interactions_down.txt" optional true into ch_nodes_down

//For Upset PLOT (only created in multiple mode)
path "UpSet_${prefix}_interactions_Promoter.txt" optional true into ch_upset_promoter
path "UpSet_${prefix}_interactions_Distal.txt" optional true into ch_upset_distal
path "UpSet_${prefix}_interactions_Promoter_genes.txt" optional true into ch_upset_promoter_genes
path "UpSet_${prefix}_interactions_Distal_genes.txt" optional true into ch_upset_distal_genes

script:
if (params.mode == 'basic')
  """
  network_preprocessing_basic.py ${interactions_annotated} ${interactions_annotated_not_aggregated} ${genes} --prefix ${prefix} --sample ${sample} --network_mode ${network_mode} --promoter_promoter ${promoter_promoter} --complete ${complete}
  """

else if (params.mode == 'multiple')
  """
  #!/usr/bin/env python

  import pandas as pd
  import numpy as np

  #Loading input file
  anchors_peaks_anno = pd.read_table("${interactions_annotated_not_aggregated}", index_col=0)
  interactions_anno = pd.read_table("${interactions_annotated}", index_col=0)

  # Aggregating interaction file to only incude one row per interaction
  interactions_anno = interactions_anno.iloc[:,np.r_[0:5,8:16, 19:len(anchors_peaks_anno.columns)]]
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
  promoter_promoter="${promoter_promoter}"

  if network_mode == "factor":
  #Filter edges based on factor
      Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target'])]
      if promoter_promoter =="true":
          Promoter_Promoter_filt_f = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter['Target'])]
          Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Target'])]
      else:
          Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Source'])]

  elif network_mode == "genes":
    #Filter edges based on gene
    genes = pd.read_table("${genes}", header=None)
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
      Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
    else:
      Edges =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()

  elif network_mode == "genes":
    if promoter_promoter =="true":
      Edges =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()
    else:
      Edges =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()

  Edges.to_csv("Network_Edges_${prefix}_interactions.txt", index=False, sep='\t' )


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
                                       (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))
    else:
      Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                   (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Factor_Promoter['Target']), 'Promoter',
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
                                     (np.where(Nodes['Node'].isin(Distal_Promoter_filt_g['Target'])  | Nodes['Node'].isin(Factor_Promoter_filt_g['Target']), 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene_filt_g['Target']), 'Gene', np.nan)))))))

  Nodes.to_csv("Network_Nodes_${prefix}_interactions.txt", index=False, sep='\t' )

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

else if (params.mode == 'differential')
  """
  #!/usr/bin/env python

  import pandas as pd
  import numpy as np

  #Loading input file
  anchors_peaks_anno = pd.read_table("${interactions_annotated_not_aggregated}", index_col=0)
  interactions_anno = pd.read_table("${interactions_annotated}", index_col=0)

  skip_expression = "${skip_expression}"
  if skip_expression == "false":
    expression = pd.read_table("${expression}", index_col=0).sort_index()
    expression = expression.iloc[:, [${expression_log2FC_column}-2, ${expression_padj_column}-2]]
    expression.columns = ['log2FC', 'padj']
    expression_diff = expression.loc[(expression['padj'] <= ${padj}) & (abs(expression['log2FC']) >= ${log2FC}),:]

  # Aggregating interaction file to only incude one row per interaction
  interactions_anno = interactions_anno.iloc[:,np.r_[0:5,9:17, 21:len(anchors_peaks_anno.columns)]]
  interactions_anno['Anchor1'] = interactions_anno["chr1"].map(str) +':'+ (interactions_anno["s1"]).map(str) +'-'+ interactions_anno["e1"].map(str)
  interactions_anno['Anchor2'] = interactions_anno["chr2"].map(str) +':'+ (interactions_anno["s2"]).map(str) +'-'+ interactions_anno["e2"].map(str)
  interactions_anno = pd.concat([interactions_anno['Anchor1'], interactions_anno.iloc[:,3:5], interactions_anno['Anchor2'],interactions_anno.iloc[:,8:(len(interactions_anno.columns)-2)]], axis=1)

  ### Creating edge table for cytoscape
  #Factor-Interaction
  Factor_Interaction = anchors_peaks_anno.copy(deep=True)
  Factor_Interaction.loc[Factor_Interaction.Overlap_1 == 1, 'Peak1'] = "${sample}"
  Factor_Interaction.loc[Factor_Interaction.Overlap_2 == 1, 'Peak2'] = "${sample}"
  Factor_Interaction = Factor_Interaction[['chr1', 's1', 'e1','Gene_Name_1', 'Peak1','Peak1_ID','Peak1_score', 'log2FC_1', 'padj_1','chr2', 's2', 'e2',  'Gene_Name_2','Peak2','Peak2_ID','Peak2_score','log2FC_2', 'padj_2', 'TSS_1', 'TSS_2']]
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

  Factor_Distal_1_diff = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 0) & (Factor_Interaction['TSS_2'] == 1) & (Factor_Interaction['padj_1'] <= ${padj}) & (abs(Factor_Interaction['log2FC_1']) >= ${log2FC}), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Distal_1_diff.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal_2_diff = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['TSS_2'] == 0) & (Factor_Interaction['padj_2'] <= ${padj}) & (abs(Factor_Interaction['log2FC_2']) >= ${log2FC}), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Distal_2_diff.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal_diff = Factor_Distal_1_diff.append(Factor_Distal_2_diff)
  Factor_Distal_diff['Edge_type'] = 'Factor-Distal'

  Factor_Distal_1_up = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 0) & (Factor_Interaction['TSS_2'] == 1) & (Factor_Interaction['padj_1'] <= ${padj}) & (Factor_Interaction['log2FC_1'] >= ${log2FC}), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Distal_1_up.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal_2_up = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['TSS_2'] == 0) & (Factor_Interaction['padj_2'] <= ${padj}) & (Factor_Interaction['log2FC_2'] >= ${log2FC}), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Distal_2_up.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal_up = Factor_Distal_1_up.append(Factor_Distal_2_up)
  Factor_Distal_up['Edge_type'] = 'Factor-Distal'

  Factor_Distal_1_down = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 0) & (Factor_Interaction['TSS_2'] == 1) & (Factor_Interaction['padj_1'] <= ${padj}) & (Factor_Interaction['log2FC_1'] <= -${log2FC}), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Distal_1_down.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal_2_down = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['TSS_2'] == 0) & (Factor_Interaction['padj_2'] <= ${padj}) & (Factor_Interaction['log2FC_2'] <= -${log2FC}), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Distal_2_down.columns = ['Source', 'Target', 'Edge_score']
  Factor_Distal_down = Factor_Distal_1_down.append(Factor_Distal_2_down)
  Factor_Distal_down['Edge_type'] = 'Factor-Distal'

  #Factor-Promoter
  Factor_Promoter_1 = Factor_Interaction.loc[Factor_Interaction['TSS_1'] == 1, ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Promoter_1.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_2 = Factor_Interaction.loc[Factor_Interaction['TSS_2'] == 1, ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Promoter_2.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter = Factor_Promoter_1.append(Factor_Promoter_2)
  Factor_Promoter['Edge_type'] = 'Factor-Promoter'

  Factor_Promoter_1_diff = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['padj_1'] <= ${padj}) & ((abs(Factor_Interaction['log2FC_1']) >= ${log2FC})), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Promoter_1_diff.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_2_diff = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['padj_1'] <= ${padj}) & ((abs(Factor_Interaction['log2FC_1']) >= ${log2FC})), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Promoter_2_diff.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_diff = Factor_Promoter_1_diff.append(Factor_Promoter_2_diff)
  Factor_Promoter_diff['Edge_type'] = 'Factor-Promoter'

  Factor_Promoter_1_up = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['padj_1'] <= ${padj}) & (Factor_Interaction['log2FC_1'] >= ${log2FC}), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Promoter_1_up.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_2_up = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['padj_2'] <= ${padj}) & (Factor_Interaction['log2FC_2'] >= ${log2FC}), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Promoter_2_up.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_up = Factor_Promoter_1_up.append(Factor_Promoter_2_up)
  Factor_Promoter_up['Edge_type'] = 'Factor-Promoter'

  Factor_Promoter_1_down = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['padj_1'] <= ${padj}) & (Factor_Interaction['log2FC_1'] <= -${log2FC}), ['Peak1',  'Anchor1', 'Peak1_score']].dropna(subset=['Peak1']).drop_duplicates()
  Factor_Promoter_1_down.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_2_down = Factor_Interaction.loc[(Factor_Interaction['TSS_1'] == 1) & (Factor_Interaction['padj_2'] <= ${padj}) & (Factor_Interaction['log2FC_2'] <= -${log2FC}), ['Peak2',  'Anchor2', 'Peak2_score']].dropna(subset=['Peak2']).drop_duplicates()
  Factor_Promoter_2_down.columns = ['Source', 'Target', 'Edge_score']
  Factor_Promoter_down = Factor_Promoter_1_down.append(Factor_Promoter_2_down)
  Factor_Promoter_down['Edge_type'] = 'Factor-Promoter'

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
  promoter_promoter="${promoter_promoter}"

  if network_mode == "factor":
  #Filter edges based on factor
      Distal_Promoter_filt_f = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal['Target'])]
      promoter_promoter="${promoter_promoter}"
      if promoter_promoter =="true":
          Promoter_Promoter_filt_f = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter['Target'])]
          Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Source'])| Promoter_Gene['Source'].isin(Promoter_Promoter_filt_f['Target'])]
      else:
          Promoter_Gene_filt_f = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_f['Source'])]

  elif (network_mode == "genes" or network_mode == "expression"):
      #Filter edges based on gene
      if network_mode == "genes":
          genes = pd.read_table("${genes}", header=None)

      elif network_mode == "expression":
          genes = pd.DataFrame(pd.unique(expression_diff.index.dropna().values.ravel('K')))

      Promoter_Gene_filt_g = Promoter_Gene[Promoter_Gene['Target'].isin(genes.iloc[:,0])]
      if promoter_promoter =="true":
          Promoter_Promoter_filt_g = Promoter_Promoter[Promoter_Promoter['Source'].isin(Promoter_Gene_filt_g['Source']) | Promoter_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
      Distal_Promoter_filt_g = Distal_Promoter[Distal_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
      Factor_Promoter_filt_g = Factor_Promoter[Factor_Promoter['Target'].isin(Promoter_Gene_filt_g['Source'])]
      Factor_Distal_filt_g = Factor_Distal[Factor_Distal['Target'].isin(Distal_Promoter_filt_g['Source'])]

  elif network_mode == "differential":
      #Filter edges based on differnetial peaks
      Distal_Promoter_filt_diff = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal_diff['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter_diff['Target'])]
      Distal_Promoter_filt_up = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal_up['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter_up['Target'])]
      Distal_Promoter_filt_down = Distal_Promoter[Distal_Promoter['Source'].isin(Factor_Distal_down['Target']) | Distal_Promoter['Target'].isin(Factor_Promoter_down['Target'])]

      if promoter_promoter =="true":
          Promoter_Promoter_filt_diff = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter_diff['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter_diff['Target'])]
          Promoter_Gene_filt_diff = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_diff['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_diff['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_diff['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_diff['Target'])]
          Promoter_Promoter_filt_up = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter_up['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter_up['Target'])]
          Promoter_Gene_filt_up = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_up['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_up['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_up['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_up['Target'])]
          Promoter_Promoter_filt_down = Promoter_Promoter[Promoter_Promoter['Source'].isin(Factor_Promoter_down['Target']) | Promoter_Promoter['Target'].isin(Factor_Promoter_down['Target'])]
          Promoter_Gene_filt_down = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_down['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_down['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_down['Source']) | Promoter_Gene['Source'].isin(Promoter_Promoter_filt_down['Target'])]
      else:
          Promoter_Gene_filt_diff = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_diff['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_diff['Source'])]
          Promoter_Gene_filt_up = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_up['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_up['Source'])]
          Promoter_Gene_filt_down = Promoter_Gene[Promoter_Gene['Source'].isin(Factor_Promoter_down['Target']) | Promoter_Gene['Source'].isin(Distal_Promoter_filt_down['Source'])]

  ### Creating edge table for cytoscape
  if network_mode == "all":
    if promoter_promoter =="true":
        Egdes = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Promoter, Promoter_Gene]).drop_duplicates()
    else:
        Egdes = Factor_Distal.append([Factor_Promoter, Distal_Promoter, Promoter_Gene]).drop_duplicates()

  elif network_mode == "factor":
    if promoter_promoter =="true":
        Egdes =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()
    else:
        Egdes =  Factor_Distal.append([Factor_Promoter, Distal_Promoter_filt_f, Promoter_Gene_filt_f]).drop_duplicates()

  elif (network_mode == "genes" or network_mode == "expression"):
    if promoter_promoter =="true":
        Egdes =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()
    else:
        Egdes =  Factor_Distal_filt_g.append([Factor_Promoter_filt_g, Distal_Promoter_filt_g, Promoter_Gene_filt_g]).drop_duplicates()

  elif network_mode == "differential":
    if promoter_promoter =="true":
        Egdes =  Factor_Distal_diff.append([Factor_Promoter_diff, Distal_Promoter_filt_diff, Promoter_Promoter_filt_diff, Promoter_Gene_filt_diff]).drop_duplicates()
        Egdes_up =  Factor_Distal_up.append([Factor_Promoter_up, Distal_Promoter_filt_up, Promoter_Promoter_filt_up, Promoter_Gene_filt_up]).drop_duplicates()
        Egdes_down =  Factor_Distal_down.append([Factor_Promoter_down, Distal_Promoter_filt_down, Promoter_Promoter_filt_down, Promoter_Gene_filt_down]).drop_duplicates()
    else:
        Egdes =  Factor_Distal_diff.append([Factor_Promoter_diff, Distal_Promoter_filt_diff, Promoter_Gene_filt_diff]).drop_duplicates()
        Egdes_up =  Factor_Distal_up.append([Factor_Promoter_up, Distal_Promoter_filt_up, Promoter_Gene_filt_up]).drop_duplicates()
        Egdes_down =  Factor_Distal_down.append([Factor_Promoter_down, Distal_Promoter_filt_down, Promoter_Gene_filt_down]).drop_duplicates()

    Egdes_up.to_csv("Network_Edges_${prefix}_interactions_up.txt", index=False, sep='\t' )
    Egdes_down.to_csv("Network_Edges_${prefix}_interactions_down.txt", index=False, sep='\t' )

  Egdes.to_csv("Network_Edges_${prefix}_interactions.txt", index=False, sep='\t' )

  ### Creating node table for cytoscape

  #Peak padj and log2FC
  Factor_stat_1 = Factor_Interaction.loc[:, ['Anchor1', 'Peak1_ID',  'padj_1', 'log2FC_1']].dropna(subset=['Peak1_ID']).drop_duplicates()
  Factor_stat_1.columns = ['Anchor', 'Peak_ID',  'padj', 'log2FC']
  Factor_stat_2 = Factor_Interaction.loc[:, ['Anchor2', 'Peak2_ID',  'padj_2', 'log2FC_2']].dropna(subset=['Peak2_ID']).drop_duplicates()
  Factor_stat_2.columns = ['Anchor', 'Peak_ID',  'padj', 'log2FC']
  Factor_stat = Factor_stat_1.append(Factor_stat_2).drop_duplicates()

  if network_mode == "all":
    #Specifying node type for all nodes
    Nodes = pd.DataFrame(pd.unique(Egdes[['Source', 'Target']].dropna().values.ravel('K')))
    Nodes.columns=['Node']
    if promoter_promoter =="true":
        Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Distal_Promoter['Target']) | Nodes['Node'].isin(Promoter_Promoter['Source']) | Nodes['Node'].isin(Promoter_Promoter['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))
    else:
        Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter['Source']), 'Distal',
                                     (np.where(Nodes['Node'].isin(Distal_Promoter['Target']) | Nodes['Node'].isin(Promoter_Promoter['Source']), 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene['Target']), 'Gene', np.nan)))))))

  elif network_mode == "factor":
    # Specifying node type for all nodes that are associated with factor binding
    Nodes = pd.DataFrame(pd.unique(Egdes[['Source', 'Target']].dropna().values.ravel('K')))
    Nodes.columns=['Node']
    if promoter_promoter =="true":
        Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_f['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))
    else:
        Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal['Source']) | Nodes['Node'].isin(Factor_Promoter['Source']), 'Factor',
                                (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Source']), 'Distal',
                                   (np.where(Nodes['Node'].isin(Distal_Promoter_filt_f['Target']) | Nodes['Node'].isin(Factor_Promoter['Source']) , 'Promoter',
                                      (np.where(Nodes['Node'].isin(Promoter_Gene_filt_f['Target']), 'Gene', np.nan)))))))

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

  elif network_mode == "differential":
      # Specifying node type for all nodes that are associated with differntial factor binding
      Nodes = pd.DataFrame(pd.unique(Egdes[['Source', 'Target']].dropna().values.ravel('K')))
      Nodes.columns=['Node']
      Nodes_up = pd.DataFrame(pd.unique(Egdes_up[['Source', 'Target']].dropna().values.ravel('K')))
      Nodes_up.columns=['Node']
      Nodes_down = pd.DataFrame(pd.unique(Egdes_down[['Source', 'Target']].dropna().values.ravel('K')))
      Nodes_down.columns=['Node']

      if promoter_promoter =="true":
          Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_diff['Source']) | Nodes['Node'].isin(Factor_Promoter_diff['Source']), 'Factor',
                                    (np.where(Nodes['Node'].isin(Distal_Promoter_filt_diff['Source']), 'Distal',
                                       (np.where(Nodes['Node'].isin(Distal_Promoter_filt_diff['Target']) | Nodes['Node'].isin(Promoter_Promoter_filt_diff['Source']) | Nodes['Node'].isin(Promoter_Promoter_filt_diff['Target']), 'Promoter',
                                          (np.where(Nodes['Node'].isin(Promoter_Gene_filt_diff['Target']), 'Gene', np.nan)))))))

          Nodes_up['Node_type'] = np.where(Nodes_up['Node'].isin(Factor_Distal_up['Source']) | Nodes_up['Node'].isin(Factor_Promoter_up['Source']), 'Factor',
                                      (np.where(Nodes_up['Node'].isin(Distal_Promoter_filt_up['Source']), 'Distal',
                                         (np.where(Nodes_up['Node'].isin(Distal_Promoter_filt_up['Target']) | Nodes_up['Node'].isin(Promoter_Promoter_filt_up['Source']) | Nodes_up['Node'].isin(Promoter_Promoter_filt_up['Target']), 'Promoter',
                                            (np.where(Nodes_up['Node'].isin(Promoter_Gene_filt_up['Target']), 'Gene', np.nan)))))))

          Nodes_down['Node_type'] = np.where(Nodes_down['Node'].isin(Factor_Distal_down['Source']) | Nodes_down['Node'].isin(Factor_Promoter_down['Source']), 'Factor',
                                      (np.where(Nodes_down['Node'].isin(Distal_Promoter_filt_down['Source']), 'Distal',
                                         (np.where(Nodes_down['Node'].isin(Distal_Promoter_filt_down['Target']) | Nodes_down['Node'].isin(Promoter_Promoter_filt_down['Source']) | Nodes_down['Node'].isin(Promoter_Promoter_filt_down['Target']), 'Promoter',
                                            (np.where(Nodes_down['Node'].isin(Promoter_Gene_filt_down['Target']), 'Gene', np.nan)))))))

      else:
          Nodes['Node_type'] = np.where(Nodes['Node'].isin(Factor_Distal_diff['Source']) | Nodes['Node'].isin(Factor_Promoter_diff['Source']), 'Factor',
                                  (np.where(Nodes['Node'].isin(Distal_Promoter_filt_diff['Source']), 'Distal',
                                     (np.where(Nodes['Node'].isin(Distal_Promoter_filt_diff['Target'])  | Nodes['Node'].isin(Factor_Promoter_diff['Target']), 'Promoter',
                                        (np.where(Nodes['Node'].isin(Promoter_Gene_filt_diff['Target']), 'Gene', np.nan)))))))

          Nodes_up['Node_type'] = np.where(Nodes_up['Node'].isin(Factor_Distal_up['Source']) | Nodes_up['Node'].isin(Factor_Promoter_up['Source']), 'Factor',
                                      (np.where(Nodes_up['Node'].isin(Distal_Promoter_filt_up['Source']), 'Distal',
                                         (np.where(Nodes_up['Node'].isin(Distal_Promoter_filt_up['Target'])  | Nodes_up['Node'].isin(Factor_Promoter_up['Target']), 'Promoter',
                                            (np.where(Nodes_up['Node'].isin(Promoter_Gene_filt_up['Target']), 'Gene', np.nan)))))))

          Nodes_down['Node_type'] = np.where(Nodes_down['Node'].isin(Factor_Distal_down['Source']) | Nodes_down['Node'].isin(Factor_Promoter_down['Source']), 'Factor',
                                      (np.where(Nodes_down['Node'].isin(Distal_Promoter_filt_down['Source']), 'Distal',
                                         (np.where(Nodes_down['Node'].isin(Distal_Promoter_filt_down['Target']) | Nodes_down['Node'].isin(Factor_Promoter_down['Target']), 'Promoter',
                                            (np.where(Nodes_down['Node'].isin(Promoter_Gene_filt_down['Target']), 'Gene', np.nan)))))))


      Factor_stat_up = Factor_stat.loc[(Factor_stat['log2FC'] >= ${log2FC}) & (Factor_stat['padj'] <= ${padj}),:]
      Factor_stat_up = Factor_stat_up.sort_values('padj').drop_duplicates(subset=['Anchor'],keep='first').sort_index()
      Nodes_peak_anno_up = Nodes_up[Nodes_up.Node_type.isin(('Distal', 'Promoter'))].merge(Factor_stat_up.loc[:,['Anchor', 'padj', 'log2FC']], how='left', right_on='Anchor', left_on='Node')
      Nodes_gene_anno_up = Nodes_up[Nodes_up.Node_type =='Gene'].merge(expression.loc[:,['padj', 'log2FC']], how='left', right_index=True, left_on='Node')
      Nodes_anno_up = Nodes_peak_anno_up.append(Nodes_gene_anno_up)
      Nodes_up=Nodes_up.merge(Nodes_anno_up[['Node','padj', 'log2FC']], on='Node', how='left')
      Nodes_up.to_csv("Network_Nodes_${prefix}_interactions_up.txt", index=False, sep='\t' )

      Factor_stat_down = Factor_stat.loc[(Factor_stat['log2FC'] <= -${log2FC}) & (Factor_stat['padj'] <= ${padj}),:]
      Factor_stat_down = Factor_stat_down.sort_values('padj').drop_duplicates(subset=['Anchor'],keep='first').sort_index()
      Nodes_peak_anno_down = Nodes_down[Nodes_down.Node_type.isin(('Distal', 'Promoter'))].merge(Factor_stat_down.loc[:,['Anchor', 'padj', 'log2FC']], how='left', right_on='Anchor', left_on='Node')
      Nodes_gene_anno_down = Nodes_down[Nodes_down.Node_type =='Gene'].merge(expression.loc[:,['padj', 'log2FC']], how='left', right_index=True, left_on='Node')
      Nodes_anno_down = Nodes_peak_anno_down.append(Nodes_gene_anno_down)
      Nodes_down=Nodes_down.merge(Nodes_anno_down[['Node','padj', 'log2FC']], on='Node', how='left')
      Nodes_down.to_csv("Network_Nodes_${prefix}_interactions_down.txt", index=False, sep='\t' )

  #Adding node stats
  Factor_stat = Factor_stat.sort_values('padj').drop_duplicates(subset=['Anchor'],keep='first').sort_index()
  Nodes_peak_anno = Nodes[Nodes.Node_type.isin(('Distal', 'Promoter'))].merge(Factor_stat.loc[:,['Anchor', 'padj', 'log2FC']], how='left', right_on='Anchor', left_on='Node')

  skip_expression = "${skip_expression}"
  if skip_expression == "false":
    Nodes_gene_anno = Nodes[Nodes.Node_type =='Gene'].merge(expression.loc[:,['padj', 'log2FC']], how='left', right_index=True, left_on='Node')
    Nodes_anno = Nodes_peak_anno.append(Nodes_gene_anno)
  else:
    Nodes_anno = Nodes_peak_anno

  Nodes=Nodes.merge(Nodes_anno[['Node','padj', 'log2FC']], on='Node', how='left')
  Nodes.to_csv("Network_Nodes_${prefix}_interactions.txt", index=False, sep='\t' )
  """
}

if ({params.network | params.complete} && params.network_mode !="differential" ){
  ch_nodes_up = file("No_nodes_up")
  ch_nodes_down = file("No_nodes_down")
  ch_edges_up = file("No_edges_up")
  ch_edges_down = file("No_edges_down")

}

/*
 * 11. NETWORK VISULAIZATION
 */
process NETWORK_VISUALIZATION {
  publishDir "${params.outdir}/tmp/process11", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Network/Visualization", mode: 'copy'

  when:
  params.network | params.complete

  input:
  path nodes from ch_nodes
  path edges from ch_edges
  val log2FC from Channel.value(params.log2FC)
  val padj from Channel.value(params.padj)
  val expression_log2FC from Channel.value(params.expression_log2FC)
  val expression_padj from Channel.value(params.expression_padj)
  val mode from Channel.value(params.mode)
  val network_mode from Channel.value(params.network_mode)
  val use_peakscore from Channel.value(params.use_peakscore)

  path nodes_up from ch_nodes_up
  path nodes_down from ch_nodes_down
  path edges_up from ch_edges_up
  path edges_down from ch_edges_down

  output:
  path  "Network.pdf" into ch_network_pdf
  path  "Network_up.pdf" optional true into ch_network_up_pdf
  path  "Network_down.pdf" optional true into ch_network_down_pdf

  path "Network.xgmml" optional true into ch_network_xgmml
  path "Network_up.xgmml" optional true into ch_network_up_xgmml
  path "Network_down.xgmml" optional true into ch_network_down_xgmml


  script:
  """
  #!/usr/local/bin/Rscript --vanilla

  require(RCy3)
  require(circlize)
  require(viridisLite)

  #Loading and organizing data
  mode="${mode}"
  network_mode="${network_mode}"

  if (mode=="differential"){
    nodes <- read.table("${nodes}", header=TRUE, sep="\t", col.names=c("id", "type", "padj", "log2FC"))
    edges <- read.table("${edges}", header=TRUE, sep="\t", col.names=c("source", "target", "score", "type"))
  } else{
    nodes <- read.table("${nodes}", header=TRUE, sep="\t", col.names=c("id", "type"))
    edges <- read.table("${edges}", header=TRUE, sep="\t", col.names=c("source", "target","score", "type"))
  }
  edges[,"interaction"] <- "interacts"
  edges[,"name"] <- paste(edges\$source, "(interacts)", edges\$target, sep=" ")
  createNetworkFromDataFrames(nodes,edges, title="Network", collection="Networks" )

  #Creating a style for the network
  style.name = "style1"
  defaults <- list(NETWORK_BACKGROUND_PAINT="#ffffff",
                   NODE_BORDER_PAINT="#000000",
                   NODE_BORDER_WIDTH=1,
                   EDGE_WIDTH=1)
  nodeShape <- mapVisualProperty('Node Shape','type','d',c("Factor","Distal", "Promoter", "Gene"), c("ELLIPSE","ROUND_RECTANGLE","ROUND_RECTANGLE", "ELLIPSE"))
  nodeFills <- mapVisualProperty('Node Fill Color','type','d',c("Factor","Distal", "Promoter", "Gene"), c("#ffffff   ","#ededed","#8e8e8e", "#3e3e3e"))
  nodeHeights <- mapVisualProperty('Node Height', 'type', 'd', c("Factor","Distal", "Promoter", "Gene"), c(120, 16, 16, 50))
  nodeWidths <- mapVisualProperty('Node Width', 'type', 'd', c("Factor","Distal", "Promoter", "Gene"), c(120, 120, 120, 50))
  nodeBorderWidths <- mapVisualProperty('Node Border Width', 'type', 'd', c("Factor","Distal", "Promoter", "Gene"), c(6, 1, 1, 1))
  nodeZ <- mapVisualProperty('Node Z Location', 'type', 'd', c("Factor","Distal", "Promoter", "Gene"), c(4, 2, 1, 3))
  nodeLabels <- mapVisualProperty('Node Label','name','p')
  nodeLabelSize <- mapVisualProperty('Node Label Font Size', 'type', 'd', c("Factor","Distal", "Promoter", "Gene"), c(32, 8, 8, 10))
  nodeLabelColor <- mapVisualProperty('Node Label Color', 'type', 'd',c("Factor","Distal", "Promoter", "Gene"), c("#000000", "#000000", "#000000", "#ffffff"))
  createVisualStyle(style.name, defaults, list(nodeShape, nodeFills, nodeHeights, nodeWidths, nodeLabels,nodeLabelSize, nodeLabelColor, nodeBorderWidths, nodeZ))
  setVisualStyle(style.name)
  lockNodeDimensions(FALSE, style.name)

  #Set edge width based on q-value for interactions (and peak score for peaks if the argument use_peakscore is set to true)
  use_peakscore = "${use_peakscore}"
  if (use_peakscore=="true"){
    edges_interaction <- edges[edges\$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
    edges_interaction\$score <- 9*((edges_interaction\$score-min(edges_interaction\$score))/(max(edges_interaction\$score)-min(edges_interaction\$score)))+1
    edges_factor <- edges[edges\$type %in% c("Factor-Distal", "Factor-Promoter"),]
    edges_factor\$score <- 9*((edges_factor\$score-min(edges_factor\$score))/(max(edges_factor\$score)-min(edges_factor\$score)))+1
    edges_score <- rbind(edges_interaction, edges_factor)
  } else{
    edges_score <- edges[edges\$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
    edges_score\$score <- 9*((edges_score\$score-min(edges_score\$score))/(max(edges_score\$score)-min(edges_score\$score)))+1
  }
  setEdgePropertyBypass(edge.names=edges_score\$name, new.values=edges_score\$score, visual.property='EDGE_WIDTH',bypass = TRUE)


  #For differntial mode:color by log2FC of peaks and GENES
  if (mode=="differential"){
    map2color<-function(x,pal,limits=NULL){
      if(is.null(limits)) limits=range(x)
      pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
    }
    pal_up <- viridisLite::rocket(100, begin = 0.5, end = 0.25, direction = 1)
    pal_down <- viridisLite::mako(100, begin = 0.6, end = 0.35, direction = -1)
    nodes_up_peak <- nodes[nodes\$log2FC >=${log2FC} & !is.na(nodes\$log2FC) & nodes\$padj <=${padj} & !is.na(nodes\$padj) & nodes\$type %in% c("Distal", "Promoter"),]
    nodes_down_peak <- nodes[nodes\$log2FC <= -${log2FC} & !is.na(nodes\$log2FC) & nodes\$padj <=${padj} & !is.na(nodes\$padj) & nodes\$type %in% c("Distal", "Promoter"),]
    nodes_up_gene <- nodes[nodes\$log2FC >=${expression_log2FC} & !is.na(nodes\$log2FC) & nodes\$padj <=${expression_padj} & !is.na(nodes\$padj) & nodes\$type == "Gene",]
    nodes_down_gene <- nodes[nodes\$log2FC <= -${expression_log2FC} & !is.na(nodes\$log2FC) & nodes\$padj <=${expression_padj} & !is.na(nodes\$padj) & nodes\$type == "Gene",]
    nodes_up <- rbind(nodes_up_peak, nodes_up_gene)
    nodes_down <- rbind(nodes_down_peak, nodes_down_gene)
    if (nrow(nodes_up)>0){
      nodes_up[,"color"] <- map2color(log(nodes_up[,"log2FC"]),pal_up)
    }
    if (nrow(nodes_down)>0){
      nodes_down[,"color"] <- map2color(log(abs(nodes_down[,"log2FC"])),pal_down)
    }
    nodes_diff <- rbind(nodes_up, nodes_down)
    if (ncol(nodes_diff)==5){
      nodes_diff[,"color"] <- gsub("FF\$", "", nodes_diff[,"color"])
      setNodePropertyBypass(nodes_diff\$id,nodes_diff\$color,'NODE_FILL_COLOR',bypass = TRUE)
    }
  }
  toggleGraphicsDetails()
  exportImage("Network.pdf", 'PDF')
  exportNetwork("Network.xgmml", type= 'xGMML')

  if (mode=="differential" & network_mode=="differential"){
  #Up
    up_nodes <- read.table("${nodes_up}", header=TRUE, sep="\t", col.names=c("id", "type", "padj", "log2FC"))
    up_edges <- read.table("${edges_up}", header=TRUE, sep="\t", col.names=c("source", "target", "score", "type"))
    up_edges[,"interaction"] <- "interacts"
    up_edges[,"name"] <- paste(up_edges\$source, "(interacts)", up_edges\$target, sep=" ")
    createNetworkFromDataFrames(up_nodes,up_edges, title="Network_up", collection="Networks" )
    setVisualStyle(style.name)
    lockNodeDimensions(FALSE, style.name)

    use_peakscore = "${use_peakscore}"
    if (use_peakscore=="true"){
      up_edges_interaction <- up_edges[up_edges\$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
      up_edges_interaction\$score <- 9*((up_edges_interaction\$score-min(up_edges_interaction\$score))/(max(up_edges_interaction\$score)-min(up_edges_interaction\$score)))+1
      up_edges_factor <- up_edges[up_edges\$type %in% c("Factor-Distal", "Factor-Promoter"),]
      up_edges_factor\$score <- 9*((up_edges_factor\$score-min(up_edges_factor\$score))/(max(up_edges_factor\$score)-min(up_edges_factor\$score)))+1
      up_edges_score <- rbind(up_edges_interaction, up_edges_factor)
    } else{
      up_edges_score <- up_edges[up_edges\$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
      up_edges_score\$score <- 9*((up_edges_score\$score-min(up_edges_score\$score))/(max(up_edges_score\$score)-min(up_edges_score\$score)))+1
    }
    setEdgePropertyBypass(edge.names=up_edges_score\$name, new.values=up_edges_score\$score, visual.property='EDGE_WIDTH',bypass = TRUE)

    up_nodes_up_peak <- up_nodes[up_nodes\$log2FC >=${log2FC} & !is.na(up_nodes\$log2FC) & up_nodes\$padj <=${padj} & !is.na(up_nodes\$padj) & up_nodes\$type %in% c("Distal", "Promoter"),]
    up_nodes_down_peak <- up_nodes[up_nodes\$log2FC <= -${log2FC} & !is.na(up_nodes\$log2FC) & up_nodes\$padj <=${padj} & !is.na(up_nodes\$padj) & up_nodes\$type %in% c("Distal", "Promoter"),]
    up_nodes_up_gene <- up_nodes[up_nodes\$log2FC >=${expression_log2FC} & !is.na(up_nodes\$log2FC) & up_nodes\$padj <=${expression_padj} & !is.na(up_nodes\$padj) & up_nodes\$type == "Gene",]
    up_nodes_down_gene <- up_nodes[up_nodes\$log2FC <= -${expression_log2FC} & !is.na(up_nodes\$log2FC) & up_nodes\$padj <=${expression_padj} & !is.na(up_nodes\$padj) & up_nodes\$type == "Gene",]
    up_nodes_up <- rbind(up_nodes_up_peak, up_nodes_up_gene)
    up_nodes_down <- rbind(up_nodes_down_peak, up_nodes_down_gene)
    if (nrow(up_nodes_up)>0){
      up_nodes_up[,"color"] <- map2color(log(up_nodes_up[,"log2FC"]),pal_up)
    }
    if (nrow(up_nodes_down)>0){
      up_nodes_down[,"color"] <- map2color(log(abs(up_nodes_down[,"log2FC"])),pal_down)
    }
    up_nodes_diff <- rbind(up_nodes_up, up_nodes_down)
    if (ncol(up_nodes_diff)==5){
      up_nodes_diff[,"color"] <- gsub("FF\$", "", up_nodes_diff[,"color"])
      setNodePropertyBypass(up_nodes_diff\$id,up_nodes_diff\$color,'NODE_FILL_COLOR',bypass = TRUE)
    }
  toggleGraphicsDetails()
  exportImage("Network_up.pdf", 'PDF')
  exportNetwork("Network_up.xgmml", type= 'xGMML')

  #Down
    down_nodes <- read.table("${nodes_down}", header=TRUE, sep="\t", col.names=c("id", "type", "padj", "log2FC"))
    down_edges <- read.table("${edges_down}", header=TRUE, sep="\t", col.names=c("source", "target", "score", "type"))
    down_edges[,"interaction"] <- "interacts"
    down_edges[,"name"] <- paste(down_edges\$source, "(interacts)", down_edges\$target, sep=" ")
    createNetworkFromDataFrames(down_nodes,down_edges, title="Network_up", collection="Networks" )
    setVisualStyle(style.name)
    lockNodeDimensions(FALSE, style.name)

    use_peakscore = "${use_peakscore}"
    if (use_peakscore=="true"){
      down_edges_interaction <- down_edges[down_edges\$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
      down_edges_interaction\$score <- 9*((down_edges_interaction\$score-min(down_edges_interaction\$score))/(max(down_edges_interaction\$score)-min(down_edges_interaction\$score)))+1
      down_edges_factor <- down_edges[down_edges\$type %in% c("Factor-Distal", "Factor-Promoter"),]
      down_edges_factor\$score <- 9*((down_edges_factor\$score-min(down_edges_factor\$score))/(max(down_edges_factor\$score)-min(down_edges_factor\$score)))+1
      down_edges_score <- rbind(down_edges_interaction, down_edges_factor)
    } else{
      down_edges_score <- down_edges[down_edges\$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
      down_edges_score\$score <- 9*((down_edges_score\$score-min(down_edges_score\$score))/(max(down_edges_score\$score)-min(down_edges_score\$score)))+1
    }
    setEdgePropertyBypass(edge.names=down_edges_score\$name, new.values=down_edges_score\$score, visual.property='EDGE_WIDTH',bypass = TRUE)

    down_nodes_up_peak <- down_nodes[down_nodes\$log2FC >=${log2FC} & !is.na(down_nodes\$log2FC) & down_nodes\$padj <=${padj} & !is.na(down_nodes\$padj) & down_nodes\$type %in% c("Distal", "Promoter"),]
    down_nodes_down_peak <- down_nodes[down_nodes\$log2FC <= -${log2FC} & !is.na(down_nodes\$log2FC) & down_nodes\$padj <=${padj} & !is.na(down_nodes\$padj) & down_nodes\$type %in% c("Distal", "Promoter"),]
    down_nodes_up_gene <- down_nodes[down_nodes\$log2FC >=${expression_log2FC} & !is.na(down_nodes\$log2FC) & down_nodes\$padj <=${expression_padj} & !is.na(down_nodes\$padj) & down_nodes\$type == "Gene",]
    down_nodes_down_gene <- down_nodes[down_nodes\$log2FC <= -${expression_log2FC} & !is.na(down_nodes\$log2FC) & down_nodes\$padj <=${expression_padj} & !is.na(down_nodes\$padj) & down_nodes\$type == "Gene",]
    down_nodes_up <- rbind(down_nodes_up_peak, down_nodes_up_gene)
    down_nodes_down <- rbind(down_nodes_down_peak, down_nodes_down_gene)
    if (nrow(down_nodes_up)>0){
      down_nodes_up[,"color"] <- map2color(log(down_nodes_up[,"log2FC"]),pal_up)
    }
    if (nrow(down_nodes_down)>0){
      down_nodes_down[,"color"] <- map2color(log(abs(down_nodes_down[,"log2FC"])),pal_down)
    }
    down_nodes_diff <- rbind(down_nodes_up, down_nodes_down)
    if (ncol(down_nodes_diff)==5){
      down_nodes_diff[,"color"] <- gsub("FF\$", "", down_nodes_diff[,"color"])
      setNodePropertyBypass(down_nodes_diff\$id,down_nodes_diff\$color,'NODE_FILL_COLOR',bypass = TRUE)
    }
  toggleGraphicsDetails()
  exportImage("Network_down.pdf", 'PDF')
  exportNetwork("Network_down.xgmml", type= 'xGMML')
  }
  """
}

if ({params.upset_plot | params.circos_plot | params.complete} && !params.filter_genes ){
  ch_upset_promoter_genes = file("No_input_upset_promoter_genes")
  ch_upset_distal_genes = file("No_input_upset_distal_genes")
}

/*
 * 12. UPSET PLOT FOR FACTOR BINDING IN ANCHOR POINTS - ALL INTERACTIONS WITH FACTOR AND INTERACTIONS FILTERES FOR GENES
 */
process UPSET_PLOT {
  publishDir "${params.outdir}/tmp/process12", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Co_occupancy/Upset", mode: 'copy', pattern: 'Upset_plot_*.pdf', enabled: params.upset_plot | params.complete

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
  path "Circos_peaks_${prefix}_interactions.txt" optional true into ch_circos_f
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
    circos_f.to_csv("Circos_peaks_${prefix}_interactions.txt", index=False, sep='\t' )

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
 * 13. CIRCOS PLOTS - ALL INTERACTIONS WITH FACTOR AND INTERACTIONS FILTERED FOR GENES
 */
process CIRCOS_PLOT {
  publishDir "${params.outdir}/tmp/process13", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Co_occupancy/Circos", mode: 'copy', enabled: params.circos_plot | params.complete

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
  pdf("Circos_plot_peaks.pdf")
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


/*
 * 14. DIFFERTIAL PEAKS ASSOCIATES WITH DIFFERENTIAL GENE EXPRESSION
 */
process DIFFERENTIAL_EXPRESSION_ASSOCIATED_PEAKS {
  publishDir "${params.outdir}/tmp/process14", mode: 'copy', enabled: params.save_tmp
  publishDir "${params.outdir}/Differential_expression_associated_peaks", mode: 'copy'

  when:
  params.mode == 'differential' && !params.skip_expression

  input:
  set val(sample), file(annotated_peaks) from ch_peak_PLACseq_annotated_for_differntial_expression
  path expression from ch_expression_2
  val prefix from Channel.value(params.prefix)
  val log2FC_column from Channel.value(params.log2FC_column)
  val padj_column from Channel.value(params.padj_column)
  val log2FC from Channel.value(params.log2FC)
  val padj from Channel.value(params.padj)
  val expression_padj from Channel.value(params.expression_padj)
  val expression_log2FC from Channel.value(params.expression_log2FC)
  val expression_padj_column from Channel.value(params.expression_padj_column)
  val expression_log2FC_column from Channel.value(params.expression_log2FC_column)


  output:
  path "${sample}_${prefix}_annotated_differential_expression.txt" into ch_peak_PLACseq_annotated_differntial_expression
  path "${sample}_${prefix}_annotated_differential_expression_proximal_activating.txt" into ch_peak_PLACseq_annotated_differntial_expression_proximal_activating
  path "${sample}_${prefix}_annotated_differential_expression_proximal_repressive.txt" into ch_peak_PLACseq_annotated_differntial_expression_proximal_repressive
  path "${sample}_${prefix}_annotated_differential_expression_distal_activating.txt" into ch_peak_PLACseq_annotated_differntial_expression_distal_activating
  path "${sample}_${prefix}_annotated_differential_expression_distal_repressive.txt" into ch_peak_PLACseq_annotated_differntial_expression_distal_repressive
  script:
  """
  #!/usr/bin/env python

  import pandas as pd
  import numpy as np

  #Loading data
  annotated_peaks = pd.read_table("${annotated_peaks}", index_col=0).sort_index()
  expression = pd.read_table("${expression}", index_col=0).sort_index()
  expression = expression.iloc[:, [${expression_log2FC_column}-2, ${expression_padj_column}-2]]
  expression.columns = ['Gene_log2FC', 'Gene_padj']

  #Adding gene stats to annotated peaks
  annotated_peaks_expression = annotated_peaks.reset_index().merge(expression, how='left', left_on='Gene', right_on='symbol').set_index('Peak')

  #Finding differntial peaks associated with changes in expression (sepeartion based on.annotation proximal/distal and activating/prepressive)
  annotated_peaks_expression_proximal_activating = annotated_peaks_expression[(annotated_peaks_expression.Annotation.isin(['Promoter', 'Proximal_anno'])) & (annotated_peaks_expression.padj < ${padj}) & (annotated_peaks_expression.Gene_padj < ${expression_padj}) & (((annotated_peaks_expression.log2FC > ${log2FC}) & (annotated_peaks_expression.Gene_log2FC > ${expression_log2FC})) | ((annotated_peaks_expression.log2FC < -${log2FC}) & (annotated_peaks_expression.Gene_log2FC < -${expression_log2FC})))]
  annotated_peaks_expression_proximal_repressive = annotated_peaks_expression[(annotated_peaks_expression.Annotation.isin(['Promoter', 'Proximal_anno'])) & (annotated_peaks_expression.padj < ${padj}) & (annotated_peaks_expression.Gene_padj < ${expression_padj}) & (((annotated_peaks_expression.log2FC < -${log2FC}) & (annotated_peaks_expression.Gene_log2FC > ${expression_log2FC})) | ((annotated_peaks_expression.log2FC > ${log2FC}) & (annotated_peaks_expression.Gene_log2FC < -${expression_log2FC})))]
  annotated_peaks_expression_distal_activating = annotated_peaks_expression[(annotated_peaks_expression.Annotation == 'Interaction_anno') & (annotated_peaks_expression.padj < ${padj}) & (annotated_peaks_expression.Gene_padj < ${expression_padj}) & (((annotated_peaks_expression.log2FC >${log2FC}) & (annotated_peaks_expression.Gene_log2FC > ${expression_log2FC})) | ((annotated_peaks_expression.log2FC < -${log2FC}) & (annotated_peaks_expression.Gene_log2FC < -${expression_log2FC})))]
  annotated_peaks_expression_distal_repressive = annotated_peaks_expression[(annotated_peaks_expression.Annotation =='Interaction_anno') & (annotated_peaks_expression.padj < ${padj}) & (annotated_peaks_expression.Gene_padj < ${expression_padj}) & (((annotated_peaks_expression.log2FC < -${log2FC}) & (annotated_peaks_expression.Gene_log2FC > ${expression_log2FC})) | ((annotated_peaks_expression.log2FC > ${log2FC}) & (annotated_peaks_expression.Gene_log2FC < -${expression_log2FC})))]

  #Adding activatin/repressive funktion to all annotated peaks
  annotated_peaks_expression['Activating_or_repressive'] = np.where((annotated_peaks_expression.index.isin(annotated_peaks_expression_proximal_activating.index)) | (annotated_peaks_expression.index.isin(annotated_peaks_expression_distal_activating.index)), 'Activating', np.where((annotated_peaks_expression.index.isin(annotated_peaks_expression_proximal_repressive.index)) | (annotated_peaks_expression.index.isin(annotated_peaks_expression_distal_repressive.index)), 'Repressive', ''))

  #Save files
  annotated_peaks_expression.to_csv("${sample}_${prefix}_annotated_differential_expression.txt", index=True, sep='\t' )
  annotated_peaks_expression_proximal_activating.to_csv("${sample}_${prefix}_annotated_differential_expression_proximal_activating.txt", index=True, sep='\t' )
  annotated_peaks_expression_proximal_repressive.to_csv("${sample}_${prefix}_annotated_differential_expression_proximal_repressive.txt", index=True, sep='\t' )
  annotated_peaks_expression_distal_activating.to_csv("${sample}_${prefix}_annotated_differential_expression_distal_activating.txt", index=True, sep='\t' )
  annotated_peaks_expression_distal_repressive.to_csv("${sample}_${prefix}_annotated_differential_expression_distal_repressive.txt", index=True, sep='\t' )
  """
}
