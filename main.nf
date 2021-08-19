params{
  // General parameters
  bed2D = false
  peaks = false
  genome = false
  genes = "No genes specified"
  bed2D_anno = false
  mode = "basic" //basic, multiple, differntial
  outdir = "$baseDir/result"
  homerdir = "$baseDir/homer"
  prefix = "PLACseq"
  proximity_unannotated = false
  multiple_anno = "concentrate"  //keep, concentrate, qvalue
  skip_anno = false
  annotate_interactions = false
  network = false
  network_mode = "factor" //Options: all (for all plac-seq interactions), factor (default, for all interactions have factor binding in at least on anchor point) or genes (show only interactions realted to specified gene list)
  complete = false
  save_tmp = false
  help = false

  // Multiple mode specific
  upset_plot = false
  circos_plot = false
  filter_genes = false

  // Differntial mode specific
  peak_differential="No differntial peak file specified"
  log2FC_column=3
  padj_column=9
  log2FC=1.5
  padj=0.05
  expression=false
}

//TODO:
//General:
  //-Plot network
  //-allow differnt bin size

// MULIPLE:
  //-circos details

//Differntial

//Github:
  //- ??? in textfile
  //- name
  //- Figures
  //- profile (picture, licence etc)
  //-output
