# PLACseq_anno
**Nextflow pipeline for PLAC-seq based peak annotation**

## Introduction
Long range chromatin interaction imortant for gene regualtion (source)???. This pipeline aim to imrpove the peak annotaiton (e.g. ChIP-seq, Cut&Run, ATAC-seq data) compared to tranditional proximity based annotation methods, by taking advantage of known genomic interaction from Proximity Ligation-Assisted ChIP-seq (PLAC-seq) data. Even though the pipeline is designed for PLAC-seq data, it is also possible to use it with data gerneated by other methods that capture long-ange chromatin interactions (e.g. Hi-C).  

More about the idea, how it works and figure???

## Annotation modes
To fashiliate anootation of differnt types of input data, the pipeline can be run in three differnt modes: basic, multiple and differential.

### Basic mode

### Multiple mode

### DIfferntial mode

## Pipeline summary

The pipeline consist of the following processes:

1. BED2D_SPLIT - 2D-bed file is splitted for annotation of the two anchor points separatly
2. ANNOTATE_INTERACTION - Annotation of 2D-bed anchor points using [`HOMER`](http://homer.ucsd.edu/homer/)
3. JOIN_ANNOTATED_INTERACTIONS - Recreating the original 2D-bed interactions with annotations
4. ANNOTATE_PEAKS - Annotation of provided peak file(s) using [`HOMER`](http://homer.ucsd.edu/homer/)
5. SPLIT_ANNOTATED_INTERACTIONS - Annotated 2D-bed file splitted for peak intersection
6. PEAK_INTERACTION_INTERSECT - Peak centered intersection of provided peak file(s) with genomic interactions
7. PEAK_INTERACTION_BASED_ANNOTATION - Performing PLAC-seq based annotation of provided peak file(s)
8. [OPTIONAL] INTERACTION_PEAK_INTERSECT - Interaction centered intersection of provided peak file(s) with genomic interactions
9. [OPTIONAL] ANNOTATE_INTERACTION_WITH_PEAKS - Performing interaction based annotation of 2D-bed file with provided peak files(s)
10. [OPTIONAL] NETWORK - Creating files for network visualization of peak annotation in [`Cytoscape`](https://cytoscape.org/)
11. [OPTIONAL] UPSET_PLOT - Creating Upset plots for overlap of peak files in promoter and distal regions (only available in Multiple mode)
12. [OPTIONAL] CIRCOS PLOT - Creatig Circos plot representing peak overlap in genomic interactions (only available in Multiple mode)


## Run the pipeline

Download the plac_anno_env.yml file and create a conda environemnt that contain all packages neccisary to run the pieline.
```bash
conda create -f plac_anno_env.yaml
```

Dowload the pipeline (including main.nf & nextflow.config). To avoid having specify the path to the config file, make sure to place the two files in the same directory.
```bash
nextflow run PLAC_anno_2.nf --bed2D interactions.bed  --genome mm10 --peaks peaks.txt
```


