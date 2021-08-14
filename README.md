# PLACseq_anno
**Nextflow pipeline for PLAC-seq based peak annotation**

## Introduction
Long range chromatin interaction imortant for gene regualtion (source)???. This pipeline aim to imrpove the peak annotaiton (e.g. ChIP-seq, Cut&Run, ATAC-seq data) compared to tranditional proximity based annotation methods, by taking advantage of known genomic interaction from Proximity Ligation-Assisted ChIP-seq (PLAC-seq) data. Even though the pipeline is designed for PLAC-seq data, it is also possible to use it with data gerneated by other methods that capture long-ange chromatin interactions (e.g. Hi-C).  



## Annotation modes
To fashiliate anootation of differnt types of input data, the pipeline can be run in three differnt modes: basic, multiple and differential.

### Basic mode

### Multiple mode

### DIfferntial mode

## Pipeline summary

The pipeline consist of the following processes:

1. BED2D_SPLIT
2. ANNOTATE_INTERACTION
3. MERGE_ANNOTATED_INTERACTIONS
4. ANNOTATE_PEAKS
5. SPLIT_ANNOTATED_INTERACTIONS
6. PEAK_INTERACTION_INTERSECT
7. PEAK_ANNOTATION
8. INTERACTION_PEAK_INTERSECT
9. ANNOTATE_INTERACTION_WITH_PEAKS
10. NETWORK
11. UPSET_PLOT
12. CIRCOS PLOT


## Run the pipeline

Download the plac_anno_env.yml file and create a conda environemnt that contain all packages neccisary to run the pieline.
```bash
conda create -f plac_anno_env.yaml
```

Dowload the pipeline (including main.nf & nextflow.config). To avoid having specify the path to the config file, make sure to place the two files in the same directory.
```bash
nextflow run PLAC_anno_2.nf --bed2D interactions.bed  --genome mm10 --peaks peaks.txt
```


