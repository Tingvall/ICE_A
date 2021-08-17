# PLACseq_anno
**Nextflow pipeline for PLAC-seq based peak annotation**

## Introduction
Long range chromatin interaction imortant for gene regualtion (source)???. This pipeline aim to imrpove the peak annotaiton (e.g. ChIP-seq, Cut&Run, ATAC-seq data) compared to tranditional proximity based annotation methods, by taking advantage of known genomic interaction from Proximity Ligation-Assisted ChIP-seq (PLAC-seq) data. Even though the pipeline is designed for PLAC-seq data, it is also possible to use it with data gerneated by other methods that capture long-ange chromatin interactions (e.g. Hi-C).  

More about the idea, how it works and figure???

### Annotation modes
To fashiliate anootation of differnt types of input data, the pipeline can be run in three differnt modes: basic, multiple and differential.

#### Basic mode
The basic annotation mode performs PLAC-seq based annotation of a single peak file. Examples include annoation of ATAC-seq peaks or a single  factor from ChIP-seq/Cut&Run. The main output is a text file that for each peak in the input bed file, provides PLAC_seq based annotation.    A genelist, contianing all genes that peaks are annoated is also provided. In addition to the peak-cenetered annotation, it is a slo          possible to perform interaction-centered annotation and create input files for netowrk visualization in Cytoscape.

#### Multiple mode
The multiple annotation mode is designed simmultanous annotation of multiple peak files. This mode is usuful in situtations where multiple peak files exist and where overlap between these peaks are of interest (e.g. CHIP-seq/Cut&Run data for multiple factors in the same cell types). In addition to peak-centered PLAC-seq based annotation for each peak file (identical to basic mode), the multiple mode also provide the option to investigate peak overlap in the form of Netowrk visualization, Upset plot and Circos plot. 

#### Differntial mode
Add???

## Pipeline summary

The pipeline consist of the following processes:

1. **BED2D_SPLIT** - 2D-bed file is splitted for annotation of the two anchor points separatly
2. **ANNOTATE_INTERACTION** - Annotation of 2D-bed anchor points using [`HOMER`](http://homer.ucsd.edu/homer/)
3. **JOIN_ANNOTATED_INTERACTIONS** - Recreating the original 2D-bed interactions with annotations
4. **ANNOTATE_PEAKS** - Annotation of provided peak file(s) using [`HOMER`](http://homer.ucsd.edu/homer/)
5. **SPLIT_ANNOTATED_INTERACTIONS** - Annotated 2D-bed file splitted for peak intersection
6. **PEAK_INTERACTION_INTERSECT** - Peak-centered intersection of provided peak file(s) with genomic interactions
7. **PEAK_INTERACTION_BASED_ANNOTATION** - Performing PLAC-seq based annotation of provided peak file(s)
8. **INTERACTION_PEAK_INTERSECT [Optional]** - Interaction-centered intersection of peak file(s) with genomic interactions
9. **ANNOTATE_INTERACTION_WITH_PEAKS [Optional]** - Interaction based annotation of 2D-bed file with provided peak files(s)
10. **NETWORK [Optional]** - Creating files for network visualization of peak annotation in [`Cytoscape`](https://cytoscape.org/)
11. **UPSET_PLOT [Optional]** - Upset plots for overlap of peak files in promoter and distal regions (only available in Multiple mode)
12. **CIRCOS PLOT [Optional]** - Circos plot representing peak overlap in genomic interactions (only available in Multiple mode)


## Installation

**Create conda environment:**
Download the plac_anno_env.yml file and create a conda environemnt that contain all packages neccisary to run the pieline:
```bash
conda create -f plac_anno_env.yaml
```
Activate conda environment:
```bash
conda activate plac_anno_env
```
**Configure HOMER:**
Before you use the pipeline, HOMER must be configured correctly. For each genome you intend to use, download the reference fasta file by using the following command (for available genomes visit: http://homer.ucsd.edu/homer):
```bash
perl <path to conda>/envs/plac_anno_env/share/homer*/configureHomer.pl -install <genome (e.g. mm10, hg19)>
```
Next, download the promoter annotations for the organism you are working with (for available genomes visit: http://homer.ucsd.edu/homer):
```bash
perl <path to conda>/envs/plac_anno_env/share/homer*/configureHomer.pl -install <organism-p (e.g. mouse-p, human-p)>
```

**Launce pipeline:**
Dowload the pipeline (including main.nf & nextflow.config). To avoid having specify the path to the config file, make sure to place the two files in the same directory. Try to launce the pipeline:
```bash
nextflow run PLAC_anno_2.nf --help
```

## Inputs

#### Required input
| Input | Description |
| --- | --- |
| `--peaks` | Text file specifying the name and path of the peak file(s), see example file: ???|
| `--bed2D` | Path to chromain interaction from PLAC-seq (or any similar method that captures long-range interactions) in 2D-bed format. Currently the pipeline is designed for 2D-bed files created by ???. This input can be replcaed by `bed2D_anno` if an alredy annotated 2D-bed file is available and the arguement `skip_anno` is used. |
| `--genome` | Specification of genome for annotation (e.g. mm10). Currently the annotation is performed by [`HOMER`](http://homer.ucsd.edu/homer/), visit documentation for details and available genomes: http://homer.ucsd.edu/homer |

#### Optional input
| Input | Description |
| --- | --- |
| `--genes` | Only used when the option `--filtering_genes` is specified or if `--network_mode` is set to `genes`. Textfile with genenames (specify???), that is used for filtering of interactions associated with the specified genes. The filitering is performed during plotting of Upset and Circos plot (if `--filtering_genes` is specified) and for network visulaization (if `--network_mode` is set to `genes`). |
| `--bed2D_anno` | Specifies path to annotated 2D-bed file if `--skip_anno` is used. |

## Running the pipeline

### General 
The typical command for running the pipeline is as follows:

```bash
nextflow run PLAC_anno_2.nf --bed2D interactions.bed  --genome mm10 --peaks peaks.txt
```
The default mode is basic, to run the pipeline in another mode specify it with the argumnet `mode`.

#### Arguments
| Argumnet | Description |
| --- | --- |
| `--mode` | Define which mode to run the pipeline in. The options are basic (default), multiple or differntial. |
| `--outdir` | Specifies the output driectory (default: ./results). |
| `--prefix` | Prefix used for interactions (default: PLACseq).|
| `--proximity_unannotated` | Specifies if unannotated distal peaks should be annotated by proximity annotation (default: false) |
| `--multiple_anno` | Defines how to handle peaks annotated to more than one promoter. Options are keep (all anotations are kept with one row for each annotation), concetrate (the annotated peak file is concetrated to only incude one row per peak but information about all annotations are kept) and qvalue (only one annotation per peak is kept. The annotation is decided by the interaction with the lowset qvalue). Default is: concentrate.|
| `--skip_anno` | If you already have an annotated 2D-bed file from a previous run, you can skip the HOMER annotation of the interactions by using this argumnet. Requires specification of path to annotated 2D-bed by using the argumnet `--bed2D_anno`. |
| `--annotate_interactions` | Specifes if interaction-centered annotation with peak overlap should be performed. Only valid if `--complete` is set to false. |
| `--network` | Specifes if files for network visualization in Cytoskape should be created. Only valid if `--complete` is set to false. |
| `--network_mode` | Defines mode network. Options are all (all interaction in the 2D-bed file), factor (all interaction with at least on peak overlap either anchor point) or genes (interactions associates with a genelist, provided by `--genes`). |
| `--complete` | If set to true (default), all available processes for the selected mode and provided inputs are run.|
| `--save_tmp` | If used, all intermediate files are saved in the directory ./tmp. Can be useful for investigating promblems. Default: false.
| `--help` | Help message is shown. |

### Multiple mode specific arguments
When the pipleine is run in multiple mode, some aditional processes based on peak overlap are available. The specific parameters for these processes are listed below:

#### Arguments
| Argumnet | Description |
| --- | --- |
| `--upset_plot` | If specified, Upset plot of peak overlap will be created. Only valid if `--complete` is set to false. |
| `--circos_plot` | If specified, Circos plot of peak overlap will be created. Only valid if `--complete` is set to false. |
| `--filter_genes` | Specifies if additional plot (Upset and/or Circos plots) should be created based on interactions filtered by provided genelist (default: false). This option requires that a genelist is provided with the argument `--genes`. |

### Differntial mode
Add text here???

#### Arguments
| Argumnet | Description |
| --- | --- |
| `` | |

## Outputs and interpretation
