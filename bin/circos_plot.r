#!/usr/bin/env Rscript

### PROCESS 13 MULTIPLE: CIRCOS PLOTS - ALL INTERACTIONS WITH FACTOR AND INTERACTIONS FILTERED FOR GENES ###

# Load libraries
require(circlize)
require(viridis)
require(stringr)
require(mgsub)
require(optparse)

# Arguments
option_list <- list(make_option(c("--circos_f"), type="character", default=NULL, help="Input circos plot data for all interaction with factor overlap.", metavar="path"),
                    make_option(c("--circos_g"), type="character", default=NULL, help="Input circos plot data filtered for for genes.", metavar="path"),
                    make_option(c("--filter_genes"), type="character", default='NULL', help="", metavar="string"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Prepare and plot circos plot
## Regions with factor in at least one anchor point
circos_data_all <- read.table(opt$circos_f, header=TRUE, sep="\t")
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

cols_all <- c(viridis_pal(begin = 0, end = 1, option="D", alpha=1)(np_all), viridis_pal(begin = 0, end = 1, option="D", alpha=1)(nd_all))

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
filter_genes <- opt$filter_genes
if (filter_genes == 'true'){
  circos_data_genes <- read.table(opt$circos_g, header=TRUE, sep="\t")

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

  cols_genes <- c(viridis_pal(begin = 0, end = 1, option="D", alpha=1)(np_genes), viridis_pal(begin = 0, end = 1, option="D", alpha=1)(nd_genes))

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
