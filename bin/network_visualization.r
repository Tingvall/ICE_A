#!/usr/bin/env Rscript

### PROCESS 13 MULTIPLE: CIRCOS PLOTS - ALL INTERACTIONS WITH FACTOR AND INTERACTIONS FILTERED FOR GENES ###

# Load libraries
require(RCy3)
require(circlize)
require(viridisLite)
require(optparse)


# Arguments
option_list <- list(make_option(c("--nodes"), type="character", default='NULL', help="Network nodes.", metavar="path"),
                    make_option(c("--edges"), type="character", default='NULL', help="Network edges.", metavar="path"),
                    make_option(c("--nodes_up"), type="character", default='NULL', help="Network nodes filtered for upregulated peaks.", metavar="path"),
                    make_option(c("--nodes_down"), type="character", default='NULL', help="Network nodes filtered for downregulated peaks.", metavar="path"),
                    make_option(c("--edges_up"), type="character", default='NULL', help="Network edges filtered for upregulated peaks.", metavar="path"),
                    make_option(c("--edges_down"), type="character", default='NULL', help="Network edges filtered for downregulated peaks.", metavar="path"),
                    make_option(c("--log2FC"), type="double", default='NULL', help="Log2FC threshold for differential peaks.", metavar="value"),
                    make_option(c("--padj"), type="double", default='NULL', help="Padj threshold for differential peaks.", metavar="value"),
                    make_option(c("--expression_log2FC"), type="double", default='NULL', help="Log2FC column for differential expression.", metavar="value"),
                    make_option(c("--expression_padj"), type="double", default='NULL', help="Padj column for differential expression.", metavar="character"),
                    make_option(c("--mode"), type="character", default='NULL', help="Define which mode to run the pipeline in. The options are basic (default), multiple or differential", metavar="character"),
                    make_option(c("--network_mode"), type="character", default='NULL', help="Defines mode network. Options are all (all interaction in the 2D-bed file), factor (all interaction with at least on peak overlap either anchor point) or genes (interactions associates with a gene list, provided by --genes).", metavar="value"),
                    make_option(c("--use_peakscore"), type="character", default='NULL', help="If set to true, peak scores will be used to set edge width in network visualization. Default: false.", metavar="character"),
                    make_option(c("--outdir"), type="character", default="./", help="Path for output files", metavar="path"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Loading and organizing data
mode=opt$mode
network_mode=opt$network_mode

if (mode=="differential"){
  nodes <- read.table(opt$nodes, header=TRUE, sep="\t", col.names=c("id", "type", "padj", "log2FC"))
  edges <- read.table(opt$edges, header=TRUE, sep="\t", col.names=c("source", "target", "score", "type"))
} else{
  nodes <- read.table(opt$nodes, header=TRUE, sep="\t", col.names=c("id", "type"))
  edges <- read.table(opt$edges, header=TRUE, sep="\t", col.names=c("source", "target","score", "type"))
}

edges[,"interaction"] <- "interacts"
edges[,"name"] <- paste(edges$source, "(interacts)", edges$target, sep=" ")
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

#Set edge width based on interaction score for interactions (and peak score for peaks if the argument use_peakscore is set to true)
use_peakscore = opt$use_peakscore
if (use_peakscore=="true"){
  edges_interaction <- edges[edges$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
  if (sum(! is.na(edges_interaction$score)) < 1){
    edges_interaction[is.na(edges_interaction$score), "score"] <- 1
  }
  else {
    edges_interaction[is.na(edges_interaction$score), "score"] <- min(edges_interaction[!is.na(edges_interaction$score), "score"])/2
    edges_interaction$score <- 9*((edges_interaction$score-min(edges_interaction$score))/(max(edges_interaction$score)-min(edges_interaction$score)))+1
  }
  edges_factor <- edges[edges$type %in% c("Factor-Distal", "Factor-Promoter"),]
  edges_factor$score <- 9*((edges_factor$score-min(edges_factor$score))/(max(edges_factor$score)-min(edges_factor$score)))+1
  edges_score <- rbind(edges_interaction, edges_factor)
} else{
  edges_score <- edges[edges$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
  if (sum(! is.na(edges_score$score)) < 1){
    edges_score[is.na(edges_score$score), "score"] <- 1
  }
  else {
    edges_score[is.na(edges_score$score), "score"] <- min(edges_score[!is.na(edges_score$score), "score"])/2
    edges_score$score <- 9*((edges_score$score-min(edges_score$score))/(max(edges_score$score)-min(edges_score$score)))+1
  }
}
setEdgePropertyBypass(edge.names=edges_score$name, new.values=edges_score$score, visual.property='EDGE_WIDTH',bypass = TRUE)


#For differntial mode:color by log2FC of peaks and GENES
if (mode=="differential"){
  map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  }
  pal_up <- viridisLite::rocket(100, begin = 0.5, end = 0.25, direction = 1)
  pal_down <- viridisLite::mako(100, begin = 0.6, end = 0.35, direction = -1)
  nodes_up_peak <- nodes[nodes$log2FC >= opt$log2FC & !is.na(nodes$log2FC) & nodes$padj <= opt$padj & !is.na(nodes$padj) & nodes$type %in% c("Distal", "Promoter"),]
  nodes_down_peak <- nodes[nodes$log2FC <= - opt$log2FC & !is.na(nodes$log2FC) & nodes$padj <= opt$padj & !is.na(nodes$padj) & nodes$type %in% c("Distal", "Promoter"),]
  nodes_up_gene <- nodes[nodes$log2FC >= opt$expression_log2FC & !is.na(nodes$log2FC) & nodes$padj <= opt$expression_padj & !is.na(nodes$padj) & nodes$type == "Gene",]
  nodes_down_gene <- nodes[nodes$log2FC <= - opt$expression_log2FC & !is.na(nodes$log2FC) & nodes$padj <= opt$expression_padj & !is.na(nodes$padj) & nodes$type == "Gene",]
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
    nodes_diff[,"color"] <- gsub("FF$", "", nodes_diff[,"color"])
    setNodePropertyBypass(nodes_diff$id,nodes_diff$color,'NODE_FILL_COLOR',bypass = TRUE)
  }
}
toggleGraphicsDetails()
exportImage("Network.pdf", 'PDF')
exportNetwork("Network.xgmml", type= 'xGMML')

#exportImage(paste(opt$outdir, "/Network.pdf", sep=""), 'PDF')
#exportNetwork(paste(opt$outdir, "/Network.xgmml", sep=""), type= 'xGMML')

if (mode=="differential" & network_mode=="differential"){
#Up
  up_nodes <- read.table(opt$nodes_up, header=TRUE, sep="\t", col.names=c("id", "type", "padj", "log2FC"))
  up_edges <- read.table(opt$edges_up, header=TRUE, sep="\t", col.names=c("source", "target", "score", "type"))
  up_edges[,"interaction"] <- "interacts"
  up_edges[,"name"] <- paste(up_edges$source, "(interacts)", up_edges$target, sep=" ")
  createNetworkFromDataFrames(up_nodes,up_edges, title="Network_up", collection="Networks" )
  setVisualStyle(style.name)
  lockNodeDimensions(FALSE, style.name)

  if (use_peakscore=="true"){
    up_edges_interaction <- up_edges[up_edges$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
    if (sum(! is.na(up_edges_interaction$score)) < 1){
      up_edges_interaction[is.na(up_edges_interaction$score, "score"), "score"] <- 1
    } else {
      up_edges_interaction[is.na(up_edges_interaction$score), "score"] <- min(up_edges_interaction[!is.na(up_edges_interaction$score)], "score")/2
      up_edges_interaction$score <- 9*((up_edges_interaction$score-min(up_edges_interaction$score))/(max(up_edges_interaction$score)-min(up_edges_interaction$score)))+1}
    up_edges_factor <- up_edges[up_edges$type %in% c("Factor-Distal", "Factor-Promoter"),]
    up_edges_factor$score <- 9*((up_edges_factor$score-min(up_edges_factor$score))/(max(up_edges_factor$score)-min(up_edges_factor$score)))+1
    up_edges_score <- rbind(up_edges_interaction, up_edges_factor)
  } else{
    up_edges_score <- up_edges[up_edges$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
    if (sum(! is.na(up_edges_score$score)) < 1){
      up_edges_score[is.na(up_edges_score$score), "score"] <- 1
    } else {
      up_edges_score[is.na(up_edges_score$score), "score"] <- min(up_edges_score[!is.na(up_edges_score$score)], "score")/2
      up_edges_score$score <- 9*((up_edges_score$score-min(up_edges_score$score))/(max(up_edges_score$score)-min(up_edges_score$score)))+1}
  }
  setEdgePropertyBypass(edge.names=up_edges_score$name, new.values=up_edges_score$score, visual.property='EDGE_WIDTH',bypass = TRUE)

  up_nodes_up_peak <- up_nodes[up_nodes$log2FC >= opt$log2FC & !is.na(up_nodes$log2FC) & up_nodes$padj <= opt$padj & !is.na(up_nodes$padj) & up_nodes$type %in% c("Distal", "Promoter"),]
  up_nodes_down_peak <- up_nodes[up_nodes$log2FC <= - opt$log2FC & !is.na(up_nodes$log2FC) & up_nodes$padj <= opt$padj & !is.na(up_nodes$padj) & up_nodes$type %in% c("Distal", "Promoter"),]
  up_nodes_up_gene <- up_nodes[up_nodes$log2FC >= opt$expression_log2FC & !is.na(up_nodes$log2FC) & up_nodes$padj <= opt$expression_padj & !is.na(up_nodes$padj) & up_nodes$type == "Gene",]
  up_nodes_down_gene <- up_nodes[up_nodes$log2FC <= - opt$expression_log2FC & !is.na(up_nodes$log2FC) & up_nodes$padj <= opt$expression_padj & !is.na(up_nodes$padj) & up_nodes$type == "Gene",]
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
    up_nodes_diff[,"color"] <- gsub("FF$", "", up_nodes_diff[,"color"])
    setNodePropertyBypass(up_nodes_diff$id,up_nodes_diff$color,'NODE_FILL_COLOR',bypass = TRUE)
  }
toggleGraphicsDetails()
exportImage("Network_up.pdf", 'PDF')
exportNetwork("Network_up.xgmml", type= 'xGMML')

#Down
  down_nodes <- read.table(opt$nodes_down, header=TRUE, sep="\t", col.names=c("id", "type", "padj", "log2FC"))
  down_edges <- read.table(opt$edges_down, header=TRUE, sep="\t", col.names=c("source", "target", "score", "type"))
  down_edges[,"interaction"] <- "interacts"
  down_edges[,"name"] <- paste(down_edges$source, "(interacts)", down_edges$target, sep=" ")
  createNetworkFromDataFrames(down_nodes,down_edges, title="Network_up", collection="Networks" )
  setVisualStyle(style.name)
  lockNodeDimensions(FALSE, style.name)

  if (use_peakscore=="true"){
    down_edges_interaction <- down_edges[down_edges$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
    if (sum(! is.na(down_edges_interaction$score)) < 1){
      down_edges_interaction[is.na(down_edges_interaction$score), "score"] <- 1
    } else {
      down_edges_interaction[is.na(down_edges_interaction$score), "score"] <- min(down_edges_interaction[!is.na(down_edges_interaction$score)], "score")/2
      down_edges_interaction$score <- 9*((down_edges_interaction$score-min(down_edges_interaction$score))/(max(down_edges_interaction$score)-min(down_edges_interaction$score)))+1}
    down_edges_factor <- down_edges[down_edges$type %in% c("Factor-Distal", "Factor-Promoter"),]
    down_edges_factor$score <- 9*((down_edges_factor$score-min(down_edges_factor$score))/(max(down_edges_factor$score)-min(down_edges_factor$score)))+1
    down_edges_score <- rbind(down_edges_interaction, down_edges_factor)
  } else{
    down_edges_score <- down_edges[down_edges$type %in% c("Distal-Promoter", "Promoter-Promoter"),]
    if (sum(! is.na(down_edges_score$score)) < 1){
      down_edges_score[is.na(down_edges_score$score), "score"] <- 1
    } else {
      down_edges_score[is.na(down_edges_score$score), "score"] <- min(down_edges_score[!is.na(down_edges_score$score)], "score")/2
      down_edges_score$score <- 9*((down_edges_score$score-min(down_edges_score$score))/(max(down_edges_score$score)-min(down_edges_score$score)))+1}
  }
  setEdgePropertyBypass(edge.names=down_edges_score$name, new.values=down_edges_score$score, visual.property='EDGE_WIDTH',bypass = TRUE)

  down_nodes_up_peak <- down_nodes[down_nodes$log2FC >= opt$log2FC & !is.na(down_nodes$log2FC) & down_nodes$padj <= opt$padj & !is.na(down_nodes$padj) & down_nodes$type %in% c("Distal", "Promoter"),]
  down_nodes_down_peak <- down_nodes[down_nodes$log2FC <= - opt$log2FC & !is.na(down_nodes$log2FC) & down_nodes$padj <= opt$padj & !is.na(down_nodes$padj) & down_nodes$type %in% c("Distal", "Promoter"),]
  down_nodes_up_gene <- down_nodes[down_nodes$log2FC >= opt$expression_log2FC & !is.na(down_nodes$log2FC) & down_nodes$padj <= opt$expression_padj & !is.na(down_nodes$padj) & down_nodes$type == "Gene",]
  down_nodes_down_gene <- down_nodes[down_nodes$log2FC <= - opt$expression_log2FC & !is.na(down_nodes$log2FC) & down_nodes$padj <= opt$expression_padj & !is.na(down_nodes$padj) & down_nodes$type == "Gene",]
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
    down_nodes_diff[,"color"] <- gsub("FF$", "", down_nodes_diff[,"color"])
    setNodePropertyBypass(down_nodes_diff$id,down_nodes_diff$color,'NODE_FILL_COLOR',bypass = TRUE)
  }
toggleGraphicsDetails()
exportImage("Network_down.pdf", 'PDF')
exportNetwork("Network_down.xgmml", type= 'xGMML')
}
