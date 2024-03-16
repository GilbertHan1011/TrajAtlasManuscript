#== function and package------------
library(tidyverse)
source("../function/anno_integrate.R")
source("../function//stratified_wilcoxon_functions.R")

#== read parameters-------------------------

parameter_list = jsonlite::read_json("script/annotation/config/parameters_wtintegrate.json")
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})


assay_markers = parameter_list$assay_markers
assay_slot = parameter_list$assay_slot
test.use = parameter_list$test.use
logfc.threshold = parameter_list$logfc.threshold
min.pct =parameter_list$min.pct
min.diff.pct = parameter_list$min.diff.pct
max.cells.per.ident = parameter_list$max.cells.per.ident
min.cells.feature = parameter_list$min.cells.feature
min.cells.group =  parameter_list$min.cells.group
base = parameter_list$base
only.pos = parameter_list$only.pos
use_stratified=parameter_list$use_stratified
batch_var = parameter_list$batch_var
test.use = parameter_list$test.use
message("test.use ",test.use," ... ",parameter_list$test.use)
if(is.null(base)){base = 2}


# load clusters
wtIntegrate_clusterint  = data.table::fread("processed_data/4.26_leiden_316.csv",data.table = F,)
#wtIntegrate_clusterint  = wtIntegrate_clusterint [,c(1,3:ncol(wtIntegrate_clusterint ))]


# load seurat
#curated_seurat_object = readRDS(paste0(parameter_list$harmonization_folder_path,parameter_list$new_name_suffix,"_curated",".rds"))

full_seurat = FindNeighbors(full_seurat,
                                      reduction="X_scANVI",
                                      dims = 1:ncol(full_seurat@reductions[["X_scANVI"]]@cell.embeddings),
                                      nn.method="annoy", graph.name = "SNN_scANVI",verbose=TRUE)

# get nn
neighbor_nn = as.Neighbor(full_seurat@graphs$SNN_scANVI)
nn_idx = neighbor_nn@nn.idx

#== load cluseter----------------------
full_seurat@meta.data$Cell_ID <- colnames(full_seurat)
wtIntegrate_clusterint$Cell_ID <-  colnames(full_seurat)
colnames(wtIntegrate_clusterint)[1] <- "cellid"
wtIntegrate_clusterint$Cell_ID <- colnames(full_seurat)
temp_meta = dplyr::left_join(full_seurat@meta.data,wtIntegrate_clusterint,by="Cell_ID")
rownames(temp_meta) = temp_meta$Cell_ID
full_seurat@meta.data = temp_meta
cluster_column = colnames(wtIntegrate_clusterint)[ncol(wtIntegrate_clusterint)]

wtIntegrate_clusterint_clear = apply(wtIntegrate_clusterint[,2:ncol(wtIntegrate_clusterint)],2,clear_clustering,min_cells = parameter_list$min_cells_valid,nn_idx = nn_idx) %>% as.data.frame()
wtIntegrate_clusterint_clear <- subset(wtIntegrate_clusterint_clear,select=-Cell_ID)
#length(table(wtIntegrate_clusterint_clear$leiden_clusters_0.1))
sort(apply(wtIntegrate_clusterint_clear,2,function(x){length(table(x))}))

# I drop th first level because it contains only 1 cluster
#wtIntegrate_clusterint_clear = wtIntegrate_clusterint_clear[,2:ncol(wtIntegrate_clusterint_clear)]

# add cell id
wtIntegrate_clusterint_clear = dplyr::bind_cols(Cell_ID = wtIntegrate_clusterint$Cell_ID,wtIntegrate_clusterint_clear)


# make summary of clusters per level
wtIntegrate_clusterint_long = wtIntegrate_clusterint_clear %>% dplyr::select(-Cell_ID) %>% tidyr::gather(key="cluster_level",value="cluster_id")
wtIntegrate_clusterint_long_summary = wtIntegrate_clusterint_long %>% dplyr::group_by(cluster_level) %>%
  dplyr::summarise(n=length(unique(cluster_id))) %>% dplyr::arrange(n)
wtIntegrate_clusterint_long_summary$cluster_level = factor(wtIntegrate_clusterint_long_summary$cluster_level,levels = wtIntegrate_clusterint_long_summary$cluster_level)

full_seurat$merge_id_level3

#== MODIFY cluster------------------
DimPlot(full_seurat,group.by = "leiden_clusters_0.1")
DimPlot(full_seurat,group.by = "leiden_clusters_0.25")
DimPlot(full_seurat,group.by = "leiden_clusters_0.75")
DimPlot(full_seurat,group.by = "leiden_clusters_2")
DimPlot(full_seurat,group.by = "leiden_clusters_5")
p1 <-  DimPlot(full_seurat,group.by = "mix_level_2")
p2 <- DimPlot(full_seurat,group.by = "leiden_clusters_0.05")
p1+p2

Idents(full_seurat) <- full_seurat$leiden_clusters_0.01
DimPlot(full_seurat,group.by = "leiden_clusters_0.05")
DefaultAssay(full_seurat) <- "RNA"
full_seurat@reductions$umap@assay.used <- "RNA"
full_seurat@reductions$draw_graph_fa@assay.used <- "RNA"
full_seurat <- subset(full_seurat,ident="0")
DimPlot(full_seurat,reduction = "umap",group.by = "leiden_clusters_0.05")
full_seurat$leiden_clusters_level_1 <- full_seurat$mix_level_1
full_seurat$leiden_clusters_level_2 <- full_seurat$mix_level_2
full_seurat$leiden_clusters_level_3 <- full_seurat$leiden_clusters_0.1
full_seurat$leiden_clusters_level_4 <- full_seurat$leiden_clusters_0.25
full_seurat$leiden_clusters_level_5 <- full_seurat$leiden_clusters_0.75
full_seurat$leiden_clusters_level_6 <- full_seurat$leiden_clusters_2
full_seurat$leiden_clusters_level_7 <- full_seurat$leiden_clusters_5



cluster_cols=c("leiden_clusters_level_1","leiden_clusters_level_2","leiden_clusters_level_3","leiden_clusters_level_4",
               "leiden_clusters_level_5","leiden_clusters_level_6","leiden_clusters_level_7")
mrtree_input_labels = full_seurat@meta.data[,cluster_cols]

apply(mrtree_input_labels,2,function(x){length(table(x))})
mrtree_input_labels$leiden_clusters_level_1
x_factor <- factor(mrtree_input_labels$leiden_clusters_level_1, levels = unique(mrtree_input_labels$leiden_clusters_level_1))
x_numeric <- as.numeric(x_factor)
mrtree_input_labels$leiden_clusters_level_1 <- x_numeric
x_factor <- factor(mrtree_input_labels$leiden_clusters_level_2, levels = unique(mrtree_input_labels$leiden_clusters_level_2))
x_numeric <- as.numeric(x_factor)
mrtree_input_labels$leiden_clusters_level_2 <- x_numeric

full_seurat@meta.data = full_seurat@meta.data[,!grepl("leiden_clusters_",colnames(full_seurat@meta.data))]
#add
full_seurat@meta.data = cbind(full_seurat@meta.data,mrtree_input_labels)

data.table::fwrite(mrtree_input_labels,"result/4.26_annotation/4.29_anno_mrtree_input_label.csv")


