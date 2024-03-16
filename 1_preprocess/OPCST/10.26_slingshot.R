slingshot_lineage
library(slingshot)
library(tidyverse)
slingshot_lineage_la <- slingshot(full_seurat@reductions$X_draw_graph_fa@cell.embeddings,
                                clusterLabels = full_seurat$start,
                                start.clus = 'LimbAdult_Fibro',end.clus=c("Osteoblast"))
lineage_la <- slingshot_lineage_la@assays@data$weights[,2]
full_seurat$lineage_la <- lineage_la
FeaturePlot(full_seurat,"lineage_la",reduction = "X_draw_graph_fa",raster = T)

slingshot_lineage_head <- slingshot(full_seurat@reductions$X_draw_graph_fa@cell.embeddings,
                                  clusterLabels = full_seurat$start,
                                  start.clus = 'Head_Mesenchyme',end.clus=c("Osteoblast"))
full_seurat$lineage_head <- slingshot_lineage_head@assays@data$weights[,1]

FeaturePlot(full_seurat,"ob",reduction = "X_draw_graph_fa",raster = T)
absorbProb2 <- read.csv("../temp_data/5.24_lineage_la.csv")

full_seurat$absorb_la <- absorbProb2$X1
FeaturePlot(full_seurat,"absorb_la",raster = T,,reduction = "X_draw_graph_fa")

laFibro_cells <- colnames(full_seurat)[full_seurat$lineage_la>0.1& full_seurat$Organ=="Limb_adult"]


Idents(full_seurat) <- full_seurat$C7_named
DimPlot(full_seurat,cells = laFibro_cells,raster = T,reduction = "X_draw_graph_fa")

slingshot_lineage_bmsc <- slingshot(full_seurat@reductions$X_draw_graph_fa@cell.embeddings,
                                    clusterLabels = full_seurat$start,
                                    start.clus = 'LimbAdult_BMSC',end.clus=c("Osteoblast"))

slingshot_lineage_bmsc@metadata$lineages
full_seurat$lineage_bmsc <- slingshot_lineage_bmsc@assays@data$weights[,2]
FeaturePlot(full_seurat,"lineage_bmsc",raster = T,reduction = "X_draw_graph_fa")



slingshot_lineage_hc <- slingshot(full_seurat@reductions$X_draw_graph_fa@cell.embeddings,
                                    clusterLabels = full_seurat$start,
                                    start.clus = 'Limb_Hyperchondrocyte',end.clus=c("Osteoblast"))

slingshot_lineage_hc@metadata$lineages
full_seurat$lineage_hc <- slingshot_lineage_hc@assays@data$weights[,2]
FeaturePlot(full_seurat,"lineage_hc",raster = T,reduction = "X_draw_graph_fa")
lineage_infer_slingshot <- full_seurat@meta.data%>%select(starts_with("lineage"))
write.csv(lineage_infer_slingshot,"../important_processed_data/10.27_slingshot_lineage.csv")
