slingshotLienage <- read.csv("../important_processed_data/10.27_slingshot_lineage.csv")
library(tidyverse)
full_seurat$lineage_la <- slingshotLienage$lineage_la
laFibro_cells <- colnames(full_seurat)[full_seurat$lineage_la>0& full_seurat$Organ=="Limb_adult"&full_seurat$anno_level_2!="LA_Diaphyseal MSC"]
DimPlot(full_seurat,cells.highlight = laFibro_cells,reduction ="X_draw_graph_fa" )
lineageMeta <- read.csv("../important_processed_data/6.22_lineageMeta.csv")
lineageMeta$lineage_laFibro <- FALSE
lineageMeta <- lineageMeta%>%column_to_rownames("X")
lineageMeta[laFibro_cells,]$lineage_laFibro <- TRUE

lineageMeta_0 <- read.csv("../important_processed_data/6.10_lineageMeta.csv")
lineageMeta$lineage_lepr

Idents(full_seurat) <- full_seurat$anno_level_3
leprCell <- colnames(full_seurat)[full_seurat$anno_level_3%in%c("LA_Lepr+ MSC","LA_Ppbp+ MSC")]
leprCell2 <- colnames(full_seurat)[full_seurat$anno_level_4%in%c("LA_Angpt4.Metaphysis MSC")]
DimPlot(full_seurat,cells.highlight = colnames(full_seurat)[full_seurat$anno_level_3%in%c("LA_Ppbp+ MSC")],reduction = "X_draw_graph_fa")

DimPlot(full_seurat,cells.highlight = colnames(full_seurat)[full_seurat$anno_level_4%in%c("LE_Satb2.Osteoblast")],reduction = "X_draw_graph_fa")
DimPlot(full_seurat,cells.highlight = colnames(full_seurat)[full_seurat$anno_level_4%in%c("LA_Cytl1.Metaphysis MSC")],reduction = "X_draw_graph_fa")
osteoCell <- colnames(full_seurat)[full_seurat$anno_level_3%in%c("LA_Osteoblast")]
osteoCell1 <- sample(osteoCell,size = 11000)

osteoCell2 <- setdiff(osteoCell,osteoCell1)
leprCell <- c(leprCell,leprCell2,osteoCell)
DimPlot(full_seurat,cells.highlight = leprCell,reduction = "X_draw_graph_fa")
lineageMeta$lineage_lepr <- FALSE
lineageMeta[leprCell,]$lineage_lepr <- TRUE
chondroCell <- rownames(lineageMeta)[lineageMeta$lineage_chondro==TRUE]
mesenchymeCell <- rownames(lineageMeta)[lineageMeta$lineage_mesenchyme==TRUE]
p1 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = chondroCell)+ggtitle("Chondro_lineage")
p2 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = laFibro_cells)+ggtitle("Fibro_lineage")
p3 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = leprCell)+ggtitle("BMSC_lineage")
p4 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = mesenchymeCell)+ggtitle("Mesenchyme_lineage")
p1+p2+p3+p4
