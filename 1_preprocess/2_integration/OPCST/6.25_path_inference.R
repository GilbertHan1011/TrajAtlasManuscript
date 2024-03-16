library(Seurat)

setwd("../../../3.9_wt_integrate/")
fullSeurat <- readRDS("../important_processed_data/5.4_wtintegrate_fullSeurat.Rds")
DimPlot(fullSeurat)

lineageMeta <- read.csv("../important_processed_data/6.10_lineageMeta.csv",row.names = 1)

lineageMeta_622 <- read.csv("../important_processed_data/6.22_lineageMeta.csv",row.names = 1)

laCell <- rownames(lineageMeta_622)[lineageMeta_622$lineage_laFibro]

absorbProb <- read.csv("../unimportant_processed_data/5.24_lineage.csv")
fullSeurat$absorb <- absorbProb$X1
FeaturePlot(fullSeurat,"absorb")
mesenchyme_cells <- colnames(fullSeurat)[fullSeurat$absorb>0.2& fullSeurat$Organ!="Limb_adult"]
DimPlot(fullSeurat,cells.highlight = mesenchyme_cells,reduction = "X_draw_graph_fa")


Idents(fullSeurat) <- fullSeurat$anno_level_3
leprCell <- colnames(fullSeurat)[fullSeurat$anno_level_3%in%c("LA_Lepr+ MSC","LA_Ppbp+ MSC")]
leprCell2 <- colnames(fullSeurat)[fullSeurat$anno_level_4%in%c("LA_Angpt4.Metaphysis MSC")]
DimPlot(fullSeurat,cells.highlight = colnames(fullSeurat)[fullSeurat$anno_level_3%in%c("LA_Ppbp+ MSC")],reduction = "X_draw_graph_fa")

DimPlot(fullSeurat,cells.highlight = colnames(fullSeurat)[fullSeurat$anno_level_4%in%c("LE_Satb2.Osteoblast")],reduction = "X_draw_graph_fa")
DimPlot(fullSeurat,cells.highlight = colnames(fullSeurat)[fullSeurat$anno_level_4%in%c("LA_Cytl1.Metaphysis MSC")],reduction = "X_draw_graph_fa")
osteoCell <- colnames(fullSeurat)[fullSeurat$anno_level_3%in%c("LA_Osteoblast")]
osteoCell1 <- sample(osteoCell,size = 11000)

osteoCell2 <- setdiff(osteoCell,osteoCell1)
leprCell <- c(leprCell,leprCell2,osteoCell)


laFibro_cells <- c(laCell,osteoCell)


DimPlot(fullSeurat,cells.highlight = laFibro_cells,reduction = "X_draw_graph_fa")
DimPlot(fullSeurat,cells.highlight = laFibro_cells)


absorbProb2 <- read.csv("../unimportant_processed_data/5.24_lineage_la.csv")
fullSeurat$absorb_la <- absorbProb2$X1
FeaturePlot(fullSeurat,"absorb_la")
mesenchyme_cells <- colnames(fullSeurat)[fullSeurat$absorb>0.2& fullSeurat$Organ!="Limb_adult"]
DimPlot(fullSeurat,cells.highlight = mesenchyme_cells,reduction = "X_draw_graph_fa")


DimPlot(fullSeurat,cells.highlight = leprCell,reduction = "X_draw_graph_fa")
laFibro_cells <- colnames(fullSeurat)[fullSeurat$absorb_la>0.2& fullSeurat$Organ=="Limb_adult"&fullSeurat$anno_level_2!="LA_Diaphyseal MSC"]
laFibro_cells
DimPlot(fullSeurat,cells.highlight = laFibro_cells,reduction = "X_draw_graph_fa")
DimPlot(fullSeurat,cells.highlight = laFibro_cells)



Idents(fullSeurat) <- fullSeurat$C19_named
DimPlot(fullSeurat,label = T)


chondroCell <- colnames(fullSeurat)[fullSeurat$C19_named%in%c("Mature Chondro","HC","Tex14+ HC")&fullSeurat$Organ!="Head"]
chondroCell <- c(chondroCell,colnames(fullSeurat)[fullSeurat$C36_named=="Bglap3.Ob"&fullSeurat$Organ!="Head"])


mesenchyme_cells <- lineageMeta$lineage_mesenchyme
DimPlot()
DimPlot(fullSeurat,reduction = "X_draw_graph_fa",cells.highlight = chondroCell)+ggtitle("Chondro_lineage")
fullSeurat$lineage_chondro <- FALSE
fullSeurat$lineage_chondro[chondroCell] <- TRUE
fullSeurat$lineage_laFibro <- FALSE
fullSeurat$lineage_laFibro[laFibro_cells] <- TRUE
fullSeurat$lineage_lepr <- FALSE
fullSeurat$lineage_lepr[leprCell] <- TRUE
fullSeurat$lineage_mesenchyme <- FALSE
fullSeurat$lineage_mesenchyme[mesenchyme_cells] <- TRUE
p1 <- DimPlot(fullSeurat,reduction = "X_draw_graph_fa",cells.highlight = chondroCell)+ggtitle("Chondro_lineage")
p2 <- DimPlot(fullSeurat,reduction = "X_draw_graph_fa",cells.highlight = laFibro_cells)+ggtitle("Fibro_lineage")
p3 <- DimPlot(fullSeurat,reduction = "X_draw_graph_fa",cells.highlight = leprCell)+ggtitle("BMSC_lineage")
p4 <- DimPlot(fullSeurat,reduction = "X_draw_graph_fa",cells.highlight = mesenchyme_cells)+ggtitle("Mesenchyme_lineage")
p1+p2+p3+p4
ggsave("result/5.23_pesduotime/6.22_lineage_plot.pdf",width = 12,height = 8)
lineage_meta <-fullSeurat@meta.data[grep("^lineage", colnames(fullSeurat@meta.data))]
write.csv(lineage_meta,"../important_processed_data/6.25_lineageMeta.csv")

UpSetR::upset(fromList(list(mes = mesenchyme_cells,laFibro = laFibro_cells,  chondroCell= chondroCell,lepr=leprCell)))


FeaturePlot(fullSeurat,c("Alpl","Cpe","Apoe","Meg3"),ncol = 2,cols = c("grey","red"),reduction = "X_draw_graph_fa")
ggsave("result/6.20_traj_program/gene_program_expression.pdf",width = 10,height = 8)
