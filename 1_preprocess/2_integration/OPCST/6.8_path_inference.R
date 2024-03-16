library(Seurat)

absorbProb <- read.csv("../unimportant_processed_data/5.24_lineage.csv")
full_seurat$absorb <- absorbProb$X1
FeaturePlot(full_seurat,"absorb")
mesenchyme_cells <- colnames(full_seurat)[full_seurat$absorb>0.2& full_seurat$Organ!="Limb_adult"]
DimPlot(full_seurat,cells.highlight = mesenchyme_cells,reduction = "X_draw_graph_fa")


Idents(full_seurat) <- full_seurat$anno_level_3
leprCell <- colnames(full_seurat)[full_seurat$anno_level_3%in%c("LA_Lepr+ MSC","LA_Ppbp+ MSC")]
leprCell2 <- colnames(full_seurat)[full_seurat$anno_level_4%in%c("LA_Angpt4.Metaphysis MSC")]
DimPlot(full_seurat,cells.highlight = colnames(full_seurat)[full_seurat$anno_level_3%in%c("LA_Ppbp+ MSC")],reduction = "X_draw_graph_fa")

DimPlot(full_seurat,cells.highlight = colnames(full_seurat)[full_seurat$anno_level_4%in%c("LE_Satb2.Osteoblast")],reduction = "X_draw_graph_fa")
DimPlot(full_seurat,cells.highlight = colnames(full_seurat)[full_seurat$anno_level_4%in%c("LA_Cytl1.Metaphysis MSC")],reduction = "X_draw_graph_fa")
osteoCell <- colnames(full_seurat)[full_seurat$anno_level_3%in%c("LA_Osteoblast")]
osteoCell1 <- sample(osteoCell,size = 11000)
osteoCell2 <- setdiff(osteoCell,osteoCell1)
leprCell <- c(leprCell,leprCell2,osteoCell1)

DimPlot(full_seurat,cells.highlight = leprCell,reduction = "X_draw_graph_fa")
laFibro_cells <- colnames(full_seurat)[full_seurat$absorb>0.2& full_seurat$Organ=="Limb_adult"&full_seurat$anno_level_2!="LA_Chondrocyte"&
                                         !full_seurat$anno_level_3%in%c("LA_Lepr+ MSC","LA_Ppbp+ MSC","LA_Osteoblast","LA_Hypertrophic.Chondrocyte")&
                                         !full_seurat$anno_level_4%in%c("LA_Vcan.Metaphysis MSC","LA_Angpt4.Metaphysis MSC")]
laFibro_cells <- c(laFibro_cells,osteoCell2)

DimPlot(full_seurat,cells.highlight = laFibro_cells,reduction = "X_draw_graph_fa")
DimPlot(full_seurat,cells.highlight = laFibro_cells)


absorbProb2 <- read.csv("../unimportant_processed_data/5.24_lineage_la.csv")
full_seurat$absorb_la <- absorbProb2$X1
FeaturePlot(full_seurat,"absorb_la")
mesenchyme_cells <- colnames(full_seurat)[full_seurat$absorb>0.2& full_seurat$Organ!="Limb_adult"]
DimPlot(full_seurat,cells.highlight = mesenchyme_cells,reduction = "X_draw_graph_fa")


DimPlot(full_seurat,cells.highlight = leprCell,reduction = "X_draw_graph_fa")
laFibro_cells <- colnames(full_seurat)[full_seurat$absorb_la>0.2& full_seurat$Organ=="Limb_adult"&full_seurat$anno_level_2!="LA_Diaphyseal MSC"]
laFibro_cells
DimPlot(full_seurat,cells.highlight = laFibro_cells,reduction = "X_draw_graph_fa")
DimPlot(full_seurat,cells.highlight = laFibro_cells)



Idents(full_seurat) <- full_seurat$C19_named
DimPlot(full_seurat,label = T)


chondroCell <- colnames(full_seurat)[full_seurat$C19_named%in%c("Mature Chondro","HC","Tex14+ HC")&full_seurat$Organ!="Head"]
chondroCell <- c(chondroCell,colnames(full_seurat)[full_seurat$C36_named%in%c("Bglap3.Ob")]&full_seurat$Organ!="Head"&full_seurat$Organ!="Limb_adult")
DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = chondroCell)+ggtitle("Chondro_lineage")
full_seurat$lineage_chondro <- FALSE
full_seurat$lineage_chondro[chondroCell] <- TRUE
full_seurat$lineage_laFibro <- FALSE
full_seurat$lineage_laFibro[laFibro_cells] <- TRUE
full_seurat$lineage_lepr <- FALSE
full_seurat$lineage_lepr[leprCell] <- TRUE
full_seurat$lineage_mesenchyme <- FALSE
full_seurat$lineage_mesenchyme[mesenchyme_cells] <- TRUE
p1 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = chondroCell)+ggtitle("Chondro_lineage")
p2 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = laFibro_cells)+ggtitle("Fibro_lineage")
p3 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = leprCell)+ggtitle("BMSC_lineage")
p4 <- DimPlot(full_seurat,reduction = "X_draw_graph_fa",cells.highlight = mesenchyme_cells)+ggtitle("Mesenchyme_lineage")
p1+p2+p3+p4
ggsave("result/5.23_pesduotime/6.22_lineage_plot.pdf",width = 12,height = 8)
lineage_meta <-full_seurat@meta.data[grep("^lineage", colnames(full_seurat@meta.data))]
write.csv(lineage_meta,"../important_processed_data/6.22_lineageMeta.csv")

UpSetR::upset(fromList(list(mes = mesenchyme_cells,laFibro = laFibro_cells,  chondroCell= chondroCell,lepr=leprCell)))


FeaturePlot(full_seurat,c("Alpl","Cpe","Apoe","Meg3"),ncol = 2,cols = c("grey","red"),reduction = "X_draw_graph_fa")
ggsave("result/6.20_traj_program/gene_program_expression.pdf",width = 10,height = 8)
            