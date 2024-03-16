setwd("../3.9_wt_integrate/")

Ob <- full_seurat[,full_seurat$C7_named=="Ob"]
FeaturePlot(Ob,"Ccnj",order = T,split.by = "Stage")
FeaturePlot(Ob,"Amt",order = T,split.by = "Stage")
VlnPlot(Ob,group.by = "Stage",features = "Tra2b")
FeaturePlot(Ob,"Tra2b",order = T,split.by = "Stage")
FeaturePlot(Ob,"Morc4",order = T,split.by = "Stage")
VlnPlot(Ob,group.by = "Stage",features = "Morc4")
FeaturePlot(Ob,"Lrrc17",order = T,split.by = "Stage")
FeaturePlot(Ob,"Lrrc17",order = T,split.by = "Stage")

full_seurat$celltype_Project <- paste0(full_seurat$paper_label,"-",full_seurat$Project)
dput(unique(full_seurat$celltype_Project))
chondroOsteoLong <- readRDS("../data/cleandata/")
chondroOsteoLongMeta <- readRDS("../1.10_coreQC/result/ChondroOsteo_Long_meta_bk.Rds")
head(colnames(chondroSub))

unique(chondroOsteoLongMeta$paper_label)
head(colnames(chondroOsteoLong))
chondroCellName <- colnames(chondroSub)
chondroCellName <- gsub("ChondroOsteo_Long_","",chondroCellName)
chondroCellName <- gsub("-[1-9]$","",chondroCellName)
chondroPaper <- chondroOsteoLongMeta[chondroCellName]$paper_label
metaName <- gsub("SeuratProject_","",rownames(chondroOsteoLongMeta))
rownames(chondroOsteoLongMeta) <- metaName
chondroPaper <- chondroOsteoLongMeta[chondroCellName,]$paper_label

chondroSub <- full_seurat[,full_seurat$Project=="ChondroOsteo_Long"]
full_seurat$paper_label <- factor(full_seurat$paper_label,levels = c(levels(full_seurat$paper_label),as.character(unique(chondroPaper))))


full_seurat$paper_label[colnames(chondroSub)] <- chondroPaper

Idents(full_seurat) <- full_seurat$celltype_Project
chondro1 <- WhichCells(full_seurat, idents = c("Chondrocytes-progenitors-Bmsc2019_Regev"))
chondro2 <- WhichCells(full_seurat, idents = c("Chondrocytes-Bmsc2019_Regev"))
chondro3 <- WhichCells(full_seurat, idents = c("Chondrocytes-proliferating-resting-Bmsc2019_Regev"))
chondro4 <- WhichCells(full_seurat, idents = c("Chondrocytes-hypertrophic-Bmsc2019_Regev"))
chondro5 <- WhichCells(full_seurat, idents = c("Chondrocytes-prehypertrophic-Bmsc2019_Regev"))
b_chondro1 <- WhichCells(full_seurat, idents = c("Hypertrophic Chondrocytes 1-ChondroOsteo_Long"))
b_chondro2 <- WhichCells(full_seurat, idents = c("Hypertrophic Chondrocytes 2-ChondroOsteo_Long"))
b_chondro3 <- WhichCells(full_seurat, idents = c("Hypertrophic Chondrocytes 3-ChondroOsteo_Long"))
b_chondro4 <- WhichCells(full_seurat, idents = c("Terminal-hypertrophic Chondrocytes 1-ChondroOsteo_Long"))
b_chondro5 <- WhichCells(full_seurat, idents = c("Terminal-hypertrophic Chondrocytes 2-ChondroOsteo_Long"))
c_chondro <-  WhichCells(full_seurat, idents = c("Growth plate-LimbMouse2019_Kelly"))
chondrolist=list(chondro1, chondro2,chondro3,chondro4,chondro5,b_chondro1,b_chondro2,b_chondro3,b_chondro4,b_chondro5)
names(chondrolist) <- c("Chondrocytes-progenitors-Bmsc2019_Regev","Chondrocytes-Bmsc2019_Regev","Chondrocytes-proliferating-resting-Bmsc2019_Regev",
                        "Chondrocytes-hypertrophic-Bmsc2019_Regev","Chondrocytes-prehypertrophic-Bmsc2019_Regev",
                        "Hypertrophic Chondrocytes 1-ChondroOsteo_Long","Hypertrophic Chondrocytes 2-ChondroOsteo_Long","Hypertrophic Chondrocytes 3-ChondroOsteo_Long",
                        "Terminal-hypertrophic Chondrocytes 1-ChondroOsteo_Long","Terminal-hypertrophic Chondrocytes 2-ChondroOsteo_Long")
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(11)
DimPlot(full_seurat, cells.highlight= chondrolist, cols.highlight = mycolor, cols= "grey",raster = F)
a_chondro1 <- WhichCells(full_seurat, idents = c("Superficial ACC-GrowthplateSox9_Abdul"))
a_chondro2 <- WhichCells(full_seurat, idents = c("Deep ACC-GrowthplateSox9_Abdul"))
a_chondro3 <- WhichCells(full_seurat, idents = c("Proliferating ACC-GrowthplateSox9_Abdul"))
a_chondro4 <- WhichCells(full_seurat, idents = c("Pre/Hypertrophic GPC-GrowthplateSox9_Abdul"))
a_chondro6 <- WhichCells(full_seurat, idents = c("Hypertrophic/terminal GPC-GrowthplateSox9_Abdul"))
a_chondro5 <- WhichCells(full_seurat, idents = c("Chondrocytes-GrowthplateSox9_Abdul"))
a_chondrolist=list(a_chondro1, a_chondro2,a_chondro3,a_chondro4,a_chondro5,a_chondro6)
names(a_chondrolist) <- c("Superficial ACC-GrowthplateSox9_Abdul","Deep ACC-GrowthplateSox9_Abdul","Proliferating ACC-GrowthplateSox9_Abdul",
                          "Pre/Hypertrophic GPC-GrowthplateSox9_Abdul","Chondrocytes-GrowthplateSox9_Abdul","Hypertrophic/terminal GPC-GrowthplateSox9_Abdul")
mycolor2<-colorRampPalette(brewer.pal(8,'Spectral'))(6)
DimPlot(full_seurat, cells.highlight= a_chondrolist, cols.highlight = mycolor2, cols= "grey",raster = F)
chondrolist=list(chondro1, chondro2,chondro3,chondro4,chondro5,b_chondro1,b_chondro2,b_chondro3,b_chondro4,b_chondro5,a_chondro2,a_chondro5,a_chondro6)
names(chondrolist) <- c("Chondrocytes-progenitors-Bmsc2019_Regev","Chondrocytes-Bmsc2019_Regev","Chondrocytes-proliferating-resting-Bmsc2019_Regev",
                        "Chondrocytes-hypertrophic-Bmsc2019_Regev","Chondrocytes-prehypertrophic-Bmsc2019_Regev",
                        "Hypertrophic Chondrocytes 1-ChondroOsteo_Long","Hypertrophic Chondrocytes 2-ChondroOsteo_Long","Hypertrophic Chondrocytes 3-ChondroOsteo_Long",
                        "Terminal-hypertrophic Chondrocytes 1-ChondroOsteo_Long",
                        "Terminal-hypertrophic Chondrocytes 2-ChondroOsteo_Long","Deep ACC-GrowthplateSox9_Abdul","Chondrocytes-GrowthplateSox9_Abdul",
                        "Hypertrophic/terminal GPC-GrowthplateSox9_Abdul")
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(length(chondrolist))
#mycolor<-colorRampPalette( pal_npg("nrc")(10))(length(chondrolist))
DimPlot(full_seurat, cells.highlight= chondrolist, cols.highlight = mycolor, cols= "grey",raster = F)
ggsave("result/10.3_umap_vis/chondro_highRes.pdf",width = 10,height = 8)
DimPlot(full_seurat, cells.highlight= chondrolist, cols.highlight = mycolor, cols= "grey",raster = T)
ggsave("result/10.3_umap_vis/chondro_highRes_raster.pdf",width = 12,height = 8)

lepr1 <-  WhichCells(full_seurat, idents = c("Diaphyseal MSCs 2-Septoclasts_Kishor"))
lepr2 <-  WhichCells(full_seurat, idents = c("Diaphyseal MSCs 1-Septoclasts_Kishor"))
lepr3 <-  WhichCells(full_seurat, idents = c("Metaphyseal MSCs-Septoclasts_Kishor"))
lepr4 <-  WhichCells(full_seurat, idents = c("Proliferating BMSCs-BMSC-Specification_Kishor"))
leprList <-list(lepr1, lepr2,lepr3,lepr4)
names(leprList) <- c("Diaphyseal MSCs 1-Septoclasts_Kishor","Diaphyseal MSCs 2-Septoclasts_Kishor","Metaphyseal MSCs-Septoclasts_Kishor","Proliferating BMSCs-Septoclasts_Kishor")
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(length(leprList))
DimPlot(full_seurat, cells.highlight= leprList, cols.highlight = mycolor, cols= "grey",raster = T)
ggsave("result/10.3_umap_vis/lepr_highRes_raster.pdf",width = 12,height = 8)
