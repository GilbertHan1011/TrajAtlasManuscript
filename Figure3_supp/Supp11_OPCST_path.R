

lineage <- read.csv("../important_processed_data/11.2_dpt_lineage_infer.csv",row.names = 1)

lineage2 <- read.csv("../important_processed_data/11.",row.names = 1)
wt_integrate$C36_named <- Idents(wt_integrate)
wt_integrate$C36_named <- as.character(wt_integrate$C36_named)
wt_integrate$C36_named[wt_integrate$C36_named=="Ihh.Mature Chondro"]="Ihh+ HC"

wt_integrate$C36_named[wt_integrate$C36_named=="Lipc.HC"]="Osteocyte"

DimPlot(wt_integrate,cells.highlight = colnames(chondro),reduction = "X_draw_graph_fa")
ggsave("result/24.2.17_supp_chondro/chondro_lineage.pdf",width = 6,height = 4)
Idents(chondro) <- chondro$C36_named
chondro$C36_named[chondro$C36_named=="Lipc.HC"]="Osteocyte"
DimPlot(chondro,reduction = "X_draw_graph_fa",label = T)
ggsave("result/24.2.17_supp_chondro/chondro_label.pdf",width = 6,height = 4)



chondro <- wt_integrate[,rownames(lineage)[lineage$lineage_chondro=="True"]]

DimPlot(wt_integrate,reduction = "X_draw_graph_fa",label = T)
FeaturePlot(wt_integrate,c("Sost","Lipc",),reduction = "X_draw_graph_fa",label = T)
FeaturePlot(wt_integrate,c("Spp1"),reduction = "X_draw_graph_fa")
FeaturePlot(wt_integrate,c("Sost", "Ackr3", "Fbln7", "Dmp1", "Irx5"),reduction = "X_draw_graph_fa")

FeaturePlot(wt_integrate,c("Ddit4l","Loxl4","Fam20c","Cds1","F13a1","Spp1","Ibsp"),reduction = "X_draw_graph_fa")
FeaturePlot(wt_integrate,c("Ihh","Fxyd2","Il17b"),reduction = "X_draw_graph_fa")
FeaturePlot(wt_integrate,c("Bmp8a","Pla2g10"),reduction = "X_draw_graph_fa")
FeaturePlot(wt_integrate,c("Col10a1","Bmp8a","Ddit4l"),reduction = "X_draw_graph_fa",ncol = 3)
SCpubr::do_FeaturePlot(wt_integrate,features = c("Col10a1","Bmp8a","Ddit4l"),reduction = "X_draw_graph_fa",
                       ncol = 3,raster = T,raster.dpi = 600,pt.size = 2)
ggsave("result/24.2.17_supp_chondro/chondro_marker.pdf",width = 12,height = 5)
SCpubr::do_DimPlot(chondro,reduction = "X_draw_graph_fa",label = T)
SCpubr::do_DimPlot(wt_integrate,cells.highlight = colnames(chondro),reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/chondro_lineage.pdf",width = 4,height = 4)

SCpubr::do_DimPlot(wt_integrate,cells.highlight = colnames(wt_integrate)[wt_integrate$Project=="ChondroOsteo_Long"],reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/Col10a1_tracked.pdf",width = 4,height = 4)

SCpubr::do_FeaturePlot(wt_integrate,features = c("Col10a1"),reduction = "X_draw_graph_fa",raster = T,raster.dpi = 1000)
ggsave("result/24.2.17_supp_chondro/chondro_marker.pdf",width = 12,height = 4)
FeaturePlot(wt_integrate,c("Sost", "Ackr3", "Fbln7", "Dmp1", "Irx5"),reduction = "X_draw_graph_fa")



DimPlot(wt_integrate,cells.highlight = colnames(wt_integrate)[wt_integrate$Project=="ChondroOsteo_Long"],reduction = "X_draw_graph_fa")+ggtitle("Col10a1-tracked cells")
ggsave("result/24.2.17_supp_chondro/Col10a1_tracked.pdf",width = 6,height = 4)

wt_integrate$chondro_label <- wt_integrate$C36_named
indexes <- rownames(lineage)[lineage$lineage_chondro == "True"]
wt_integrate$chondro_label[!colnames(wt_integrate)%in%indexes]=NA


DimPlot(wt_integrate,reduction = "X_draw_graph_fa",group.by = "chondro_label")

SCpubr::do_DimPlot(wt_integrate,group.by = "chondro_label",reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/chondro_labeled.pdf",width = 6,height = 6)

#== Lepr BMSC

FeaturePlot(wt_integrate,c("Cxcl12","Angptl4","Ifitm3"),reduction = "X_draw_graph_fa")

SCpubr::do_FeaturePlot(wt_integrate,features = c("Cxcl12","Angptl4","Ifitm3"),reduction = "X_draw_graph_fa",
                       ncol = 3,raster = T,raster.dpi = 600,pt.size = 2)
ggsave("result/24.2.17_supp_chondro/lepr_marker.pdf",width = 12,height = 5)

SCpubr::do_DimPlot(wt_integrate,cells.highlight = rownames(lineage)[lineage$lineage_lepr=="True"],reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/lepr_lineage.pdf",width = 4,height = 4)


SCpubr::do_FeaturePlot(wt_integrate,features = c("Postn","Ctsk"),reduction = "X_draw_graph_fa",
                       ncol = 3,raster = T,raster.dpi = 600,pt.size = 2)
ggsave("result/24.2.17_supp_chondro/Cxcl12_tracked.pdf",width = 4,height = 4)

wt_integrate$lepr_label <- wt_integrate$C19_named
indexes <- rownames(lineage)[lineage$lineage_lepr == "True"]
wt_integrate$lepr_label[!colnames(wt_integrate)%in%indexes]=NA

SCpubr::do_DimPlot(wt_integrate,group.by = "lepr_label",reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/lepr_labeled.pdf",width = 6,height = 6)





#==Fibroblast---------------
DimPlot(wt_integrate,cells.highlight = rownames(lineage)[lineage$lineage_laFibro=="True"],reduction = "X_draw_graph_fa")
SCpubr::do_DimPlot(wt_integrate,cells.highlight =  rownames(lineage)[lineage$lineage_laFibro=="True"],reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/fibro_lineage.pdf",width = 4,height = 4)
SCpubr::do_FeaturePlot(wt_integrate,features = c("Postn","Ctsk"),reduction = "X_draw_graph_fa",
                        ncol = 2,raster = T,raster.dpi = 600,pt.size = 2,max.cutoff = c(4,4))
ggsave("result/24.2.17_supp_chondro/fibro_marker.pdf",width = 8,height = 5)


wt_integrate$fibro_label <- wt_integrate$C19_named
indexes <- rownames(lineage)[lineage$lineage_laFibro == "True"]
wt_integrate$fibro_label[!colnames(wt_integrate)%in%indexes]=NA

SCpubr::do_DimPlot(wt_integrate,group.by = "fibro_label",reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/fibro_labeled.pdf",width = 6,height = 6)




#== MSC-------------------
unique(wt_integrate$C19_named)

FeaturePlot(wt_integrate,c("Twist1","Prrx1"),reduction = "X_draw_graph_fa")

SCpubr::do_DimPlot(wt_integrate,cells.highlight =  rownames(lineage)[lineage$lineage_mesenchyme=="True"],reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/MSC_lineage.pdf",width = 4,height = 4)
SCpubr::do_FeaturePlot(wt_integrate,features = c("Prrx1","Twist1"),reduction = "X_draw_graph_fa",
                       ncol = 2,raster = T,raster.dpi = 600,pt.size = 2)
ggsave("result/24.2.17_supp_chondro/MSC_marker.pdf",width = 8,height = 5)

preob <- wt_integrate[,wt_integrate$C19_named=="Pre-ob"]
Idents(preob) <- preob$Organ

identity <- c("MSC OPCST","MSC OPCST","Lepr+ BMSC OPCST")
names(identity) <- levels(preob)
preob <- RenameIdents(preob, identity)
preobMarker <- FindMarkers(preob,ident.1  = "MSC OPCST")


wt_integrate$msc_label <- wt_integrate$C19_named
indexes <- rownames(lineage)[lineage$lineage_mesenchyme == "True"]
wt_integrate$msc_label[!colnames(wt_integrate)%in%indexes]=NA

SCpubr::do_DimPlot(wt_integrate,group.by = "msc_label",reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.17_supp_chondro/msc_labeled.pdf",width = 6,height = 6)


#== Osteocyte
FeaturePlot(wt_integrate,c("Sost", "Ackr3", "Irx5"),reduction = "X_draw_graph_fa")
SCpubr::do_FeaturePlot(wt_integrate,features = c("Irx5"),reduction = "X_draw_graph_fa",
                       raster = T,raster.dpi = 600,pt.size = 2)
ggsave("result/24.2.17_supp_chondro/Irx5.pdf",width = 6,height = 6)


