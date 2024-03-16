library(MuDataSeurat)
trajMap <- ReadH5MU("../important_processed_data/12.26_trajMap_v2.h5mu")
Idents(trajMap) <- trajMap$Lineage
DimPlot(trajMap)
mypal <- pal_npg("nrc", alpha = 0.7)(4)
names(mypal) <- unique(Idents(trajMap))
SCpubr::do_DimPlot(sample = trajMap,pt.size = 5,reduction = "UMAP",colors.use = mypal)
ggsave('result/12.26_fig4/umap_lineage.pdf',width = 5,height = 5)

Idents(trajMap) <- trajMap$Organ
mypal <- pal_npg("nrc", alpha = 0.7)(3)
names(mypal) <- unique(Idents(trajMap))
SCpubr::do_DimPlot(sample = trajMap,pt.size = 8,reduction = "UMAP",colors.use = mypal,border.size=1.3)
ggsave('result/12.26_fig4/umap_organ.pdf',width = 5,height = 5)


#== plot TRAV heatmap-------------------
DefaultAssay(trajMap) <- "TRAV"
Idents(trajMap) <- trajMap$Lineage
travMarker <- FindAllMarkers(trajMap,group.by = "Lineage",)
coor <- read.csv("processed_data/12.9_trajMap/12.9_coorDf.csv",row.names = 1)
peak <- read.csv("processed_data/12.9_trajMap/12.9_maxDf.csv",row.names = 1)
coorMat <- as.matrix(coor)
colnames(valAll) <- shortName2

namesDesend <- names(sort(valAll["RAV138",],decreasing = T))
geneDesendDf <- RAVindex[,"RAV138"]
geneDesendDf <- geneDesendDf[geneDesendDf>0]
geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)]

Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F)
Heatmap(peak[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F)

IdentidyBar <- trajMap$Lineage[namesDesend]
my_color_lineage<- pal_npg("nrc", alpha = 0.7)(4)
names(my_color_lineage) <- unique(IdentidyBar)
haAll = HeatmapAnnotation(
  Lineage = IdentidyBar,
  col = list(
    Lineage=my_color_lineage
  )
)

geneVal <- RAVindex[,"RAV138"]
geneVal <- geneVal[geneDesend]

geneVal <- scale(geneVal)
colorpalette <- brewer.pal(3, "Blues")
row_ha = rowAnnotation(NMF =geneVal,col=list(NMF=colorRamp2(c(-2, 0,4), c("white","lightblue","blue"))))

labelsT <- c("Upk1b")
position <- c()
for (i in labelsT){
  position <- c(position, which(is.element(rownames(coor[geneDesend,namesDesend]), i) == TRUE))
}
har <- rowAnnotation(link = anno_mark(at = position, labels = labelsT, 
                                      labels_gp = gpar(fontsize = 9, fontface = "bold"), link_width = unit(1.5, "cm")))

hm <- Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F,top_annotation = haAll,
        left_annotation = row_ha,right_annotation = har,use_raster = T,col=colorRamp2(c(-1, 0, 1), c("DeepSkyBlue3", "white", "red")))
pdf("result/12.26_fig4/12.26_heatmap_Lepr_OCP.pdf",width = 8,height = 8)
hm <- draw(hm)
dev.off()
RAVindex <- RAVmodel@assays@data$RAVindex
colnames(RAVindex) <- rownames(valAll)


filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Bmsc2019_Regev_samp4.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))
pdf("result/12.26_fig4/12.26_heatmap_NMF_Bmsc2019_Regev.pdf")
Heatmap(filterMat,cluster_columns = F,cluster_rows = F,show_row_names = F,
        show_column_names = F,use_raster = T,left_annotation = row_ha,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
dev.off()
        

GOres <- readRDS("processed_data/12.12_TRAV/geneset/gsea_var_nmf/gsea_138.rds")
clusterProfiler::dotplot(GOres)


filterMat2 <- read.csv("processed_data/12.12_TRAV/filter/mesenchyme//LimbG610C_Gorrell_femurWT2.csv",row.names = 1)
filterMat2 <- filterMat2[geneDesend,]
filterMat2[is.na(filterMat2)] <- 0
filterMat2 <- t(scale(t(filterMat2)))
#pdf("result/12.26_fig4/12.26_heatmap_NMF_Bmsc2019_Regev.pdf")
Heatmap(filterMat2,cluster_columns = F,cluster_rows = F,show_row_names = F,
        show_column_names = F,use_raster = T,left_annotation = row_ha,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
dev.off()



clusterProfiler::dotplot(GOres)




i="RAV143"
namesDesend <- names(sort(valAll[i,],decreasing = T))
geneDesendDf <- RAVindex[,i]
geneDesendDf <- geneDesendDf[geneDesendDf>0]
geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)]

# Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F)
# Heatmap(peak[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F)

filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Septoclasts_Kishor_Pdgfra.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))

IdentidyBar <- trajMap$Lineage[namesDesend]
my_color_lineage<- pal_npg("nrc", alpha = 0.7)(4)
names(my_color_lineage) <- unique(IdentidyBar)
haAll = HeatmapAnnotation(
  Lineage = IdentidyBar,
  col = list(
    Lineage=my_color_lineage
  )
)

geneVal <- RAVindex[,i]
geneVal <- geneVal[geneDesend]

geneVal <- scale(geneVal)
colorpalette <- brewer.pal(3, "Blues")
row_ha = rowAnnotation(NMF =geneVal,col=list(NMF=colorRamp2(c(-2, 0,4), c("white","lightblue","blue"))))

labelsT <- c("F3","Lrp4")
position <- c()
for (i in labelsT){
  position <- c(position, which(is.element(rownames(coor[geneDesend,namesDesend]), i) == TRUE))
}
har <- rowAnnotation(link = anno_mark(at = position, labels = labelsT, 
                                      labels_gp = gpar(fontsize = 9, fontface = "bold"), link_width = unit(1.5, "cm")))

hm <- Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F,top_annotation = haAll,
              left_annotation = row_ha,use_raster = T,col=colorRamp2(c(-1, 0, 1),hcl_palette="PiYG",reverse = TRUE),raster_quality = 10)
pdf("result/12.26_fig4/12.26_heatmap_Lepr_OCP_RAV143.pdf",width = 8,height = 8)
hm <- draw(hm)
dev.off()
RAVindex <- RAVmodel@assays@data$RAVindex
colnames(RAVindex) <- rownames(valAll)





haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(filterMat)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)


pdf("result/12.26_fig4/12.26_heatmap_RAV143_Septoclasts_Kishor_Pdgfra.pdf")
hm2 <- Heatmap(filterMat,cluster_columns = F,cluster_rows = F,show_row_names = F,
        show_column_names = F,use_raster = T,right_annotation = har,top_annotation = haPseudutime,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
hm2
dev.off()

hmRes <- hm+hm2
pdf("result/12.26_fig4/OPC_specific_hm.pdf",width = 8,height = 6)
draw(hmRes)
dev.off()

write.csv(geneDesend,"../temp_data/12.28_geneDesend.csv")

GOres2 <- readRDS("processed_data/12.12_TRAV/geneset/gsea_var_nmf/gsea_143.rds")

dotplot(GOres2,showCategory=20)

GO143 <- readRDS("processed_data/12.12_TRAV/geneset/gsea_go_var/go_143.rds")

dotplot(GO143,showCategory=20)
ggsave("result/12.26_fig4/OPC_specific_GO.pdf",width = 6,height = 6)


#== diagrame----------------------


Heatmap(filterMat[1:100,],cluster_columns = F,cluster_rows = T,show_row_names = F,
        show_column_names = F,use_raster = T,show_column_dend = F,show_row_dend = F,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))

model_test <- RcppML::nmf(t(filterMat[geneDesend[300:400],]), 3, verbose = F, seed = 1234)
w <- model_test$w
d <- model_test$d
h <- model_test$h

#colnames(h) <- rownames(filterMat)

Heatmap(t(h),show_column_names = F,show_row_names = F,cluster_columns = F,show_row_dend = F,col=colorRamp2(c(0, 0.02, 0.04), c("DeepSkyBlue3", "white", "red")))
Heatmap(valAll[150:200,],show_column_names = F,show_row_names = F,cluster_columns = F,show_row_dend = F)
Heatmap(avg_load[geneDesend[1:100],1:100],show_column_names = F,show_row_names = F,
        cluster_columns = T,show_row_dend = F,col=colorRamp2(c(0, 0.001, 0.002), c("DeepSkyBlue3", "white", "red")))
