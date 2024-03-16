setwd("../3.9_wt_integrate/")
ravModel <- readRDS("processed_data/12.12_TRAV/12.26_ravmodel_nmf2.Rds")
geneList <- RAVindex(ravModel)[,143]
geneList <- sort(geneList, decreasing = TRUE)

write.csv( RAVindex(ravModel),"processed_data/12.12_TRAV/24.1.31_ravMatrix.csv")

valMat=read.csv("processed_data/12.12_TRAV/12.26_valide_origin_nmf.csv",row.names = 1)

library(MuDataSeurat)
trajMap <- ReadH5MU("../important_processed_data/12.30_trajMap.h5mu")
DefaultAssay(trajMap) <- "TRAV"
Idents(trajMap) <- trajMap$Lineage
travMarker <- FindAllMarkers(trajMap,group.by = "Lineage",test.use = "wilcox")
leprTRAVDeg <- read.csv("processed_data/12.12_TRAV/1.31_leprTRAVDeg.csv",row.names = 1)
leprTRAVName <- leprTRAVDeg$names%>%
  gsub("RAV","",.)%>%
  as.numeric()
#leprTRAVName <- 

leprGeneList <- list()
for (i in leprTRAVName){
  geneSet <- RAVindex(ravModel)[,i]
  geneSet <- sort(geneSet, decreasing = TRUE)
  geneSet <- names(geneSet)[geneSet>0]
  k = min(length(geneSet),100)
  name <- paste0("TRAV_",as.character(i))
  leprGeneList[[name]] <- geneSet[1:k]
  
}

#== select representive TRAV--------------------
goLeprTrav <- compareCluster(geneCluster = leprGeneList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
goLeprTravRes <- goLeprTrav@compareClusterResult
ckResSub <- goLeprTravRes[c("Description","Cluster","p.adjust")]


ckWide <- pivot_wider(ckResSub,names_from = Cluster,values_from = p.adjust)
ckWide <- ckWide%>%column_to_rownames("Description")
ckWideMod <- ckWide%>%
  log10()*-1
ckWideMod[is.na(ckWideMod)] <- 0
ckRowList <- list()
for (i in 1:10){
  ckRow <- ckWideMod[i]-rowMeans(ckWideMod[-1])
  ckRowList[[i]] <- ckRow
}

ckRowDf <- data.frame(do.call(cbind, ckRowList))
ckRowDf%>%
  arrange()

ckRowDf %>% arrange(desc(across(1)))%>%rownames%>%.[1:4]
leprSelectGO <- c()
for (i in 1:10){
  leprSelectGO <- c(leprSelectGO,ckRowDf %>% arrange(desc(across(i)))%>%rownames%>%.[1])
}

ckResSelect <- goLeprTravRes%>%filter(Description%in% leprSelectGO)
ckResSelect$p.adjust <- -log10(ckResSelect$p.adjust)

ckResSelect$Description <- factor(ckResSelect$Description,levels = leprSelectGO)
ckResSelect <- dplyr::arrange(ckResSelect,Description)
goLeprTrav@compareClusterResult <- ckResSelect


dotplot(goLeprTrav,)+  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 15,face = "bold"),  # Adjust text size here
  axis.text.y = element_text(size = 15,face = "bold"),  # Adjust y-axis text size
  axis.title = element_text(size = 14),  # Adjust axis title size
  legend.text = element_text(size = 12),  # Adjust legend text size
  legend.title = element_text(size = 14) # Adjust legend title size
)+scale_fill_gradient( low = "DeepSkyblue3", high = "red")
ggsave("result/24.1.31_fig5_trajmap_plot/2.1_dotPlot.pdf",width = 8,height = 6)

ckResSelect$GeneRatio <- sapply(strsplit(ckResSelect$GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

ggplot(ckResSelect, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), high = "red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")

saveRDS(goLeprTrav,"processed_data/12.12_TRAV/1.31_GOLeprTrav.Rds")
write.csv(goLeprTravRes,"processed_data/12.12_TRAV/1.31_GOLeprTable.csv")

test <- ckRowList[[1]]
sort(ckRowList[[1]])


ckWide <- ckWide[,c("Description", "Up_Up", "Up_Down", "Up_0", "Down_Down", "Down_0", "0_Down", 
                    "Down_Up")]

geneList <- RAVindex(RAVmodel)[,i]
geneList <- sort(geneList, decreasing = TRUE)
geneList <- names(geneList)[geneList>0]
k = min(length(geneList),100)





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
