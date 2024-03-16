library(MuDataSeurat)
travMatrix <- read.csv("processed_data/12.12_TRAV/12.26_valide_origin_nmf.csv",row.names = 1)
RAVindex <- read.csv("processed_data/12.12_TRAV/24.1.31_ravMatrix.csv",row.names = 1)
trajMap <- ReadH5MU("../important_processed_data/12.30_trajMap.h5mu")

rownames(travMatrix) <- paste0("T",rownames(travMatrix))
leprTRAVDeg <- read.csv("processed_data/12.12_TRAV/1.31_leprTRAVDeg.csv",row.names = 1)
leprTRAVName <- leprTRAVDeg$names
leprTRAVName <- paste0("T",leprTRAVName)
conservedTRAV <- read.csv("processed_data/12.12_TRAV/conservedtravColSumDf_trav.csv",row.names = 1)

conservedTRAV <- rownames(conservedTRAV)[conservedTRAV>43]


TRAVsplit <- rep("Other",nrow(travMatrix))
TRAVsplit[rownames(travMatrix)%in%conservedTRAV] <- "osteogensis conserved"
TRAVsplit[rownames(travMatrix)%in%leprTRAVName] <- "Lepr+ BMSC specific"
levels(TRAVsplit) <- rev(c("Other", "osteogensis conserved", "Lepr+ BMSC specific"))
trajSplit <- colnames(travMatrix)%>%strsplit("_sep_")%>%lapply(`[`, 2)%>%unlist()

travMatrix_scale <- t(scale(t(travMatrix)))
pdf("result/24.2.22_fig5_supp/TRAV_split.pdf",width = 6,height = 6)
ComplexHeatmap::Heatmap(travMatrix,row_split = TRAVsplit,column_split = trajSplit,show_row_names = F,
                        show_column_names = F,border = T,cluster_row_slices = F,raster_quality = 8,use_raster = T)
dev.off()
# ComplexHeatmap::Heatmap(travMatrix,column_split = trajSplit,show_row_names = F,
#                         show_column_names = F,border = T,cluster_row_slices = F,km = 12)


colnames(RAVindex) <- rownames(travMatrix)

i="TRAV143"
namesDesend <- names(sort(travMatrix[i,],decreasing = T))
geneDesendDf <- RAVindex[,i]
names(geneDesendDf) <- rownames(RAVindex)
geneDesendDf <- geneDesendDf[geneDesendDf>0]
geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)][1:20]

# Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F)
# Heatmap(peak[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F)

filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Septoclasts_Kishor_Pdgfra.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))
haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(filterMat)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)
pdf("result/24.2.22_fig5_supp/TRAV143_Septoclasts_Kishor_Pdgfra.pdf",width = 6,height = 6)
Heatmap(filterMat,cluster_columns = F,cluster_rows = F,show_row_names = T,
        show_column_names = F,use_raster = T,top_annotation = haPseudutime,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
dev.off()

geneVal <- RAVindex[,i]
geneVal <- geneVal[geneDesend]

geneVal <- scale(geneVal)
colorpalette <- brewer.pal(3, "Blues")
row_ha = rowAnnotation(NMF =geneVal,col=list(NMF=colorRamp2(c(-2, 0,4), c("white","lightblue","blue"))))


i="TRAV77"
#namesDesend <- names(sort(travMatrix[i,],decreasing = T))
geneDesendDf <- RAVindex[,i]
names(geneDesendDf) <- rownames(RAVindex)
geneDesendDf <- geneDesendDf[geneDesendDf>0]
geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)][1:20]

filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Ablation_Matsushita_cxcl2.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))
haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(filterMat)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)
pdf("result/24.2.22_fig5_supp/TRAV77_Ablation_Matsushita_cxcl2.pdf",width = 6,height = 6)

#== Age related TRAV----------------------------------
library(lme4)
#valAllScale <- t(scale(t(valAll)))
travLepr <- travMatrix[,trajSplit=="Lepr_BMSC"]
trajMeta <- trajMap@meta.data[trajMap$Lineage=="Lepr_BMSC",]
trajMeta$Age <-  factor(trajMeta$Age,levels = c("Postnatal", 
                                                "Young Adult","Adult",  "Old"))
trajMeta$Age_num <- as.numeric(trajMeta$Age)
trajMeta <- trajMeta[colnames(travLepr),]
form <- ~ Age_num
fit <- dream(travLepr, form, trajMeta)
fit <- eBayes(fit)
AgeTb <- topTable(fit, coef = "Age_num", number = Inf)

trajMeta$ageOPC <- paste0(trajMeta$Lineage,"_",trajMeta$Age)

form2 <- ~ Age_num + (1|Lineage)
fit2 <- dream(valAllScale, form2, trajMeta)
fit2 <- eBayes(fit2)
AgeTb2 <- topTable(fit2, coef = "Age_num", number = Inf)

trajMeta$ageOPC <- paste0(trajMeta$Lineage,"_",trajMeta$Age)


AgeRelate <- AgeTb%>%filter(adj.P.Val<0.01)%>%rownames_to_column("TRAV")%>%
  dplyr::select(TRAV)%>%unlist()
intersect(AgeRelate,leprTRAVName)

sigGene <- rownames(fdrMatrix)
i="TRAV140"
#namesDesend <- names(sort(travMatrix[i,],decreasing = T))
geneDesendDf <- RAVindex[,i]
names(geneDesendDf) <- rownames(RAVindex)
geneDesendDf <- geneDesendDf[geneDesendDf>0]
geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)][1:100]
intersect(geneDesend,sigGene)

geneList <- list(`TrajDiff gene` = sigGene, 
                 `TRAV 140 genes` = geneDesend)

library(ggvenn)
pdf("result/24.2.22_fig5_supp/vennplot.pdf")
ggvenn(geneList)
dev.off()
#ggVennDiagram(geneList, label_alpha = 0,set_color = "black")+
  # ggplot2::scale_fill_gradient(low="lightblue",high = "#FF8000")
  # 
  
trajName=travMatrix["TRAV140",]%>%unlist
trajNameDescend <- sort(trajName,,decreasing = T)
leprPostnatalName <- rownames(trajMeta)[trajMeta$Age_OPC=="Lepr_BMSC_Postnatal"]
leprAdultName <- rownames(trajMeta)[trajMeta$Age_OPC=="Lepr_BMSC_Adult"]
trajNameDescend[leprPostnatalName]
trajNameDescend[leprAdultName]
i="TRAV140"
geneDesendDf <- RAVindex[,i]
names(geneDesendDf) <- rownames(RAVindex)
geneDesendDf <- geneDesendDf[geneDesendDf>0]
geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)][1:20]
filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Metaphysis_Yang_3.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))
haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(filterMat)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)

pdf("result/24.2.22_fig5_supp/TRAV140_Metaphysis_Yang_3.pdf",width = 6,height = 6)
Heatmap(filterMat,cluster_columns = F,cluster_rows = F,show_row_names = T,
        show_column_names = F,use_raster = T,top_annotation = haPseudutime,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
dev.off()

filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Bmsc2019_Regev_bm3.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))
haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(filterMat)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)
pdf("result/24.2.22_fig5_supp/TRAV140_Bmsc2019_Regev_bm3.pdf",width = 6,height = 6)
Heatmap(filterMat,cluster_columns = F,cluster_rows = F,show_row_names = T,
        show_column_names = F,use_raster = T,top_annotation = haPseudutime,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
dev.off()
