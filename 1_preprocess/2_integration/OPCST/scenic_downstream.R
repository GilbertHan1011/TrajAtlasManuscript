#============ 
library(SCENIC)
library(SCopeLoomR)
library(scFunctions)
library(Seurat)
scenicLoomPath='../3.9_wt_integrate/9.5_scenic/9.6_auc.loom'
loom <- open_loom(scenicLoomPath)
loom
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)

wt_integrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
wt_integrate <- wt_integrate[,regulonAUC@colData@rownames]
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
cellInfo <- data.frame(seuratCluster=wt_integrate$C19_named)
cellInfoAuc <- cellInfo%>%rownames_to_column("rowname")%>%
  .[rownames(cellInfo)%in%rownames(regulonAUC@colData),]
rownames(cellInfoAuc)<- NULL
cellInfoAuc<- column_to_rownames(cellInfoAuc,"rowname")

regulonActivity_byCellType <- sapply(split(as.integer(rownames(cellInfoAuc)), cellInfoAuc$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled<- regulonActivity_byCellType_Scaled[,c(1,2,3,4,5,6)]
set.seed(1234)
ComplexHeatmap::Heatmap(na.exclude(regulonActivity_byCellType_Scaled), name="Regulon activity",row_km = 7)
col_order <- c("Common progenitor","Osteo-lineage progenitor","Osteogenic precursor","osteo-start","mine-start","Osteoblast")



auc_matrix <- as.matrix(regulonAUC@assays@data$AUC)
auc_assay <- CreateAssayObject(counts = auc_matrix)
wt_integrate[["auc"]] <- auc_assay
Idents(wt_integrate) <- wt_integrate$C19_named
DefaultAssay(wt_integrate) <- "auc"

aucMarker <- FindAllMarkers(wt_integrate,min.pct = 0.01,logfc.threshold = 0.01)
marker <- unique(aucMarker$gene)
marker <- aucMarker%>%
  filter(cluster%in%c("Lepr+ BMSC","Pre-ob","Ob"))%>%
  dplyr::select(gene)%>%unlist()%>%
  unlist()
avg <- AverageExpression(wt_integrate,return.seurat = F,assays = "auc",features = marker)
plotDataRegulon <- avg$auc
library(pheatmap)
pheatmap(plotDataRegulon,cluster_rows = T,cluster_cols = F,scale = 'row')
plotDataRegulon <- plotDataRegulon%>%t()%>%scale()%>%t()

leprLineage <- subset(wt_integrate,idents = c("Lepr+ BMSC","Pre-ob","Ob"))
marker_lepr <- aucMarker%>%
  filter(cluster%in%c("Lepr+ BMSC","Pre-ob","Ob"))%>%
  dplyr::select(gene)%>%unlist()%>%
  unlist()
avg <- AverageExpression(leprLineage,return.seurat = F,assays = "auc",features = marker_lepr)
plotDataRegulon <- avg$auc
library(pheatmap)
pheatmap(plotDataRegulon,cluster_rows = T,cluster_cols = F,scale = 'row')
plotDataRegulon <- plotDataRegulon%>%t()%>%scale()%>%t()

mesLineage <- subset(wt_integrate,idents = c("Lepr+ BMSC","Middle.MSC", "Ob", "Irx1+ MSC","Pre-ob","Late.MSC", "Early.MSC"))
marker_mes <- aucMarker%>%
  filter(cluster%in%c("Lepr+ BMSC","Middle.MSC", "Ob", "Irx1+ MSC","Pre-ob","Late.MSC", "Early.MSC"))%>%
  dplyr::select(gene)%>%unlist()%>%
  unlist()
avg <- AverageExpression(mesLineage,return.seurat = F,assays = "auc",features = marker_mes)
plotDataRegulon <- avg$auc
library(pheatmap)
plotDataRegulon <- plotDataRegulon[,c( "Irx1+ MSC","Early.MSC","Middle.MSC", "Late.MSC","Lepr+ BMSC", "Pre-ob","Ob")]
pheatmap(plotDataRegulon,cluster_rows = T,cluster_cols = F,scale = 'row')
plotDataRegulon <- plotDataRegulon%>%t()%>%scale()%>%t()



DefaultAssay(wt_integrate) <- "RNA"
mesenchymalMine<- BuildClusterTree(mesenchymalMine,features = VariableFeatures(mesenchymalMine))
PlotClusterTree(mesenchymalMine)