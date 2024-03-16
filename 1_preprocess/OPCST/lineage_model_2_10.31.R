#============ 
library(SCENIC)
library(SCopeLoomR)
library(scFunctions)
library(Seurat)
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
lfile <- connect(filename = "../3.9_wt_integrate/9.5_scenic/9.6_wt_full_prepare.loom", mode = "r+",skip.validate = T)
lfile

scenicLoomPath='../3.9_wt_integrate/9.5_scenic/10.2_dpt/10.30_auc.loom'
dpt_bin <- read.csv("../important_processed_data/11.4_dpt_bin_df.csv")
loom <- open_loom(scenicLoomPath)
loom
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)

#wt_integrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
dpt_sub=zellkonverter::readH5AD("../temp_data/10.27_dpt_sub.h5ad")
wt_integrate <- as.Seurat(dpt_sub,counts = NULL,data = "X")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

dpt_bin=read.csv("../temp_data/10.31_dpt_bin_16w.csv")
wt_integrate$dpt_bin=dpt_bin
c("chondro7.0", "chondro6.0", "chondro8.0", "chondro4.0", "chondro5.0", 
  "chondro3.0", "chondro2.0", "chondro10.0", "chondro1.0", "chondro11.0", 
  "chondro9.0", "Fibro9.0", "Fibro6.0", "Fibro8.0", "Fibro7.0", 
  "Fibro4.0", "Fibro5.0", "Fibro10.0", "Fibro2.0", "Fibro3.0", 
  "Fibro1.0", "Lepr9.0", "Lepr2.0", "Lepr10.0", "Lepr3.0", "Lepr6.0", 
  "Lepr7.0", "Lepr8.0", "Lepr5.0", "Lepr4.0", "Lepr1.0", "Lepr11.0", 
  "Mesenchyme5.0", "Mesenchyme8.0", "Mesenchyme7.0", "Mesenchyme4.0", 
  "Mesenchyme6.0", "Mesenchyme3.0", "Mesenchyme9.0", "Mesenchyme10.0", 
  "Mesenchyme2.0", "Mesenchyme1.0")


write.csv(as.data.frame(unique(dpt_bin$dpt_lineage)),"../temp_data/11.4_dpt_lineage_temp.csv")

dpt_dict <- read.csv("../temp_data/11.4_dpt_lineage_temp.csv",row.names = 1)

dpt_bin <- dpt_bin %>%
  left_join(dpt_dict, by = "dpt_lineage")

hmcell <- dpt_bin$hmLineage
names(hmcell) <- dpt_bin$X
wt_integrate$hmLineage <- hmcell[colnames(wt_integrate)]



leprCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Lepr1.0","Lepr2.0","Lepr3.0","Lepr4.0","Lepr5.0")]%>%unique()
chondroCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("chondro7.0", "chondro6.0", "chondro4.0", "chondro5.0", 
                                                  "chondro3.0", "chondro2.0", "chondro1.0")]%>%unique()
MesCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Mesenchyme2.0", "Mesenchyme1.0","Mesenchyme3.0","Mesenchyme4.0")]%>%unique()
FibroCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro2.0", "Fibro3.0", "Fibro1.0" , "Fibro4.0" )]%>%unique()
Fibro_Mes <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro5.0", "Mesenchyme5.0")]%>%unique()
Fibro_Mes_Lepr <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro6.0", "Mesenchyme6.0","Lepr6.0","Fibro7.0", "Mesenchyme7.0","Lepr7.0")]
Osteo <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro8.0", "Mesenchyme8.0","Lepr8.0","chondro8.0",
                                            "Fibro9.0", "Mesenchyme9.0","Lepr9.0","chondro9.0",
                                            "Fibro10.0", "Mesenchyme10.0","Lepr10.0","chondro10.0",
                                            "Fibro11.0","Lepr11.0","chondro11.0")]%>%unique()

wt_integrate$coverage_node <- "None"
wt_integrate$coverage_node[chondroCell] <- "Chondro"
wt_integrate$coverage_node[leprCell] <- "Lepr"
wt_integrate$coverage_node[MesCell] <- "Mesenchyme"
wt_integrate$coverage_node[FibroCell] <- "Fibro"
wt_integrate$coverage_node[Fibro_Mes] <- "Fibro_Mes"
wt_integrate$coverage_node[Fibro_Mes_Lepr] <- "Fibro_Mes_Lepr"
wt_integrate$coverage_node[Osteo] <- "Osteo"

dpt_bin_val <- dpt_bin$dpt_lineage
names(dpt_bin_val) <- dpt_bin$X
dptCell <-  dpt_bin$X
names(dptCell) <- dpt_bin$dpt_lineage
wt_integrate$dpt_bin_lineage <- "None"
wt_integrate$dpt_bin_lineage[wt_integrate$coverage_node!="Osteo"]=dpt_bin_val[colnames(wt_integrate)[wt_integrate$coverage_node!="Osteo"]]
osteo

wt_integrate$dpt_bin_lineage[dptCell[]]

Idents(wt_integrate) <- wt_integrate$hmLineage
aucMarker <- FindAllMarkers(wt_integrate,min.pct = 0.01,logfc.threshold = 0.01,only.pos = T)
aucGene <- aucMarker%>%filter(avg_log2FC>0)%>%group_by(cluster)%>%top_n(20)
marker <- unique(aucGene$gene)
# marker <- aucMarker%>%
#   filter(cluster%in%c("Lepr+ BMSC","Pre-ob","Ob"))%>%
#   dplyr::select(gene)%>%unlist()%>%
#   unlist()
avg <- AverageExpression(wt_integrate,return.seurat = F,assays = "auc",features = marker)
plotDataRegulon <- avg$auc

pheatmap(plotDataRegulon,cluster_rows = T,cluster_cols = F,scale = 'row')
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

plotDataRegulon <- scale_rows(plotDataRegulon)
library(ComplexHeatmap)
set.seed(12345)
hm <- Heatmap(plotDataRegulon,cluster_columns = T,cluster_rows = T,km = 10,show_row_names = F)
pdf("result/10.30_lineage_draw/11.4_km_parttern.pdf",width = 10,height = 10)
hm <- draw(hm)
hm
dev.off()
saveRDS(hm,"result/10.30_lineage_draw/11.14_regulonHM.Rds")

clusterlist = row_order(hm)
countData <- hm@ht_list$matrix_12@matrix
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(countData[clusterlist[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)
rownames(clu_df) <- clu_df$GeneID
write.csv(clu_df,"result/10.30_lineage_draw/hm_cluster.csv")

#regulonList <- clu_df%>%group_by(Cluster)
split_df <- split(clu_df$GeneID, clu_df$Cluster)

#wt_integrate <- AddModuleScore(wt_integrate,features = split_df,name = names(split_df),ctrl = 10)
clusterGene <- list("earlyMesenchyme"=c(split_df$cluster2),"lateMesenchyme"=split_df$cluster1,
                    "earlyChondro"=split_df$cluster7,
                    "lateChondro_1"=c(split_df$cluster9),"lateChondro_2"=c(split_df$cluster10),"Osteo"=split_df$cluster8,
                    "Lepr"=split_df$cluster4,"Fibro"=split_df$cluster6)
meta <- wt_integrate@meta.data
meta <- meta %>%
  select(-ncol(meta):-(ncol(meta)-11))
wt_integrate@meta.data <- meta
wt_integrate <- AddModuleScore(wt_integrate,features = clusterGene,name = names(clusterGene),ctrl = 10)


plotData=wt_integrate@meta.data[c("hmLineage","Mesenchyme1", 
                                  "earlyChondro2", "lateChondro_13", "lateChondro_24", "Osteo5", 
                                  "Lepr6", "Fibro7")]

columns_to_mean <- c("Mesenchyme1","Mesenchyme1", 
                     "earlyChondro2", "lateChondro_13", "lateChondro_24", "Osteo5", 
                     "Lepr6", "Fibro7")  # Add the column names you want to calculate the mean for

# Use the group_by and summarize_at functions to calculate the mean for selected columns
result <- plotData %>%
  group_by(hmLineage) %>%
  summarize_at(vars(columns_to_mean), ~ mean(., na.rm = TRUE))
result <- result%>%column_to_rownames("hmLineage")
FeaturePlot(wt_integrate,"lateChondro3",raster = T)
pdf("result/10.30_lineage_draw/score_df.pdf",width = 6,height = 8)
write.csv(result,"result/10.30_lineage_draw/11.4_avg_score_hm_base.csv")

Heatmap(result)
dev.off()
# 
# cellInfo <- data.frame(seuratCluster=wt_integrate$coverage_node)
# cellInfoAuc <- cellInfo%>%rownames_to_column("rowname")%>%
#   .[rownames(cellInfo)%in%rownames(regulonAUC@colData),]
# rownames(cellInfoAuc)<- NULL
# cellInfoAuc<- column_to_rownames(cellInfoAuc,"rowname")
# 
# regulonActivity_byCellType <- sapply(split(as.integer(rownames(cellInfoAuc)), cellInfoAuc$seuratCluster),
#                                      function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
# regulonActivity_byCellType_Scaled<- regulonActivity_byCellType_Scaled[,c(1,2,3,4,5,6)]
# set.seed(1234)
# ComplexHeatmap::Heatmap(na.exclude(regulonActivity_byCellType_Scaled), name="Regulon activity",row_km = 7)
# col_order <- c("Common progenitor","Osteo-lineage progenitor","Osteogenic precursor","osteo-start","mine-start","Osteoblast")
# 


auc_matrix <- as.matrix(regulonAUC@assays@data$AUC)
auc_matrix <- auc_matrix[,colnames(wt_integrate)]
auc_assay <- CreateAssayObject(counts = auc_matrix)
wt_integrate[["auc"]] <- auc_assay
Idents(wt_integrate) <- wt_integrate$coverage_node
DefaultAssay(wt_integrate) <- "auc"

aucMarker <- FindAllMarkers(wt_integrate,min.pct = 0.01,logfc.threshold = 0.01)
aucGene <- aucMarker%>%filter(avg_log2FC>0)%>%group_by(cluster)%>%top_n(20)
marker <- unique(aucMarker$gene)
# marker <- aucMarker%>%
#   filter(cluster%in%c("Lepr+ BMSC","Pre-ob","Ob"))%>%
#   dplyr::select(gene)%>%unlist()%>%
#   unlist()
avg <- AverageExpression(wt_integrate,return.seurat = F,assays = "auc",features = marker)
plotDataRegulon <- avg$auc

library(pheatmap)
pheatmap(plotDataRegulon,cluster_rows = T,cluster_cols = F,scale = 'row')


columnOrder <-  c( "Fibro","Mesenchyme","Fibro_Mes", "Lepr","Fibro_Mes_Lepr", "Chondro", "Osteo")
pheatmap(plotDataRegulon[aucGene$gene,columnOrder],cluster_rows = F,cluster_cols = F,scale = 'row')
plotDataRegulon <- plotDataRegulon%>%t()%>%scale()%>%t()


aucAssay=wt_integrate@assays$auc
saveRDS(aucAssay,"../important_processed_data/11.6_scenic_aucassay.Rds")
saveRDS(hm,"result/10.30_lineage_draw/11.6_lineage_hm.Rds")
rm(list=ls())
