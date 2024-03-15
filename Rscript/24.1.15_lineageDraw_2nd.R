hmCluster <- read.csv("result/10.30_lineage_draw/hm_cluster.csv")
hm <- readRDS("result/10.30_lineage_draw/11.6_lineage_hm.Rds")



plotDatahm <- hm@ht_list$matrix_12@matrix
resultDf <- plotDatahm %>%
  group_by(hmLineage) %>%
  summarize_at(vars(columns_to_mean), ~ mean(., na.rm = TRUE))
result <- result%>%column_to_rownames("hmLineage")

regulonHM <- readRDS("result/10.30_lineage_draw/11.14_regulonHM.Rds")


aucAssay <- readRDS("../important_processed_data/11.6_scenic_aucassay.Rds")
aucDpt <- CreateSeuratObject(aucAssay)



dpt_bin <- read.csv("../important_processed_data/11.4_dpt_bin_df.csv")
dpt_dict <- read.csv("../unimportant_processed_data/11.4_dpt_lineage_temp.csv",row.names = 1)

dpt_bin <- dpt_bin %>%
  left_join(dpt_dict, by = "dpt_lineage")

hmcell <- dpt_bin$hmLineage
names(hmcell) <- dpt_bin$X
aucDpt$hmLineage <- hmcell[colnames(aucDpt)]




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

aucDpt$coverage_node <- "None"
aucDpt$coverage_node[chondroCell] <- "Chondro"
aucDpt$coverage_node[leprCell] <- "Lepr"
aucDpt$coverage_node[MesCell] <- "Mesenchyme"
aucDpt$coverage_node[FibroCell] <- "Fibro"
aucDpt$coverage_node[Fibro_Mes] <- "Fibro_Mes"
aucDpt$coverage_node[Fibro_Mes_Lepr] <- "Fibro_Mes_Lepr"
aucDpt$coverage_node[Osteo] <- "Osteo"

dpt_bin_val <- dpt_bin$dpt_lineage
names(dpt_bin_val) <- dpt_bin$X
dptCell <-  dpt_bin$X
names(dptCell) <- dpt_bin$dpt_lineage
aucDpt$dpt_bin_lineage <- "None"
aucDpt$dpt_bin_lineage[aucDpt$coverage_node!="Osteo"]=dpt_bin_val[colnames(aucDpt)[aucDpt$coverage_node!="Osteo"]]


Idents(aucDpt) <- aucDpt$coverage_node

aucMarker <- FindAllMarkers(aucDpt,min.pct = 0.01,logfc.threshold = 0.01)
aucGene <- aucMarker%>%filter(avg_log2FC>0)%>%group_by(cluster)%>%top_n(20)
marker <- unique(aucMarker$gene)
# marker <- aucMarker%>%
#   filter(cluster%in%c("Lepr+ BMSC","Pre-ob","Ob"))%>%
#   dplyr::select(gene)%>%unlist()%>%
#   unlist()
avg <- AverageExpression(aucDpt,return.seurat = F,features = marker)
plotDataRegulon <- avg$RNA%>%as.data.frame()
pheatmap(plotDataRegulon,cluster_rows = T,cluster_cols = F,scale = 'row')
columnOrder <-  c( "Fibro","Mesenchyme","Fibro-Mes", "Lepr","Fibro-Mes-Lepr", "Chondro", "Osteo")
pheatmap(plotDataRegulon[aucGene$gene,columnOrder],cluster_rows = F,cluster_cols = F,scale = 'row')



#== score--------------


clusterlist = row_order(regulonHM)
countData <- regulonHM@ht_list$matrix_1@matrix
clu_df2 <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(countData[clusterlist[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)
rownames(clu_df2) <- clu_df2$GeneID

clu_df <- read.csv("result/10.30_lineage_draw/hm_cluster.csv")

#regulonList <- clu_df%>%group_by(Cluster)
split_df2 <- split(clu_df2$GeneID, clu_df2$Cluster)
split_df2 <- split(hmCluster$GeneID, hmCluster$Cluster)
#wt_integrate <- AddModuleScore(wt_integrate,features = split_df,name = names(split_df),ctrl = 10)
clusterGene <- list("Mesenchyme"=c(split_df$cluster1),
                    "earlyChondro"=split_df$cluster2,
                    "lateChondro"=c(split_df$cluster6),"Osteo"=split_df$cluster4,
                    "Lepr"=split_df$cluster8,"Fibro"=split_df$cluster3)
meta <- aucDpt@meta.data
meta <- meta %>%
  select(-ncol(meta):-(ncol(meta)-7))
aucDpt@meta.data <- meta

aucDpt <- AddModuleScore(aucDpt,features = clusterGene,name = names(clusterGene),ctrl = 10)


plotData=aucDpt@meta.data[c("hmLineage","Mesenchyme1", "earlyChondro2", "lateChondro3", 
                            "Osteo4", "Lepr5", "Fibro6")]

columns_to_mean <- c("Mesenchyme1", "earlyChondro2", "lateChondro3", 
                     "Osteo4", "Lepr5", "Fibro6")  # Add the column names you want to calculate the mean for



#== create heatmap
# Use the group_by and summarize_at functions to calculate the mean for selected columns
result <- plotData %>%
  group_by(hmLineage) %>%
  summarize_at(vars(columns_to_mean), ~ mean(., na.rm = TRUE))
result <- result%>%column_to_rownames("hmLineage")
result <- t(result)
colName <- colnames(result)
binPart <- gsub("[^0-9.]", "" , colName)
binPart <- as.numeric(binPart)
lineagePart <- gsub("([0-9.]+)", "" , colName)
FibroLineageLogic <-c( grep("Fibro",lineagePart), grep("osteo",lineagePart))
LeprLineageLogic <-c( grep("Lepr",lineagePart), grep("osteo",lineagePart))
MesLineageLogic <-c( grep("Mes",lineagePart), grep("osteo",lineagePart))
chondroLineageLogic <-c( grep("chondro",lineagePart), grep("osteo",lineagePart))
Logic=rep("FALSE", length(lineagePart))
FibroLog <- LeprLog <- MesLog <- chondroLog <- Logic
FibroLog[FibroLineageLogic]="TRUE"
LeprLog[LeprLineageLogic]="TRUE"
MesLog[MesLineageLogic]="TRUE"
chondroLog[chondroLineageLogic]="TRUE"

logicColor <- c("darkgreen","orange")

names(logicColor) <- unique(chondroLog)
timeColor <- colorRamp2(c(0, 5, 10),hcl_palette="Spectral",reverse = TRUE)
timeColor <- colorRamp2(seq(from=0,to=10,length.out=20),viridis::viridis(20,option="H"),reverse = TRUE)

viridisCol <- viridis::viridis(20,option="H")

ha = HeatmapAnnotation(
  Fibro = FibroLog,Mesenchyme=MesLog,Lepr =LeprLog, Chondro=chondroLog,Pseudotime=binPart,
  col = list(
    Fibro=logicColor,Chondro=logicColor,Lepr=logicColor,Mesenchyme=logicColor,
    Pseudotime=timeColor
  )
)


topRegulon <- lapply(split_df2, head,3)
topRegulon <- topRegulon[c("cluster1", "cluster2", "cluster6",  "cluster4", "cluster8",  "cluster3")]
names(topRegulon) <-  rownames(result)

topRegulon$Mesenchyme1 <- c("Twist1(+)","Hoxa11(+)","Foxp2(+)")
topRegulon$earlyChondro2 <- c("Sox9(+)","Sox8(+)","Hoxc8(+)")
topRegulon$lateChondro3 <- c("Sox4(+)","Tead1(+)","Tgif1(+)")
topRegulon$Osteo4 <- c("Sp7(+)","Runx2(+)","Xbp1(+)")
topRegulon$Lepr5 <- c("Foxc2(+)","Tbx1(+)","Maf(+)")
topRegulon$Fibro6 <- c("Scx(+)","Sox7(+)","Pax1(+)")

rightAnno <-  rowAnnotation(
  textbox = anno_textbox(
    rownames(result), topRegulon, 
    word_wrap = TRUE, 
    add_new_line = F,
    max_width = unit(300, "mm"),side = "left",
    background_gp = gpar(fill = "grey90", col = "#808080"))
)




rowOrderHM <- row_order(hmScore)
colOrderHM <- column_order(hmScore)
colnames(result)[colOrderHM]

rowOrder <- c( "earlyChondro2", "lateChondro3",  "Lepr5","Fibro6", "Mesenchyme1", "Osteo4")
colOrder <- c( "chondro1.0", "chondro2.0","chondro3.0","chondro4.0", "chondro5.0",  "chondro6.0", 
               "chondro7.0", "chondro8.0","chondro9.0", 
               "Lepr1.0", "Lepr2.0","Lepr3.0","Lepr4.0", 
               "Fibro1.0",  "Fibro2.0","Fibro3.0",
               "Mesenchyme1.0","Mesenchyme2.0", "Mesenchyme3.0", "Mes_Fibro_4.0", 
              "Mes_Fibro_Lepr_5.0", "Mes_Fibro_Lepr_6.0" ,  "Mes_Fibro_Lepr_7.0","Mes_Fibro_Lepr_8.0",
              "Mes_Fibro_Lepr_9.0", "osteo10.0")
hmScore <- Heatmap(result,col=colorRamp2(c(-0.1, 0, 0.1), c("Deepskyblue3", "white", "red")),
                   top_annotation = ha,left_annotation  = rightAnno,show_row_dend = F,
                   row_order = rowOrder,column_order = colOrder,show_column_names = F)
hmScore <- draw(hmScore)
pdf("result/10.30_lineage_draw/1.15_scoreDf_3.pdf",width = 10,height = 3)
hmScore
dev.off()
pdf("result/10.30_lineage_draw/1.15_score_df.pdf")
hmScore <- Heatmap(result,col=colorRamp2(c(-0.1, 0, 0.1), c("Deepskyblue3", "white", "red")))
hmScore <- draw(hmScore)
dev.off()

saveRDS(hmScore,"result/10.30_lineage_draw/24.1.16_hmScore.Rds")

FeaturePlot(wt_integrate,"lateChondro3",raster = T)
pdf("result/10.30_lineage_draw/score_df.pdf",width = 6,height = 8)
write.csv(result,"result/10.30_lineage_draw/1.15_avg_score_hm_base.csv")
