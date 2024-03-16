#== plot mes injury trajDiff----------------------------

#== Run GO-----------------------------------
library(clusterProfiler)
library(ComplexHeatmap)
library(ReactomePA)
library(simplifyEnrichment)
library(circlize)
setwd("../8.25_full_integration/")
exprMatrix <- data.table::fread("process_data/traj_diff/exprMatrix.csv",header = T)
exprMatrix <- exprMatrix%>%column_to_rownames("V1")
fdrMatrix <-  read.csv("process_data/traj_diff//fdrMatrix.csv",row.names = 1)
cpmMatrix <-  read.csv("process_data/traj_diff/cpmMatrix.csv",row.names = 1)
cpm1 <- cpmMatrix[1:94]
cpm2 <- cpmMatrix[95:188]

fdrBinaray <- fdrMatrix>1
exprBinary <- exprMatrix*fdrBinaray

start <- exprBinary[,1:31]
middle <-  exprBinary[,32:62]
end <-  exprBinary[,63:94]
startSum <- rowSums(start)
middleSum <- rowSums(middle)
endSum <- rowSums(end) 
  
#==to define different category genes
threshold <- 10
geneWholeUp <- names(startSum)[startSum>threshold&endSum>threshold]
geneStartUp <- names(startSum)[startSum>threshold&abs(endSum)<threshold]
geneEndUp <- names(startSum)[abs(startSum)<threshold&endSum>threshold]
geneWholeDown <- names(startSum)[startSum< -threshold&endSum< -threshold]
geneStartdown <- names(startSum)[startSum< -threshold& abs(endSum)< threshold]
geneEndDown <- names(startSum)[abs(startSum)< threshold&endSum< -threshold]
geneUpDown <- names(startSum)[startSum> threshold& endSum< -threshold]
geneDownUp <- names(startSum)[startSum< -threshold&endSum> threshold]

geneList <- list(geneWholeUp,geneStartUp,geneEndUp,geneWholeDown,geneStartdown,geneEndDown,geneUpDown ,geneDownUp)
names(geneList) <- c("Up_Up","Up_0","0_Up","Down_Down","Down_0","0_Down","Up_Down","Down_Up")
GOlist <- lapply(geneList,FUN = function(x){
  enrichGO(gene = x,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
})

GO_ck <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
dotplot(GO_ck)

ckRes <- GO_ck@compareClusterResult
ckResSub <- ckRes[c("Description","Cluster","p.adjust")]

ckWide <- pivot_wider(ckResSub,names_from = Cluster,values_from = p.adjust)

ckWide <- ckWide[,c("Description", "Up_Up", "Up_Down", "Up_0", "Down_Down", "Down_0", "0_Down", 
                   "Down_Up")]
geneMiddleUp <- names(startSum)[middleSum> threshold&abs(startSum)< threshold&abs(endSum)< threshold]
geneMiddleDown <- names(startSum)[middleSum< -threshold&abs(startSum)< threshold&abs(endSum)< threshold]

write.csv(ckWide,"process_data/traj_diff/1.30_ck_res.csv")
write.csv(ckWide,"../8.25_full_integration/process_data//traj_diff/1.30_ck_res.csv")
saveRDS(GO_ck,"../8.25_full_integration/process_data/traj_diff/")
#== select representive GO term------------------

GO_term <- read.table("process_data/traj_diff/GO_select",sep = "\n")%>%unlist()
ckResSelect <- ckRes[ckRes$Description%in%GO_term,]

GO_ck@compareClusterResult <- ckResSelect
dotplot(GO_ck)+  theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 15,face = "bold"),  # Adjust text size here
  axis.text.y = element_text(size = 15,face = "bold"),  # Adjust y-axis text size
  axis.title = element_text(size = 14),  # Adjust axis title size
  legend.text = element_text(size = 12),  # Adjust legend text size
  legend.title = element_text(size = 14)  # Adjust legend title size
)
ggsave("result/1.1_fig6/1.30_dotplot_ck.pdf",width = 8,height = 10)


#== heatmap-----------------------

geneGroup <-  c(geneWholeUp,geneStartUp,geneEndUp,geneWholeDown,geneStartdown,geneEndDown,geneUpDown ,geneDownUp)
geneCat <- c(
  rep(names(geneList)[[1]],length(geneList[[1]])),
  rep(names(geneList)[[2]],length(geneList[[2]])),
  rep(names(geneList)[[3]],length(geneList[[3]])),
  rep(names(geneList)[[4]],length(geneList[[4]])),
  rep(names(geneList)[[5]],length(geneList[[5]])),
  rep(names(geneList)[[6]],length(geneList[[6]])),
  rep(names(geneList)[[7]],length(geneList[[7]])),
  rep(names(geneList)[[8]],length(geneList[[8]]))
)
geneSplit <- data.frame("gene"=geneGroup,"cat"=geneCat)

geneSplit$cat <- factor(geneSplit$cat,levels = c( "Up_Up", "Up_0", "0_Up", "Up_Down","Down_Down", "Down_0", "0_Down", 
                        "Down_Up"))

hm1 <- Heatmap(exprMatrix[geneSplit$gene,], name = "logChange",
        col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, 
         column_title = NULL,cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE,split = geneSplit$cat,
        border =T, cluster_row_slices = F, row_gap = unit(2, "mm"), row_title_rot = 0,
        column_title = "Diff_Expression",raster_quality = 5)


hm_cpm1 <- Heatmap(cpm1[geneSplit$gene,], name = "cpm1",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, 
                   cluster_rows = F,cluster_columns = F,
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   split = geneSplit$cat,border =T, cluster_row_slices = F, 
                   row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "MSC(Injury)",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
#hm_cpm1 <- draw(hm_cpm1)

hm_cpm2 <- Heatmap(cpm2[geneSplit$gene,], name = "cpm2",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, 
                   cluster_rows = F,cluster_columns = F,row_title = NULL,
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   split = geneSplit$cat,border =T, cluster_row_slices = F, 
                   row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "MSC(UnInjury)",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
#hm_cpm2 <- draw(hm_cpm2)

hm_expr <- Heatmap(exprMatrix[geneSplit$gene,], name = "logChange",
                   col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),
                   show_row_names = FALSE, show_column_names = FALSE, cluster_rows = T,cluster_columns = F,
                   show_row_dend = FALSE, show_column_dend = FALSE,split = geneSplit$cat,
                   border =T, cluster_row_slices = F, row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "Diff_Expression",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
#hm_expr <- draw(hm_expr)

hm_fdr <- Heatmap(fdrMatrix[geneSplit$gene,], name = "fdr",
                  col = colorRamp2(c(0, 10, 20),hcl_palette="Spectral",reverse = TRUE),
                  show_row_names = FALSE, show_column_names = FALSE, cluster_rows = T,cluster_columns = F,
                  show_row_dend = FALSE, show_column_dend = FALSE,split = geneSplit$cat,
                  border =T, cluster_row_slices = F, row_gap = unit(2, "mm"),column_title = "FDR",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)

hmList <- hm_cpm1+hm_cpm2+hm_expr+hm_fdr
hmList
pdf("result/1.1_fig6/1.30_geneDiff4Panel.pdf",width = 8,height = 10)
draw(hmList)
dev.off()

pseudotime <- 1:ncol(fdrMatrix)
colFun <- colorRamp2(c(0, ncol(fdrMatrix)/2, ncol(fdrMatrix)),hcl_palette="Spectral",reverse = TRUE)
ha = HeatmapAnnotation(pseudotime=pseudotime, col=list(psedotime=colFun))
haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(fdrMatrix)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Reds",reverse = T)
  )
)
Heatmap(fdrMatrix[geneSplit$gene,], name = "fdr",
        col = colorRamp2(c(0, 10, 20),hcl_palette="Spectral",reverse = TRUE),
        show_row_names = FALSE, show_column_names = FALSE, cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE,split = geneSplit$cat,
        border =T, cluster_row_slices = F, row_gap = unit(2, "mm"),column_title = "FDR",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)

Heatmap(exprMatrix[geneStartUp,], name = "logChange",
        col = colorRamp2(c(-10, 0, 10), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, 
        row_title = NULL, column_title = NULL,cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(exprMatrix[geneEndUp,], name = "logChange",
        col = colorRamp2(c(-3, 0, 3), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, 
        row_title = NULL, column_title = NULL,cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(exprMatrix[geneEndDown,], name = "logChange",
        col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, 
        row_title = NULL, column_title = NULL,cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(exprMatrix[geneWholeDown,], name = "logChange",
        col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, 
        row_title = NULL, column_title = NULL,cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE)
Heatmap(exprMatrix[geneStartdown,], name = "logChange",
        col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, 
        row_title = NULL, column_title = NULL,cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE)



Heatmap(fdrMatrix[geneStartUp,], name = "logChange",
        col = colorRamp2(c(-10, 0, 10), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, 
        row_title = NULL, column_title = NULL,cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE)




#exprMatrix <- exprMatrix[, -ncol(exprMatrix)]
label <- read.csv("process_data/traj_diff/1.6_label.csv")
kmean=label$X0
names(kmean) <- label$X
kmean <- kmean[rownames(exprMatrix)]
geneKmeans <- split(label,label$X0)
geneKmeans <- lapply(geneKmeans,`[[`,1)

GOlist <- lapply(geneKmeans,FUN = function(x){
  enrichGO(gene = x,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
})

geneIDs <- lapply(geneKmeans,bitr,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
geneIDs <- lapply(geneIDs,`[[`,2)
Reactomelist <- lapply(geneIDs,FUN = function(x){
  enrichPathway(gene = x,organism = "mouse")
})
dotplot(Reactomelist[[1]])

pipeFilter <- function(x){
  x <- x@result
  go <- x%>%filter(p.adjust<0.05)%>%
    dplyr::select(ID)%>%
    unlist
  return(go)
}
goDesList <-  lapply(GOlist, pipeFilter)
names(goDesList) <- 0:9
lenGO <- lapply(goDesList,length)%>%unlist()
goDesList <- goDesList[lenGO!=0]
ha_GO= rowAnnotation(go = anno_word_cloud_from_GO(kmean, goDesList, max_words = 30))
ha_GO2= rowAnnotation(go = anno_word_cloud_from_GO(kmean, goDesList, max_words = 15))



