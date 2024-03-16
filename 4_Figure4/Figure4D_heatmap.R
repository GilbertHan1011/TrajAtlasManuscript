#== fig4_replot-------------------------------
#2.2

exprMatrix <- data.table::fread("processed_data/12.6_leprDev_vs_Steady/exprMatrix.csv")
exprMatrix <- exprMatrix%>%column_to_rownames("V1")
exprMatrix <- exprMatrix[, -ncol(exprMatrix)]
fdrMatrix <-  read.csv("processed_data/12.6_leprDev_vs_Steady/fdrMatrix.csv",row.names = 1)
cpmMatrix <-  read.csv("processed_data/12.6_leprDev_vs_Steady/cpmMatrix.csv",row.names = 1)
cpm1 <- cpmMatrix[1:99]
cpm2 <- cpmMatrix[100:198]


fdrBinaray <- fdrMatrix>1
exprBinary <- exprMatrix*fdrBinaray

start <- exprBinary[,1:33]
middle <-  exprBinary[,34:66]
end <-  exprBinary[,67:99]
startSum <- rowSums(start)
middleSum <- rowSums(middle)
endSum <- rowSums(end) 

#==to define different category genes-------------------
threshold <- 10
geneWholeUp <- names(startSum)[startSum>threshold&endSum>threshold]
geneStartUp <- names(startSum)[startSum>threshold&abs(endSum)<threshold]
geneEndUp <- names(startSum)[abs(startSum)<threshold&endSum>threshold]
geneWholeDown <- names(startSum)[startSum< -threshold&endSum< -threshold]
geneStartdown <- names(startSum)[startSum< -threshold& abs(endSum)< threshold]
geneEndDown <- names(startSum)[abs(startSum)< threshold&endSum< -threshold]
geneUpDown <- names(startSum)[startSum> threshold& endSum< -threshold]
geneDownUp <- names(startSum)[startSum< -threshold&endSum> threshold]
geneMiddleUp <- names(startSum)[middleSum> threshold&abs(startSum)< threshold&abs(endSum)< threshold]
geneMiddleDown <- names(startSum)[middleSum< -threshold&abs(startSum)< threshold&abs(endSum)< threshold]



geneList <- list(geneWholeUp,geneStartUp,geneEndUp,geneWholeDown,geneStartdown,geneEndDown,geneUpDown ,geneDownUp)
names(geneList) <- c("Up_Up","Up_0","0_Up","Down_Down","Down_0","0_Down","Up_Down","Down_Up")
GOlist <- lapply(geneList,FUN = function(x){
  enrichGO(gene = x,OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL")
})
GO_ck <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")



ckRes <- GO_ck@compareClusterResult
ckResSub <- ckRes[c("Description","Cluster","p.adjust")]

ckWide <- pivot_wider(ckResSub,names_from = Cluster,values_from = p.adjust)

ckWide <- ckWide[,c("Description", "Up_Up","Up_0","0_Up","Down_Down","Down_0","0_Down","Up_Down","Down_Up")]


ckWide <- ckWide%>%column_to_rownames("Description")
ckWideMod <- ckWide%>%
  log10()*-1
ckWideMod[is.na(ckWideMod)] <- 0
ckRowList <- list()
for (i in 1:8){
  ckRow <- ckWideMod[i]-rowMeans(ckWideMod[-1])
  ckRowList[[i]] <- ckRow
}

ckRowDf <- data.frame(do.call(cbind, ckRowList))

leprSelectGO <- c()
for (i in 1:8){
  leprSelectGO <- c(leprSelectGO,ckRowDf %>% arrange(desc(across(i)))%>%rownames%>%.[1:2])
}

ckResSelect <- ckRes%>%filter(Description%in% leprSelectGO)
ckResSelect$p.adjust <- -log10(ckResSelect$p.adjust)

ckResSelect$Description <- factor(ckResSelect$Description,levels = leprSelectGO)
ckResSelect <- dplyr::arrange(ckResSelect,Description)
GO_ck@compareClusterResult <- ckResSelect




dotplot(GO_ck,)+  theme(
  axis.text.x = element_text(angle = 90, hjust = 1, size = 15,face = "bold"),  # Adjust text size here
  axis.text.y = element_text(size = 15,face = "bold"),  # Adjust y-axis text size
  axis.title = element_text(size = 14),  # Adjust axis title size
  legend.text = element_text(size = 12),  # Adjust legend text size
  legend.title = element_text(size = 14) # Adjust legend title size
)+scale_fill_gradient( low = "DeepSkyblue3", high = "red")

ggsave("result/24.2.2_fig4_replot/go_dotplot.pdf",width = 10,height = 10)
write.csv(ckRes, "processed_data/24.2.2_fig4_replot/GO_ck.csv")
saveRDS(GO_ck,"processed_data/24.2.2_fig4_replot/GO_ck.Rds")
write.csv(ckResSelect,"processed_data/24.2.2_fig4_replot/selectCk.csv")
#== draw heatmap-----------------------------


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

haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(fdrMatrix)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)

hm1 <- Heatmap(exprMatrix[geneSplit$gene,], name = "logChange",
               col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),
               show_row_names = FALSE, show_column_names = FALSE, cluster_rows = T,cluster_columns = F,
               show_row_dend = FALSE, show_column_dend = FALSE,split = geneSplit$cat,
               border =T, cluster_row_slices = F, row_gap = unit(2, "mm"), row_title_rot = 0,
               column_title = "Diff_Expression",raster_quality = 5)


hm_cpm1 <- Heatmap(cpm1[geneSplit$gene,], name = "cpm1",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, 
                   cluster_rows = F,cluster_columns = F,
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   split = geneSplit$cat,border =T, cluster_row_slices = F, 
                   row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "Young",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
#hm_cpm1 <- draw(hm_cpm1)

hm_cpm2 <- Heatmap(cpm2[geneSplit$gene,], name = "cpm2",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, 
                   cluster_rows = F,cluster_columns = F,row_title = NULL,
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   split = geneSplit$cat,border =T, cluster_row_slices = F, 
                   row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "Adult",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
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
pdf("result/24.2.2_fig4_replot//2.2_geneDiff4Panel.pdf",width = 8,height = 10)
draw(hmList)
dev.off()
