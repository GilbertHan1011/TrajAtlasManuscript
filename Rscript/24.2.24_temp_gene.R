#==
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
exprMatrix <- data.table::fread("processed_data/12.6_leprDev_vs_Steady/exprMatrix.csv")
exprMatrix <- exprMatrix%>%column_to_rownames("V1")
exprMatrix <- exprMatrix[, -ncol(exprMatrix)]
fdrMatrix <-  read.csv("processed_data/12.6_leprDev_vs_Steady/fdrMatrix.csv",row.names = 1)
cpmMatrix <-  read.csv("processed_data/12.6_leprDev_vs_Steady/cpmMatrix.csv",row.names = 1)
cpm1 <- cpmMatrix[1:99]
cpm2 <- cpmMatrix[100:198]
fdrBinaray <- fdrMatrix>1
exprBinary <- exprMatrix*fdrBinaray

tempGene <- rownames(fdrBinaray)[rowSums(fdrBinaray)<10]
fdrTemp <- fdrBinaray[tempGene,]
result=fdrTemp
for (i in 1:ncol(result)) {
  result[, i] <- result[, i] * i
}
result[result==0]=NA
tempMean <- rowMeans(result,na.rm = T)

tempMean <- sort(tempMean)
tempMean <- tempMean[!is.na(tempMean)]
geneName <- names(tempMean)
tempExpr <- exprBinary[geneName,]
tempBinary <- rowSums(tempExpr)>0
geneSplit <- rep("Adult",length(geneName))
geneSplit[tempBinary]="Young"

GeneCat <- rep("Middle",length(geneName))
GeneCat[tempMean>75] <- "Late"
GeneCat[tempMean<25] <- "Early"
GeneCat
geneSplitCat <- paste0(geneSplit,"_",GeneCat)

geneSplitCat <- factor(geneSplitCat,levels = c("Young_Early",  "Young_Middle",  "Young_Late",
                                               "Adult_Early", "Adult_Middle",  "Adult_Late"))
start <- exprBinary[,1:33]
middle <-  exprBinary[,34:66]
end <-  exprBinary[,67:99]
startSum <- rowSums(start)
middleSum <- rowSums(middle)
endSum <- rowSums(end) 

#==to define different category genes-------------------
threshold <- 10
geneWholeUp <- names(startSum)[startSum>threshold&endSum>threshold]

geneWholeDown <- names(startSum)[startSum< -threshold&endSum< -threshold]

length(geneWholeDown)
length(geneWholeUp)

transitionUp=geneName[geneSplit=="Young"]
transitionDown=geneName[geneSplit=="Adult"]
geneLength <- c(length(geneWholeDown),
                length(geneWholeUp),
                length(transitionDown),
                length(transitionUp)
                )
lengthDf <- data.frame("name"=c(
  "Persistent Adult","Persistent Young",
  "Transition Adult","Transition Young"
                         ),geneLength)
mycolor <-colorRampPalette(brewer.pal(8,'Set1'))(4)
ggplot(lengthDf,aes(x=name,y=geneLength,fill=name))+
  geom_bar(stat = "identity")+geom_text(aes(label = geneLength), vjust = -0.5, size = 3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10,face = "bold"), 
        axis.text.y = element_text(size = 10,face = "bold") )+
  scale_fill_manual(values=mycolor)
ggsave("result/24.2.20_fig4supp_tempgene/2.20_tempCount.pdf",width = 6,height = 6)

haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(fdrMatrix)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)
Heatmap(exprMatrix[geneName,], name = "logChange",
        col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),
        show_row_names = FALSE, show_column_names = FALSE, cluster_rows = T,cluster_columns = F,
        show_row_dend = FALSE, show_column_dend = FALSE,split = geneSplit$cat,
        border =T, cluster_row_slices = F, row_gap = unit(2, "mm"), row_title_rot = 0,
        column_title = "Diff_Expression",raster_quality = 5)



hm_cpm1 <- Heatmap(cpm1[geneName,], name = "cpm1",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, 
                   cluster_rows = F,cluster_columns = F,split = geneSplitCat,
                   show_row_dend = FALSE, show_column_dend = FALSE,border =T, cluster_row_slices = F, 
                   row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "Young",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
#hm_cpm1 <- draw(hm_cpm1)

hm_cpm2 <- Heatmap(cpm2[geneName,], name = "cpm2",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, split = geneSplitCat,
                   cluster_rows = F,cluster_columns = F,row_title = NULL,
                   show_row_dend = FALSE, show_column_dend = FALSE,border =T, cluster_row_slices = F, 
                   row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "Adult",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
#hm_cpm2 <- draw(hm_cpm2)

hm_expr <- Heatmap(exprMatrix[geneName,], name = "logChange",
                   col = colorRamp2(c(-15, 0, 15), hcl_palette="PiYG",reverse = T),split = geneSplitCat,
                   show_row_names = FALSE, show_column_names = FALSE, cluster_rows = F,cluster_columns = F,
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   border =T, cluster_row_slices = F, row_gap = unit(2, "mm"), row_title_rot = 0,column_title = "Diff_Expression",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)
#hm_expr <- draw(hm_expr)

hm_fdr <- Heatmap(fdrMatrix[geneName,], name = "fdr",
                  col = colorRamp2(c(0, 10, 20),hcl_palette="Spectral",reverse = TRUE),split = geneSplitCat,
                  show_row_names = FALSE, show_column_names = FALSE, cluster_rows = F,cluster_columns = F,
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  border =T, cluster_row_slices = F, row_gap = unit(2, "mm"),column_title = "FDR",raster_quality = 5,use_raster = T,top_annotation=haPseudutime)


labels <- c("Clstn2", "Flrt3", "Shisa9","Cyp11a1","Sox9", "Lhx1", "Spg11","Hes1","Stk17b")

labelsT <- intersect(labels, geneName)
position <- c()
for (i in labelsT){
  position <- c(position, which(is.element(geneName, i) == TRUE))
}

har <- rowAnnotation(link = anno_mark(at = position, labels = labelsT, 
                                      labels_gp = gpar(fontsize = 15, fontface = "bold"), link_width = unit(1.5, "cm")))
hm_fdr <- Heatmap(fdrMatrix[geneName,], name = "fdr",
                  col = colorRamp2(c(0, 10, 20),hcl_palette="Spectral",reverse = TRUE),split = geneSplitCat,
                  show_row_names = FALSE, show_column_names = FALSE, cluster_rows = F,cluster_columns = F,
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  border =T, cluster_row_slices = F, row_gap = unit(2, "mm"),column_title = "FDR",
                  raster_quality = 5,use_raster = T,top_annotation=haPseudutime,right_annotation = har)


hmList <- hm_cpm1+hm_cpm2+hm_expr+hm_fdr

hmList

pdf("result/24.2.20_fig4supp_tempgene/2.20_tempGeneHm.pdf",width = 10,height = 8)
#hmList <- hm_cpm1+hm_cpm2+hm_expr+hm_fdr
hmList
dev.off()


tempGeneDf <- data.frame(geneName,geneSplitCat)
write.csv(tempGeneDf,"processed_data/24.2.20_figure4_supp/2.20_tempGene.csv")

geneList <- split(tempGeneDf$geneName,tempGeneDf$geneSplitCat)
GO_ck <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
View(GO_ck@compareClusterResult)
ckRes <- GO_ck@compareClusterResult
ckResSub <- ckRes[c("Description","Cluster","p.adjust")]

ckWide <- pivot_wider(ckResSub,names_from = Cluster,values_from = p.adjust)

# ckWide <- ckWide[,c("Description", "Up_Up","Up_0","0_Up","Down_Down","Down_0","0_Down","Up_Down","Down_Up")]

ckWide <- ckWide%>%column_to_rownames("Description")
ckWideMod <- ckWide%>%
  log10()*-1
ckWideMod[is.na(ckWideMod)] <- 0
ckRowList <- list()
for (i in 1:5){
  ckRow <- ckWideMod[i]-rowMeans(ckWideMod[-1])
  ckRowList[[i]] <- ckRow
}

ckRowDf <- data.frame(do.call(cbind, ckRowList))

leprSelectGO <- c()
for (i in 1:5){
  leprSelectGO <- c(leprSelectGO,ckRowDf %>% arrange(desc(across(i)))%>%rownames%>%.[1:4])
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
ggsave("result/24.2.20_fig4supp_tempgene/2.20_go_temp.pdf",width = 8,height = 8)
# library(ReactomePA)
# 
# m2GMT <- read.gmt("data/m2.all.v2023.2.Mm.symbols (1).gmt")
# term <- unique(m2GMT$term)
# reactomeTerm <- grep("REACTOME",term)
# reactomeTerm <- term[reactomeTerm]
# reactomeGMT <- m2GMT[m2GMT$term%in%reactomeTerm,]
# reactome_ck <- compareCluster(geneCluster = geneList, fun = enricher,TERM2GENE=reactomeGMT,pvalueCutoff=0.1)
# 
# ckRes <- reactome_ck@compareClusterResult
# ckResSub <- ckRes[c("Description","Cluster","p.adjust")]
# 
# ckWide <- pivot_wider(ckResSub,names_from = Cluster,values_from = p.adjust)
# 
# ckWide <- ckWide%>%column_to_rownames("Description")
# ckWideMod <- ckWide%>%
#   log10()*-1
# ckWideMod[is.na(ckWideMod)] <- 0
# ckRowList <- list()
# for (i in 1:5){
#   ckRow <- ckWideMod[i]-rowMeans(ckWideMod[-1])
#   ckRowList[[i]] <- ckRow
# }
# 
# ckRowDf <- data.frame(do.call(cbind, ckRowList))
# 
# leprSelectGO <- c()
# for (i in 1:5){
#   leprSelectGO <- c(leprSelectGO,ckRowDf %>% arrange(desc(across(i)))%>%rownames%>%.[1:4])
# }
# 
# ckResSelect <- ckRes%>%filter(Description%in% leprSelectGO)
# ckResSelect$p.adjust <- -log10(ckResSelect$p.adjust)
# 
# ckResSelect$Description <- factor(ckResSelect$Description,levels = unique(leprSelectGO))
# ckResSelect <- dplyr::arrange(ckResSelect,Description)
# reactome_ck@compareClusterResult <- ckResSelect
# dotplot(reactome_ck,)+  theme(
#   axis.text.x = element_text(angle = 90, hjust = 1, size = 15,face = "bold"),  # Adjust text size here
#   axis.text.y = element_text(size = 15,face = "bold"),  # Adjust y-axis text size
#   axis.title = element_text(size = 14),  # Adjust axis title size
#   legend.text = element_text(size = 12),  # Adjust legend text size
#   legend.title = element_text(size = 14) # Adjust legend title size
# )+scale_fill_gradient( low = "DeepSkyblue3", high = "red")
