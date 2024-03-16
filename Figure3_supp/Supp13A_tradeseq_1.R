
library(zellkonverter)
seurat <- readH5AD("../important_processed_data/11.16_dpt.h5ad")

# lineage_meta2 <- read.csv("../important_processed_data/6.22_lineageMeta.csv",row.names = 1)
# sum(lineage_meta$lineage_laFibro==lineage_meta2$lineage_laFibro)
full_seurat <- as.Seurat(seurat,counts = "counts",data = "X")

varGene <- rownames(full_seurat@assays$originalexp@meta.features)[full_seurat@assays$originalexp@meta.features$highly_variable]


cellWeights <- lineage*1
cellWeights <- cellWeights[rowSums(cellWeights)!=0,]
cellWeights <- cellWeights%>%as.data.frame()
set.seed(123)
chondroSample <- sample(colnames(full_seurat)[cellWeights$lineage_chondro==1],size = 8000)
fibroSample <- sample(colnames(full_seurat)[cellWeights$lineage_laFibro==1],size = 8000)
leprSample <- sample(colnames(full_seurat)[cellWeights$lineage_lepr==1],size = 8000)
mesSample <- sample(colnames(full_seurat)[cellWeights$lineage_mesenchyme==1],size = 8000)
cellSub <- c(chondroSample,fibroSample,leprSample,mesSample)
full_seurat <- full_seurat[,colnames(full_seurat)%in%cellSub]
full_seurat <- full_seurat[varGene]
full_seurat <- DietSeurat(full_seurat)
full_seurat <- full_seurat[,colnames(full_seurat)%in%rownames(cellWeights)]
pseduo=full_seurat$pred_dpt
cellWeights=lineage[colnames(full_seurat),]
cellWeights <- cellWeights*1
pseduo=pseduo*cellWeights


BPPRARM <-BiocParallel::bpparam()
BPPRARM$workers<- 10

library(tradeSeq)
icMat <- evaluateK(counts = as.matrix(full_seurat@assays$originalexp@counts),
                   pseudotime =pseduo,
                   cellWeights = cellWeights,
                   k = 3:12,
                   parallel=TRUE,
                   nGenes=100,
                   BPPARAM=BPPRARM,
                   verbose = T)
batch <- factor(full_seurat$batch)
U <- model.matrix(~batch)

mesenchymal_gam <- fitGAM(full_seurat@assays$originalexp@counts,
                          pseudotime =pseduo,
                          cellWeights = cellWeights,
                          nknots = 6,
                          U=U,
                          parallel=TRUE,
                          genes=varGene,
                          BPPARAM=BPPRARM,
                          verbose = T)
saveRDS(mesenchymal_gam,"../important_processed_data/12.26_figgam.Rds")
table(rowData(mesenchymal_gam)$tradeSeq$converged)
rowData(mesenchymal_gam)$assocRes <- associationTest(mesenchymal_gam, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(mesenchymal_gam)$assocRes



lineage1 <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.001)
]
lineage2 <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.001)
]
lineage3 <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_3, "fdr") <= 0.001)
]
lineage4 <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_4, "fdr") <= 0.001)
]
UpSetR::upset(UpSetR::fromList(list(chondro = lineage1, fibro = lineage2,lepr=lineage3,mesenchyme=lineage4)))
geneUnion = unique(c(lineage1,lineage2,lineage3,lineage4))
yhatSmoothUnion <- predictSmooth(mesenchymal_gam, gene = unique(c(lineage1,lineage2,lineage3,lineage4)), nPoints = 50, tidy = FALSE)
heatSmoothUnion <- pheatmap(na.exclude(t(scale(t(yhatSmoothUnion[, 1:200])))),
                            cluster_cols = FALSE,
                            show_rownames = FALSE, 
                            show_colnames = FALSE
)
plotDataAll <-  na.exclude(t(scale(t(yhatSmoothUnion))))
library(circlize)
set.seed(12345)
hmAll<- Heatmap(plotDataAll,cluster_columns = F,show_column_names = F,
                show_row_names = F,row_km = 12,col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
hmAll<- draw(hmAll)

hmOrderAll <- row_order(hmAll)
hmGeneAll <- rownames(hmAll@ht_list$matrix_6@matrix)
hmGeneListAll <- lapply(hmOrderAll, function(x) hmGeneAll[x])
library(stringi)
resAll <- as.data.frame((stri_list2matrix(hmGeneListAll)))
colnames(resAll) <- names(hmGeneListAll)
write.csv(resAll,"result/6.21_traj_hm/12.26_resAll.csv")
geneOrdered <- c(hmOrderAll[["6"]],hmOrderAll[["7"]],hmOrderAll[["2"]],
                 hmOrderAll[["1"]],hmOrderAll[["5"]],hmOrderAll[["3"]],
                 hmOrderAll[["4"]],hmOrderAll[["8"]],hmOrderAll[["9"]],
                 hmOrderAll[["10"]],hmOrderAll[["12"]],hmOrderAll[["11"]])

plotDataAllOrder <- plotDataAll[geneOrdered,]
hmAllOrder <- Heatmap(plotDataAllOrder,cluster_columns = F,show_column_names = F,cluster_rows = F,
                      show_row_names = F,col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
hmAllOrder<- draw(hmAllOrder)
rowSplit <- c(rep(1, length(hmOrderAll[["6"]])),rep(2, length(hmOrderAll[["7"]])),rep(3, length(hmOrderAll[["2"]])),
              rep(4, length(hmOrderAll[["1"]])),rep(5, length(hmOrderAll[["5"]])),rep(6, length(hmOrderAll[["3"]])),
              rep(7, length(hmOrderAll[["4"]])),rep(8, length(hmOrderAll[["8"]])),rep(9, length(hmOrderAll[["9"]])),
              rep(10, length(hmOrderAll[["10"]])),rep(11, length(hmOrderAll[["12"]])),rep(12, length(hmOrderAll[["11"]])))
columnSplit = factor(rep(c("Chondro","Fibro","Lepr","Mesenchyme"),  each=50),levels = c("Chondro","Fibro","Lepr","Mesenchyme"))
symbolOrdered <- hmGeneAll[geneOrdered]
symbolOrderCluster <- cbind(symbolOrdered,rowSplit)
hmAllOrder <- Heatmap(plotDataAllOrder,cluster_columns = F,show_column_names = F,cluster_rows = F,
                      show_row_names = F,col=colorRamp2(c(-2, 0, 4), c("blue", "white", "red")),
                      row_split = rowSplit,column_split = columnSplit)
hmAllOrder<- draw(hmAllOrder)
labels <- c("Ddit3","Bcl2l11","Cited2","Bcl11b","Fgf18","Gbp7","Bst2","Notch1","Jag1","Hey1","Ntrk2","Ngf")
labelsT <- intersect(labels, geneUnion)

mes_meta2 <- full_seurat@meta.data[c("C7_named","C19_named","C36_named","C49_named","anno_level_4")]%>%as.data.frame()
mes_meta2 <- merge(pseduo,mes_meta2,by=0,all=TRUE)%>%column_to_rownames("Row.names")
mes_meta <- full_seurat@meta.data[c("C7_named","C19_named")]%>%as.data.frame()
mes_meta <- merge(pseduo,mes_meta,by=0,all=TRUE)%>%column_to_rownames("Row.names")
mes_meta <- mes_meta2
line1_meta <- na.exclude(mes_meta[c(1,5,6,7)])%>%filter(.[[1]]>0)
line2_meta <- na.exclude(mes_meta[c(2,5,6,7)])%>%filter(.[[1]]>0)
line3_meta <- na.exclude(mes_meta[c(3,5,6,7)])%>%filter(.[[1]]>0)
line4_meta <- na.exclude(mes_meta[c(4,5,6,7)])%>%filter(.[[1]]>0)

#update function----------------
library(Hmisc) 
# meta <- line1_meta
# meta <- meta%>%
#   arrange(Lineage1)%>%
#   mutate(rnkVal = 1:length(meta[,1])) %>%
#   mutate(rnk=as.numeric(cut2(rnkVal, g=100)))
# meta$rnk <- as.numeric(cut2(meta$Lineage4, g=100))

makeHa <- function(meta){
  meta_plot <-  meta%>%
    mutate(rnk=floor(.[[1]]/(max(.[[1]])+0.001) * 50)+1)
  tabMetaIdent<-table(meta_plot[c("C7_named","rnk")])%>%as.matrix()
  tabMetaTime<-table(meta_plot[c("C19_named","rnk")])%>%as.matrix()
  tabMetaC36<-table(meta_plot[c("C36_named","rnk")])%>%as.matrix()
  maxVal <- apply(tabMetaIdent,2,function(x) names(which.max(x)))
  maxValTime <- apply(tabMetaTime,2,function(x) names(which.max(x)))
  maxValC36 <- apply(tabMetaC36,2,function(x) names(which.max(x)))
  return(list(maxVal,maxValTime,maxValC36))
}
line1_max <- makeHa(line1_meta)
line2_max <- makeHa(line2_meta)
line3_max <- makeHa(line3_meta)
line4_max <- makeHa(line4_meta)
my_color_c7<-brewer.pal(5,'Spectral')
my_color_c19<-brewer.pal(9,'Set1')
my_color_c36<-colorRampPalette(brewer.pal(8,'Set2'))(13)

line2_max[[1]] <- c(rep("Fibroblast",17),line2_max[[1]][17:49])
line2_max[[2]] <- c(rep("Fibro",17),line2_max[[2]][17:49])
line2_max[[3]] <- c(rep("Aspn.Fibro",17),line2_max[[3]][17:49])

test=line1_meta%>%mutate(rnk=floor(.[[1]]/(max(.[[1]])+0.001) * 50)+1)

identHa <- c(line1_max[[1]],line2_max[[1]],line3_max[[1]],line4_max[[1]])
timeHa <- c(line1_max[[2]],line2_max[[2]],line3_max[[2]],line4_max[[2]])
C36Ha <- c(line1_max[[3]],line2_max[[3]],line3_max[[3]],line4_max[[3]])
levels(identHa) <- c("Chondro", "MSC", "Fibroblast",  "Lepr+ BMSC","Ob")
levels(timeHa) <- c("Mature Chondro","Fibro", "Lepr+ BMSC", 
                    "Irx1+ MSC", "Early.MSC",
                    "Late.MSC", "HC", "Pre-ob", "Ob")
names(my_color_c7) <- levels(identHa)
names(my_color_c19) <- levels(timeHa) 
names(my_color_c36) <- unique(C36Ha)

haAll = HeatmapAnnotation(
  level_1 = factor(identHa,levels = levels(identHa)), level_2=factor(timeHa,levels = levels(timeHa)),level_3=C36Ha,
  col = list(
    level_1=my_color_c7,
    level_2=my_color_c19,
    level_3=my_color_c36
  )
)
hmAllOrder <- Heatmap(plotDataAllOrder,cluster_columns = F,show_column_names = F,cluster_rows = F,
                      show_row_names = F,col=colorRamp2(c(-2, 0, 3), c("DeepSkyBlue3", "white", "red")),
                      row_split = rowSplit,column_split = columnSplit,top_annotation=haAll)
hmAllOrder<- draw(hmAllOrder)
  labelsT <- intersect(labels, geneUnion)
  position <- c()
  for (i in labelsT){
    position <- c(position, which(is.element(rownames(plotDataAllOrder), i) == TRUE))
  }
  har <- rowAnnotation(link = anno_mark(at = position, labels = labelsT, 
                                        labels_gp = gpar(fontsize = 9, fontface = "bold"), link_width = unit(1.5, "cm")))
hmAllOrder <- Heatmap(plotDataAllOrder,cluster_columns = F,show_column_names = F,cluster_rows = F,
                      show_row_names = F,col=colorRamp2(c(-2, 0, 3), c("DeepSkyBlue3", "white", "red")),
                      row_split = rowSplit,column_split = columnSplit,top_annotation=haAll,right_annotation = har)
pdf("result/6.21_traj_hm/12.26_traj_hm.pdf",width = 10,height = 12)
hmAllOrder<- draw(hmAllOrder)
dev.off()

hmAllOrder <- ComplexHeatmap::Heatmap(plotDataAllOrder,cluster_columns = F,show_column_names = F,cluster_rows = F,
                      show_row_names = F,col=colorRamp2(c(-2, 0, 3), c("DeepSkyBlue3", "white", "red")),
                      row_split = rowSplit,column_split = columnSplit,top_annotation=haAll,right_annotation = har,use_raster = T,border = T)
pdf("result/6.21_traj_hm/12.26_traj_hm_raster.pdf",width = 10,height = 12)
hmAllOrder<- draw(hmAllOrder)
dev.off()
pdf("result/6.21_traj_hm/1.18_traj_hm_raster2.pdf",width = 10,height = 12)
hmAllOrder<- draw(hmAllOrder)
dev.off()


symbolOrdered <- hmGeneAll[geneOrdered]
symbolOrderCluster <- cbind(symbolOrdered,paste0("GC_",rowSplit))
write.csv(symbolOrderCluster,"result/6.21_traj_hm/12.26_symbolOrderCluster.csv")


