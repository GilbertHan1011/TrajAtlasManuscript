#== trajdiff--------------------------

# function-----------

makeGeneDf <- function(selectGene,coor, expr, peak,Name){
  selectCoor <- coor[selectGene,]
  selectPeak <- peak[selectGene,]
  selectExpr <- as.data.frame(expr)[selectGene,]
  selectCoor <- selectCoor%>%rownames_to_column("Gene")
  selectPeak <- selectPeak%>%rownames_to_column("Gene")
  selectExpr <- selectExpr%>%rownames_to_column("Gene")
  coorLong <- melt(selectCoor, id.vars = "Gene", variable.name = "Lineage", value.name = "Coor")
  peakLong <- melt(selectPeak, id.vars = "Gene", variable.name = "Lineage", value.name = "Peak")
  exprLong <- melt(selectExpr, id.vars = "Gene", variable.name = "Lineage", value.name = "Expr")
  peakLong <- peakLong %>%
    mutate(Stage = case_when(
      Peak >= 7 ~ 'End',
      Peak >= 3 ~ 'Middle',
      TRUE ~ 'Start'
    ))
  combineDf=coorLong
  combineDf["Expr"]=exprLong['Expr']
  combineDf["Peak"]=peakLong['Stage']
  combineDf=combineDf[c(2,1,3,4,5)]
  combineDf$Peak <- factor(combineDf$Peak,c("Start","Middle","End"))
  return(combineDf)
}

coorAll <- read.csv("processed_data/12.9_trajMap/12.9_coorDf.csv",row.names = 1)
exprAll <- read.csv("processed_data/12.9_trajMap/12.9_sumDf.csv",row.names = 1)
peakAll <- read.csv("processed_data/12.9_trajMap/12.9_maxDf.csv",row.names = 1)
lineageName <- colnames(coorAll)
lineageOPC <-  lapply(strsplit(lineageName, "_sep_"),`[` ,2) %>%unlist()
lineageName <-  lapply(strsplit(lineageName, "_sep_"),`[` ,1) %>%unlist()
selectOPC <- c(lineageOPC=="Lepr_BMSC" & lineageName%in%colnames(exprAttr))
coorSelect <- coorAll[selectOPC]
colnames(coorSelect) <- lineageName[selectOPC]
peakSelect <- peakAll[selectOPC]
colnames(peakSelect) <- lineageName[selectOPC]
exprSelect <- exprAll[selectOPC]
colnames(exprSelect) <- lineageName[selectOPC]






coorSelect <- coorSelect[colnames(exprAttr)]
peakSelect <- peakSelect[colnames(exprAttr)]
exprSelect <- exprSelect[colnames(exprAttr)]

write.csv(coorSelect,"processed_data/24.1.20_trajdiff_dotplot/24.1.20_coorSelect.csv")
write.csv(peakSelect,"processed_data/24.1.20_trajdiff_dotplot/24.1.20_peakSelect.csv")
write.csv(exprSelect,"processed_data/24.1.20_trajdiff_dotplot/24.1.20_exprSelect.csv")



exprSelect <- t(scale(t(exprSelect)))
coorDf_sqrt <- sqrt(coorSelect)
coorDf_sqrt_neg <- sqrt(-coorSelect)
coorDf_sqrt[is.na(coorDf_sqrt)] <- 0
coorDf_sqrt_neg[is.na(coorDf_sqrt_neg)] <- 0
coorDf_sqrt=-coorDf_sqrt_neg+coorDf_sqrt

# 
# exprAttr <- data.table::fread("processed_data/12.6_leprDev_vs_Steady/")%>%as.data.frame()%>%column_to_rownames("V1")
# coorAttr <- data.table::fread("processed_data/12.6_leprDev_vs_Steady/coor_attr_v0.csv")%>%as.data.frame()%>%column_to_rownames("V1")
# peakAttr <- data.table::fread("processed_data/12.6_leprDev_vs_Steady/peak_attr_v0.csv")%>%as.data.frame()%>%column_to_rownames("V1")

columnSplit <- c('group1', 'group1', 'group1', 'group1', 'group1', 'group1',
                 'group1', 'group1', 'group1', 'group1', 'group1', 'group1',
                 'group1', 'group1', 'group1', 'group1', 'group1', 'group1',
                 'group1', 'group1', 'group1', 'group1', 'group2', 'group2',
                 'group2', 'group2', 'group2', 'group2', 'group2', 'group2',
                 'group2', 'group2', 'group2', 'group2', 'group2', 'group2',
                 'group2', 'group2', 'group2', 'group2', 'group2')


exprMatrix <- data.table::fread("processed_data/12.6_leprDev_vs_Steady/exprMatrix.csv")
exprMatrix <- data.table::fread("processed_data/12.6_leprDev_vs_Steady/exprMatrix.csv")
#exprMatrix2 <- data.table::fread("result/12.5_lepr_dev_steady//12.5_exprMatrix.csv")
exprMatrix2 <- data.table::fread("result/12.5_lepr_dev_steady//12.5_exprMatrix.csv")
exprMatrix2 <- exprMatrix2%>%column_to_rownames("V1")
exprMatrix <- exprMatrix%>%column_to_rownames("V1")
kmean=exprMatrix$Kmean
exprMatrix <- exprMatrix[, -ncol(exprMatrix)]
geneHM <- data.frame(rownames(exprMatrix),kmean)

gene1=geneHM$rownames.exprMatrix.[geneHM$kmean==0]


sum(vasGene%in%gene1)
sum(vasGene%in%rownames(cpm1))

View(geneHM)

#== enrichment-------------------------
geneSplit <-  split(geneHM$rownames.exprMatrix., geneHM$kmean)

goEnrichApply <- lapply(geneSplit,FUN = function(x){ enrichGO(x,OrgDb = org.Mm.eg.db,
                                            keyType       = 'SYMBOL',
                                            ont           = "All",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 0.01,
                                            qvalueCutoff  = 0.05)})

goEnrichReactome <- lapply(geneSplit,FUN = function(x){ enricher(x,TERM2GENE=reactomeGMT)})


View(goEnrichApply[[1]]@result)

GO1 <- enrichGO(gene= geneHM$rownames.exprMatrix.[geneHM$kmean==0],
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "All",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)

write.csv(GO1,"processed_data/24.1.20_trajdiff_dotplot/go1.csv")
vasGene <- GO1@result$geneID[1]
vasGene <- strsplit(vasGene, "/")[[1]]
debug(makeGeneDf)
dotDf1 <- makeGeneDf(selectGene = vasGene,coor = coorDf_sqrt,expr = exprSelect,peak = peakSelect)



dot_plot(data.to.plot = dotDf1, 
         dend_y_var = c("Expr","Coor"), 
         size_var = "Expr", shape.scale =8,col_var = "Coor",shape_var="Peak",cols.use= viridis::viridis(20,option="H"),shape_use = c(  "\u25CF", "\u25C6","\u25A0" ))

#== heatmap -----------------
fdrMatrix <-  read.csv("processed_data/12.6_leprDev_vs_Steady/fdrMatrix.csv",row.names = 1)

kmean <- factor(kmean,levels = c(1,3,2,6,5,4,0,7))
hm_cpm1 <- Heatmap(cpm1, name = "cpm1",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, 
                   row_title = NULL,cluster_rows = F,cluster_columns = F, row_split = kmean,
                   show_row_dend = FALSE, show_column_dend = FALSE,column_title = "Group1_CPM",
                   raster_quality = 10,cluster_row_slices = FALSE,border = TRUE)
#hm_cpm1 <- draw(hm_cpm1)

hm_cpm2 <- Heatmap(cpm2, name = "cpm2",
                   col = colorRamp2(c(-3, 0, 3.5),hcl_palette="RdBu",reverse = TRUE),
                   show_row_names = TRUE, show_column_names = FALSE, 
                   row_title = NULL,cluster_rows = F,cluster_columns = F, row_split = kmean,
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   column_title = "Group2_CPM",raster_quality = 10,
                   cluster_row_slices = FALSE,border = TRUE)

hm_expr <- Heatmap(exprMatrix, name = "logChange",
                   col = colorRamp2(c(-15, 0, 15),hcl_palette="PiYG",reverse = TRUE),
                   show_row_names = FALSE, show_column_names = FALSE, 
                   row_title = NULL, cluster_rows = F,cluster_columns = F,
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   row_split = kmean,column_title = "Diff_expression",raster_quality = 10,
                   cluster_row_slices = FALSE,border = TRUE)
#hm_expr <- draw(hm_expr)

hm_fdr <- Heatmap(fdrMatrix, name = "fdr",
                  col = colorRamp2(c(0, 10, 20),hcl_palette="Spectral",reverse = TRUE),
                  show_row_names = FALSE, show_column_names = FALSE, 
                  row_title = NULL, cluster_rows = F,cluster_columns = F,
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  row_split = kmean,column_title = "FDR",raster_quality = 10,
                  cluster_row_slices = FALSE,border = TRUE)


#hm_cpm2 <- draw(hm_cpm2)

ht_list <- hm_cpm1+hm_cpm2+hm_expr+hm_fdr
pdf("result/24.1.20_trajdiff_dotplot//expr_cpm_fdr_hm.pdf",width = 6,height = 10)
hmDraw  <- draw(ht_list)
dev.off()





#== plot1-------------
makeTable <- function(index,goindex){
  Res1 <- goEnrichReactome[[index]]@result%>%
    mutate(logP=-log10(p.adjust))%>%
    select(ID,logP)%>%
    mutate(ID = gsub("REACTOME_", "", ID))
  Res1Sub <- Res1[goindex,]
  Res1Sub$ID <- gsub("REACTOME_","",Res1Sub$ID)
  Res1Sub <- Res1Sub %>%
    mutate(ID = reorder(ID, logP))
  return(Res1Sub)
}


Res1 <- makeTable(2,c(1, 2,4,9 ,13))
ggplot(Res1, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res1_bar.pdf",width = 5,height = 3)



Res2 <- makeTable(4,c(1,4,20,23,33))
ggplot(Res2, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res2_bar.pdf",width = 5,height = 3)



Res3 <- makeTable(3,c(2,4,19,37,61))
ggplot(Res3, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res3_bar.pdf",width = 5,height = 3)

Res4 <- makeTable(7,c(1,5,8,9,11))
ggplot(Res4, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res4_bar.pdf",width = 5,height = 3)



Res5 <- makeTable(6,c(1,3,4,5,8))
ggplot(Res5, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res5_bar.pdf",width = 5,height = 3)



enrichRes <- clusterProfiler::enricher(c(geneSplit[[1]],geneSplit[[5]]),TERM2GENE=reactomeGMT)
makeTable2 <- function(goindex){
  Res1 <- enrichRes@result%>%
    mutate(logP=-log10(p.adjust))%>%
    select(ID,logP)%>%
    mutate(ID = gsub("REACTOME_", "", ID))
  Res1Sub <- Res1[goindex,]
  Res1Sub$ID <- gsub("REACTOME_","",Res1Sub$ID)
  Res1Sub <- Res1Sub %>%
    mutate(ID = reorder(ID, logP))
  return(Res1Sub)
}
Res6 <- makeTable2(c(1,4,10,11,12))
ggplot(Res6, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res6_bar.pdf",width = 5,height = 3)

Res6_2 <- makeTable(5,c(1,2,3,4,5))
ggplot(Res6_2, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res6_2_bar.pdf",width = 5,height = 3)


write.csv(geneHM,"processed_data/24.1.20_trajdiff_dotplot/geneHM.csv")


Res7 <- makeTable(1,c(1,2,3,4,5))
ggplot(Res7, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,nrow(fibroResSub))),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res7_bar.pdf",width = 5,height = 3)


.hmRes <- readRDS("processed_data/12.6_leprDev_vs_Steady/expr_cpm_fdr_go_draw.Rds")






# hm_attr_expr <- Heatmap(exprAttr, name = "expr_attr",
#                         show_row_names = FALSE, show_column_names = FALSE,
#                         row_title = NULL, cluster_rows = F,cluster_columns = T,
#                         show_row_dend = FALSE, show_column_dend = FALSE,column_split = columnSplit,
#                         row_split = kmean,column_title = "expr_attr")
# hm_attr_peak <- Heatmap(peakAttr, name = "peak_attr",
#                         show_row_names = FALSE, show_column_names = FALSE,
#                         row_title = NULL, cluster_rows = F,cluster_columns = T,
#                         show_row_dend = FALSE, show_column_dend = FALSE,column_split = columnSplit,
#                         row_split = kmean,column_title = "peak_attr")+ha
# hm_attr_coor <- Heatmap(coorAttr, name = "coor_attr",
#                         show_row_names = FALSE, show_column_names = FALSE,
#                         row_title = NULL, cluster_rows = F,cluster_columns = T,
#                         show_row_dend = FALSE, show_column_dend = FALSE,column_split = columnSplit,
#                         row_split = kmean,column_title = "coor_attr",)