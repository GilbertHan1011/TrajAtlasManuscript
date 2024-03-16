#== select most conserved osteogenic NMF------------------
meanScore <- rowMeans(valAll)
meanSort=names(sort(meanScore,decreasing = T))

binaryVal <- valAll>0.4
sumVal <- sort(rowSums(binaryVal),decreasing = T)
names(sumVal)

go77 <- readRDS("processed_data/12.12_TRAV/geneset/gsea_go_var/go_77.rds")


i="RAV77"
namesDesend <- names(sort(valAll[i,],decreasing = T))
geneDesendDf <- RAVindex[,i]
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

geneVal <- RAVindex[,i]
geneVal <- geneVal[geneDesend]

geneVal <- scale(geneVal)
colorpalette <- brewer.pal(3, "Blues")
row_ha = rowAnnotation(NMF =geneVal,col=list(NMF=colorRamp2(c(-2, 0,4), c("white","lightblue","blue"))))

labelsT <- c("Sox6","Cd59a")
position <- c()
for (i in labelsT){
  position <- c(position, which(is.element(rownames(coor[geneDesend,namesDesend]), i) == TRUE))
}
har <- rowAnnotation(link = anno_mark(at = position, labels = labelsT, 
                                      labels_gp = gpar(fontsize = 9, fontface = "bold"), link_width = unit(1.5, "cm")))


hm <- Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F,top_annotation = haAll,
              left_annotation = row_ha,use_raster = T,col=colorRamp2(c(-1, 0, 1),hcl_palette="PiYG",reverse = TRUE),raster_quality = 10)
pdf("result/12.26_fig4/12.26_heatmap_Lepr_conserved_RAV77.pdf",width = 8,height = 8)
hm <- draw(hm)
dev.off()
namesDesend[2]
filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Ablation_Matsushita_cxcl2.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))
pdf("result/12.26_fig4/12.26_heatmap_RAV77_Ablation_Matsushita_cxcl2.pdf")

haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(filterMat)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)


hm2 <- Heatmap(filterMat,cluster_columns = F,cluster_rows = F,show_row_names = F,
        show_column_names = F,use_raster = T,right_annotation = har,top_annotation = haPseudutime,
        col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
dev.off()

dotplot(go77)

tf77 <- readRDS("processed_data/12.12_TRAV/geneset/go_C3_tf/gsea_77.rds")
dotplot(tf77)
ggsave("result/12.26_fig4/conserved_RAV77_GO.pdf")




hmList <- hm + hm2
pdf("result/12.26_fig4/12.26_heatmap_RAV77_combine.pdf",width = 8,height = 6)
draw(hmList)
dev.off()



#== mix model to find age-conserved gene parttern-----------------------------
library(lme4)


valAllScale <- t(scale(t(valAll)))
trajMeta <- trajMap@meta.data
trajMeta$Age <-  factor(trajMeta$Age,levels = c("Organogenesis stage", "Fetal stage", "Postnatal", 
                                        "Young Adult","Adult",  "Old"))
trajMeta$Age_num <- as.numeric(trajMeta$Age)
trajMeta <- trajMeta[colnames(valAllScale),]
form <- ~ Age_num + Lineage
fit <- dream(valAllScale, form, trajMeta)
fit <- eBayes(fit)
AgeTb <- topTable(fit, coef = "Age_num", number = Inf)

trajMeta$ageOPC <- paste0(trajMeta$Lineage,"_",trajMeta$Age)

form2 <- ~ Age_num + (1|Lineage)
fit2 <- dream(valAllScale, form2, trajMeta)
fit2 <- eBayes(fit2)
AgeTb2 <- topTable(fit2, coef = "Age_num", number = Inf)

trajMeta$ageOPC <- paste0(trajMeta$Lineage,"_",trajMeta$Age)



go140 <- readRDS("processed_data/12.12_TRAV/geneset/gsea_go_var/go_140.rds")

go217 <- readRDS("processed_data/12.12_TRAV/geneset/gsea_go_var/go_217.rds")
View(go217@result)





#== Age heatmap-----------------
valAllScaleHM <- t(valAllScale)%>%as.data.frame()
valAllScaleHM$ageOPC <- trajMeta$ageOPC
valAllScaleHMGroup <-valAllScaleHM%>%
  group_by(ageOPC) %>%
  summarize(across(everything(), mean))

split_df <- valAllScaleHMGroup %>%
  separate(ageOPC, into = c("Celltype", "Age"), sep = "_(?=[^_]*$)", extra = "merge", fill = "right")

AgeHMdf <- split_df%>%
  dplyr::select("Celltype", "Age","RAV140")
AgeHMdf <- AgeHMdf%>%pivot_wider(names_from=Celltype,values_from = RAV140)
AgeHMdf <- AgeHMdf%>%column_to_rownames("Age")
AgeHMdf <- t(AgeHMdf)
AgeHMdf <- AgeHMdf[c( "Mesenchyme","Fibroblast", "Lepr_BMSC","Chondro"),c("Organogenesis stage", "Fetal stage",  "Postnatal", "Young Adult", "Adult","Old"
                     )]

AgeHMdf <- t(scale(t(AgeHMdf)))
pdf("result/12.26_fig4/12.27_Age_score_RAV140_HM.pdf")
Heatmap(AgeHMdf,cluster_rows = F,cluster_columns = F, col=colorRamp2(c(-2, 0, 2), c("DeepSkyBlue3", "white", "red")))
dev.off()
#?AverageExpression



i="RAV140"
namesDesend <- names(sort(valAll[i,],decreasing = T))
geneDesendDf <- RAVindex[,i]
geneDesendDf <- geneDesendDf[geneDesendDf>0]
geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)]


Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F,split = )
Heatmap(peak[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F, show_row_names = F)

IdentidyBar <- trajMap$Lineage[namesDesend]
AgeBar= trajMap$Age[namesDesend]
AgeBar <- factor(AgeBar,levels = c("Organogenesis stage", "Fetal stage",  "Postnatal", "Young Adult", "Adult","Old"))
#IdentidyBar <- trajMap$Lineage[namesDesend]
my_color_age<-colorRampPalette(brewer.pal(6,'Spectral'))(6)
names(my_color_age) <- levels(AgeBar)
haAll_age = HeatmapAnnotation(
  Age = AgeBar[namesDesend],Lineage = IdentidyBar[namesDesend],
  col = list(
    Age=my_color_age,
    Lineage=my_color_lineage
  )
)

IdentidyBar <- factor(IdentidyBar,levels =  c(  "Lepr_BMSC","Fibroblast","Mesenchyme","Chondro"))
row_ha = rowAnnotation(NMF =geneVal,col=list(NMF=colorRamp2(c(-2, 0,4), c("white","lightblue","blue"))))

hm <- Heatmap(coor[geneDesend,namesDesend],cluster_columns = F,cluster_rows = F,show_column_names = F,
              show_row_names = F,top_annotation = haAll_age,column_split =IdentidyBar[namesDesend],left_annotation = row_ha,
              use_raster = T,col=colorRamp2(c(-1, 0, 1),hcl_palette="PiYG",reverse = TRUE),raster_quality = 10)
pdf("result/12.26_fig4/12.27_age_RAV140_hm_lineage.pdf")
hm
dev.off()

geneVal <- RAVindex[,i]
geneVal <- geneVal[geneDesend]

geneVal <- scale(geneVal)
colorpalette <- brewer.pal(3, "Blues")
row_ha = rowAnnotation(NMF =geneVal,col=list(NMF=colorRamp2(c(-2, 0,4), c("white","lightblue","blue"))))

my_color_lineage<- pal_npg("nrc", alpha = 0.7)(4)
names(my_color_lineage) <- unique(IdentidyBar)
haAll = HeatmapAnnotation(
  Lineage = IdentidyBar,
  col = list(
    Lineage=my_color_lineage
  )
)
labelsT <- c("Dmkn","Tubb3")
for (i in labelsT){
  position <- c(position, which(is.element(rownames(coor[geneDesend,namesDesend]), i) == TRUE))
}
har <- rowAnnotation(link = anno_mark(at = position, labels = labelsT, 
                                      labels_gp = gpar(fontsize = 9, fontface = "bold"), link_width = unit(1.5, "cm")))


geneVal <- RAVindex[,i]
geneVal <- geneVal[geneDesend]

geneVal <- scale(geneVal)
colorpalette <- brewer.pal(3, "Blues")
row_ha = rowAnnotation(NMF =geneVal,col=list(NMF=colorRamp2(c(-2, 0,4), c("white","lightblue","blue"))))



filterMat <- read.csv("processed_data/12.12_TRAV/filter/lepr/Metaphysis_Yang_3.csv",row.names = 1)
filterMat <- filterMat[geneDesend,]
filterMat[is.na(filterMat)] <- 0
filterMat <- t(scale(t(filterMat)))
pdf("result/12.26_fig4/12.26_heatmap_RAV77_Ablation_Matsushita_cxcl2.pdf")

haPseudutime = HeatmapAnnotation(
  Psudotime = seq(0, 1, length.out=ncol(filterMat)),
  col = list(
    Psudotime=colorRamp2(c(0,0.5,1),hcl_palette="Spectral",reverse = T)
  )
)


hm2 <- Heatmap(filterMat,cluster_columns = F,cluster_rows = F,show_row_names = F,
               show_column_names = F,use_raster = T,right_annotation = har,top_annotation = haPseudutime,
               col=colorRamp2(c(-3, 0, 4), c("DeepSkyBlue3", "white", "red")))
hm2


write.csv(geneDesend,"../temp_data/12.7_gene_desend.csv")

pdf("result/12.26_fig4/age_related_TRAV.pdf",width = 10,height = 6)
hmList <- hm + hm2
hmList
dev.off()

go140
# 
# i="RAV217"
# namesDesend <- names(sort(valAll[i,],decreasing = T))
# geneDesendDf <- RAVindex[,i]
# geneDesendDf <- geneDesendDf[geneDesendDf>0]
# geneDesend <- names(geneDesendDf)[order(geneDesendDf,decreasing = T)]
# write.csv(geneDesend,"../temp_data/12.7_gene_desend_217.csv")

dotplot(go140,showCategory=20)


#== compare cluster-------------------

geneList_140 <- RAVindex(RAVmodel)[,140]
geneList_140 <- sort(geneList_140, decreasing = TRUE)
geneList_140 <- names(geneList_140)[geneList_140>0]
geneList_77 <- RAVindex(RAVmodel)[,77]
geneList_77 <- sort(geneList_77, decreasing = TRUE)
geneList_77 <- names(geneList_77)[geneList_77>0]
geneListCompare <- list(geneList_140[1:100],geneList_77[1:100])
names(geneListCompare) <- c("RAV140","RAV77")
goCompare <- clusterProfiler::compareCluster(geneListCompare,enrichGO, keyType = "SYMBOL",ont="BP",OrgDb = "org.Mm.eg.db",
                                pvalueCutoff = 0.05)

dotplot(goCompare,showCategory=80)
ggsave("../temp_data/12.28_go_high.pdf",width = 8,height = 40)
write.csv(go140@result,"../temp_data/12.28_go140.csv")

GOselect <- c("GO:0001503", "GO:0060349", "GO:0060485", "GO:0090497", "GO:0060173", 
              "GO:0090287", "GO:0014032", "GO:0048762", "GO:0048864", "GO:0014033", 
              "GO:0097094", "GO:0016525", "GO:0030326", "GO:0035113", "GO:0050920", 
              "GO:0033688", "GO:0098657", "GO:0060411", "GO:0034620",
              "GO:0019320", "GO:0035115", "GO:0033687", "GO:0048846", "GO:1902284")

goresult <- go140@result%>%filter(ID%in%GOselect)
go140@result <- goresult
dotplot(go140)
ggsave("result/12.26_fig4/12.28_GO_140.pdf",width = 6,height = 6)
#== cross-OPC parttern-----------------------
Idents(trajMap) <- trajMap$Lineage
averageHM <- AverageExpression(trajMap,return.seurat = F)
averageHMTb <- averageHM$TRAV
averageHMTb <- t(scale(t(averageHMTb)))
pheatmap(averageHMTb)
Heatmap(averageHMTb)
Heatmap(as.matrix(trajMap@assays$TRAV@counts),show_column_names = F, show_row_names = F)

