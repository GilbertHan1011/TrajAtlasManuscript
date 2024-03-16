leprCell <- rownames(lineage)[lineage[,"Lepr+ BMSC"]]
rownames(wtMeta) <- colnames(wt_intergrate)
meta <- readRDS("../important_processed_data/9.10_wt_meta.Rds")
leprMeta <- meta[leprCell,]

leprMeta$Stage[leprMeta$Stage%in%c("Steady")&leprMeta$Age%in%c("Young Adult")] <- "Development"

groupMeta <- unique(leprMeta[c("Sample","Project","Age","Stage")])
#leprLineage <- wt_intergrate[,leprCell] 
leprMeta$Stage

wt_intergrate$lepr_stage <- wt_intergrate$Stage
wt_intergrate$lepr_stage[!colnames(wt_intergrate)%in%leprCell] <- NA
wt_intergrate$lepr_stage[!wt_intergrate$lepr_stage%in%c("Steady","Development")] <- NA
wt_intergrate$lepr_stage[wt_intergrate$lepr_stage%in%c("Steady")&wt_intergrate$Age%in%c("Young Adult")] <- "Development"

mypal_lepr <- pal_npg("nrc", alpha = 0.7)(2)
names(mypal_lepr) <- c("Steady","Development")
SCpubr::do_DimPlot(wt_intergrate,split.by  = "lepr_stage",idents.keep=c("Steady","Development"),
                   na.value = "grey90",reduction = "X_draw_graph_fa",raster = T,raster.dpi = 1000,pt.size = 2,colors.use = mypal_lepr)
SCpubr::do_DimPlot(wt_intergrate,group.by  = "lepr_stage",idents.keep=c("Steady","Development"),
                   na.value = "grey90",reduction = "X_draw_graph_fa",raster = T,raster.dpi = 1000,pt.size = 2,colors.use = mypal_lepr)
p1 <- SCpubr::do_DimPlot(wt_intergrate,group.by  = "lepr_stage",idents.keep=c("Development"),
                         na.value = "grey90",reduction = "X_draw_graph_fa",raster = T,raster.dpi = 2000,pt.size = 4,colors.use = mypal_lepr)
p2 <- SCpubr::do_DimPlot(wt_intergrate,group.by  = "lepr_stage",idents.keep=c("Steady"),
                         na.value = "grey90",reduction = "X_draw_graph_fa",raster = T,raster.dpi = 2000,pt.size = 4,colors.use = mypal_lepr)
p1+p2
ggsave("result/24.1.20_trajdiff_dotplot/24.1.21_umap_vis.pdf",width = 10,height = 5)


leprPlotfilter <- read.csv("processed_data/11.25_lepr_conposition.csv",row.names = 1)

groupMeta=groupMeta[groupMeta$Sample%in%rownames(leprPlotfilter),]
mypal_project <-  colorRampPalette(brewer.pal(8,'Set1'))(9)
names(mypal_project) <- unique(groupMeta$Project)
mypal_age<-colorRampPalette(brewer.pal(8,'Spectral'))(4)
names(mypal_age) <- c( "Postnatal", "Young Adult", "Adult","Old")

groupMeta$Stage <- as.character(groupMeta$Stage)
ha = rowAnnotation(
  Stage=groupMeta$Stage,
  Project=groupMeta$Project,
  Age=groupMeta$Age,
  col = list(
    Stage=mypal_lepr,
    Project=mypal_project,
    Age=mypal_age
  )
)

hm <- Heatmap(leprPlotfilter,cluster_rows = T,show_row_names = F,cluster_columns = F,
              right_annotation = ha,split = groupMeta$Stage,col=colorRamp2(c(-4, 0, 4), c("DeepSkyBlue3", "white", "red")))
pdf("result/24.1.20_trajdiff_dotplot/24.1.21_compositionHM.pdf",width = 8,height = 4)
draw(hm)
dev.off()

# color_Age <- makeAnno(leprSample$Age,"Pastel2")
# leprSample$Stage[leprSample$Age=="Young Adult"]="Development"
# ha = rowAnnotation(
#   Age=leprSample$Age,
#   col = list(
#     Age=color_Age
#   )
# )
# 
# 
# 
# table(leprMeta$Stage)
# 
# library(ComplexHeatmap)
# 
# designPre <- table(leprMeta$Sample,leprMeta$Stage)%>%as.data.frame%>%filter(Freq!=0)
# designPre$Var2 <- as.character(designPre$Var2)
# designPre <- designPre%>%filter(Freq>300)
# leprDFPre <- leprBinDf%>%filter(Sample%in%designPre$Var1)
# leprPlotfilter=compositionFun(leprDFPre)
# pheatmap(leprPlotfilter)
# 
# wtMeta <- readRDS("../important_processed_data/9.10_wt_meta.Rds")
# dptmeta <- dpt@meta.data
# showColumns <- c("Sample", "Project", "Organ", "Tissue", "Tissue.Specific.", 
#                  "Stage", "Age", "Age.In.Detail.","Bone.Forming.Methods","Origin","Cre")
# wtMeta <- wtMeta[,showColumns]
# metaSampleTable <- wtMeta %>%
#   group_by(Sample) %>%
#   summarise(across(everything(), first, .names = "{.col}"))
# 
# metaSampleTable <- metaSampleTable%>%filter(is.na(metaSampleTable$Cre))
# leprSample <- metaSample
# 
# 
# leprSample <- metaSample%>%
#   filter(Sample%in%rownames(leprPlotfilter))%>%
#   arrange(factor(Sample, levels = rownames(leprPlotfilter)))
# leprSample <- leprSample%>%filter(leprSample$Sample%in%metaSampleTable$Sample)
# leprPlotfilter <- leprPlotfilter[rownames(leprPlotfilter)%in%metaSampleTable$Sample,]
# 
# 
# color_Stage <- makeAnno(leprSample$Stage,"Pastel2")
# ha = rowAnnotation(
#   Stage=leprSample$Stage,
#   col = list(
#     Stage=color_Stage
#   )
# )
# hm <- Heatmap(leprPlotfilter,cluster_rows = T,show_row_names = F,
#               cluster_columns = F,right_annotation = ha,split = leprSample$Stage)
# draw(hm)
# 
# color_Age <- makeAnno(leprSample$Age,"Pastel2")
# leprSample$Stage[leprSample$Age=="Young Adult"]="Development"
# ha = rowAnnotation(
#   Age=leprSample$Age,
#   col = list(
#     Age=color_Age
#   )
# )
# pdf("result/11.24_benchmark_composition/11.24_lepr_composition_hm.pdf")
# hm <- Heatmap(leprPlotfilter,cluster_rows = T,show_row_names = F,cluster_columns = F,right_annotation = ha,split = leprSample$Stage)
# draw(hm)
# dev.off()
# write.csv(leprPlotfilter,"processed_data/11.25_lepr_conposition.csv")
# 
