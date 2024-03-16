library(zellkonverter)
#library(BPCells)
library(Seurat)
library(ggsci)
library(SCpubr)
library(clusterProfilter)
setwd("../8.25_full_integration/")
sceasy::convertFormat("../important_processed_data/9.13_full_latent.h5ad", from="anndata", to="seurat",
                      outFile='../important_processed_data/9.13_full_latent.Rds')
sceasy::convertFormat("../important_processed_data/12.30_extend_reflatent.h5ad", from="anndata", to="seurat",
                      outFile='../important_processed_data/12.30_full_latent.Rds')
sceasy::convertFormat("../important_processed_data/24.1.1_reflatent.h5ad", from="anndata", to="seurat",
                      outFile='../important_processed_data/24.1.1_reflatent.Rds')
extendAtlas <- readRDS('../important_processed_data/24.1.1_reflatent.Rds')
#EAtlas <- SeuratObject::as.Seurat(extendAtlas@assays@data$X)
extendMeta <- data.table::fread("../important_processed_data/9.20_intergrate_meta.csv")

extendAtlas@meta.data <- extendMeta
DimPlot(extendAtlas)
DimPlot(extendAtlas,group.by = "transf_ann_level_1_label")

extendAtlas$Core_or_Extand <- extendMeta$Core_or_Extand

SCpubr::do_DimPlot(extendAtlas,"Core_or_Extand",colors.use = )
DimPlot(extendAtlas,group.by = "Core_or_Extand",shuffle = T)
DimPlot(extendAtlas,group.by = "transf_ann_level_2_label")
FeaturePlot(extendAtlas,"transf_ann_level_2_uncert")

extendAtlas=extendAtlas[,extendAtlas$C7_named!="Low quality"]
DimPlot(extendAtlas,group.by = "C7_named")
Idents(extendAtlas)


extendAtlas$C7_core <- extendAtlas$C7_named
extendAtlas$C7_core[extendAtlas$Core_or_Extand=="Extand"]=NA
mypal <- pal_npg("nrc", alpha = 0.7)(7)
names(mypal) <-c("Chondro", "Fibroblast", "Lepr+ BMSC", "Ly6a+ MSC", "MSC" ,"Ob", "Pericyte")
extendAtlas$C7_core <- as.character(extendAtlas$C7_core)
Idents(extendAtlas) <- extendAtlas$Core_or_Extand
SCpubr::do_DimPlot(extendAtlas,raster = T,pt.size = 3,idents.highlight = "Core",
                   na.value = "grey90",shuffle = T,sizes.highlight = 1)
ggsave("../8.25_full_integration/result/1.1_fig6/core_expand_group.pdf",width = 10,height = 8)

SCpubr::do_DimPlot(extendAtlas,raster = T,group.by = "C7_core",pt.size = 3,na.value = "grey90",colors.use = mypal)
ggsave("../8.25_full_integration/result/1.1_fig6/core_expand.pdf",width = 10,height = 8)

uncertainty <- read.csv("../8.25_full_integration/result/9.18_novel_celltype/new_uncert_trans.csv")
rownames(uncertainty) <- uncertainty$X

uncertainty <- uncertainty[colnames(extendAtlas),]
extendAtlas$uncertain_level2 <- uncertainty$C7_named
FeaturePlot(extendAtlas,"uncertain_level2")
FeaturePlot(extendAtlas,"transf_ann_level_2_uncert")
extendAtlas$uncertain_level2[extendAtlas$Core_or_Extand=="Core"]=NA
SCpubr::do_FeaturePlot(extendAtlas,raster = T,features = "uncertain_level2",na.value = "grey90",pt.size = 5,
                       viridis.palette="H",use_viridis = T,,raster.dpi = 2048)
ggsave("../8.25_full_integration/result/1.1_fig6/uncertain_transfer.pdf",width = 8,height = 8)
SCpubr::do_FeaturePlot(extendAtlas,raster = T,features = "transf_ann_level_2_uncert",na.value = "grey90",
                       pt.size = 3,,viridis.palette="H",use_viridis = F,raster.dpi = 4084)

unique(extendAtlas$Stage)
extendAtlas$Stage_universal <- extendAtlas$Stage
extendAtlas$Stage_universal <- as.character(extendAtlas$Stage_universal)
extendAtlas$Stage_universal[extendAtlas$Stage_universal%in%c("Injury(5-FU)", "Injury(Non-Regeneration)", 
                            "Injury(Radiation)", "Injury(Regeneration)")] <- "Injury"
mypal <- pal_npg("nrc", alpha = 0.7)(7)
names(mypal) <- c( "Injury", "Development","Disease", "Steady", "Treated", "Heterotopic ossification ", 
                  "in vitro")

levels(extendAtlas$Stage_universal) <- unique(extendAtlas$Stage_universal)
SCpubr::do_DimPlot(extendAtlas,raster = T,group.by  = "Stage_universal",
                       pt.size = 3,colors.use = mypal)
ggsave("../8.25_full_integration/result/1.1_fig6/stage_universal2.pdf",width = 8,height = 8)


SCpubr::do_DimPlot(extendAtlas,raster = T,group.by  = "Stage_universal",idents.keep = "Injury",
                       pt.size = 3)

DimPlot(extendAtlas,group.by = "Stage_universal")


ggsave("../8.25_full_integration/result/1.1_fig6/stage_universal2.pdf",width = 8,height = 8)



extendAtlas$update_level2 <- as.character(extendAtlas$update_level2)
extendAtlas$update_level2[extendAtlas$update_level2=="MSC"]="typical MSC"
extendAtlas$update_level2[extendAtlas$update_level2=="Timp+ MSC"]="injury MSC"
Idents(extendAtlas) <- extendAtlas$update_level2
SCpubr::do_DimPlot(extendAtlas,raster = T,idents.keep  = c("injury MSC","typical MSC"),na.value = "grey90",
                   pt.size = 2)
ggsave("../8.25_full_integration/result/1.1_fig6/injury_typical_MSC.pdf",width = 8,height = 8)


TimpCell <- colnames(extendAtlas)[extendAtlas$update_level2=="Timp+ MSC"]




#== GO plot relevant-------------------
timpGene <- read.csv("../8.25_full_integration/result/9.18_novel_celltype/9.22_marker.csv",row.names = 1)
MarkerGene <- timpGene%>%
  head(100)%>%select(names)%>%unlist()


goTop100 <- clusterProfiler::enrichGO(MarkerGene,OrgDb = "org.Mm.eg.db",ont = "BP",keyType="SYMBOL")
dotplot(goTop100)
write.csv(goTop100@result,"result/1.1_fig6/go_res.csv")
goList=c('GO:0030199',
         'GO:0031589',
         'GO:0042060',
         'GO:0034341',
         'GO:0070372',
         'GO:2001233',
         'GO:0002040',
         'GO:1901342',
         'GO:0048144',
         "GO:0060326")
goRes <- goTop100@result
goRes <- goRes[goRes$ID%in%goList,]
goTop100@result <- goRes
dotplot(goTop100)
ggsave("result/1.1_fig6/go_res.pdf",width = 6,height = 5)



levels(extendAtlas$Stage_universal) <- unique(extendAtlas$Stage_universal)
SCpubr::do_DimPlot(extendAtlas,raster = T,group.by  = "C7_named",idents.keep = names(mypal),
                   pt.size = 3,colors.use = mypal,raster.dpi = 2048)
SCpubr::do_DimPlot(extendAtlas,raster = T,group.by  = "C7_named",idents.keep = names(mypal),
                   pt.size = 3,colors.use = mypal,raster.dpi = 20)
ggsave("result/1.1_fig6/C7_universal.pdf",width = 8,height = 8)


dotplot(goRes)+  theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 15,face = "bold"),  # Adjust text size here
  axis.text.y = element_text(size = 15,face = "bold"),  # Adjust y-axis text size
  axis.title = element_text(size = 14),  # Adjust axis title size
  legend.text = element_text(size = 12),  # Adjust legend text size
  legend.title = element_text(size = 14)  # Adjust legend title size
)

