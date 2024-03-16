library(Seurat)
extendMes <- extendAtlas[,extendAtlas$lineage_mesenchyme=="True"]
extendMes$Stage_diff <- extendMes$Stage_universal
extendMes$Stage_diff[extendMes$Stage_universal!="Injury"]="Non-injury"
DimPlot(extendMes,group.by = "Stage_diff")
SCpubr::do_DimPlot(extendMes,group.by = "Stage_diff")
extendAtlas$Stage_diff <- NA
extendAtlas$Stage_diff[colnames(extendMes)] <- extendMes$Stage_diff
SCpubr::do_DimPlot(extendAtlas,group.by = "Stage_diff",raster = T,pt.size = 8,raster.dpi = 4000,na.value = "grey90")
ggsave("result/1.1_fig6/trajdiff_umap.pdf",width = 6,height = 6)
extendAtlas$C7_named <- as.character(extendAtlas$C7_named)
mypal <- pal_npg("nrc", alpha = 0.7)(7)
names(mypal) <-c("Chondro", "Fibroblast", "Lepr+ BMSC", "Ly6a+ MSC", "MSC" ,"Ob", "Pericyte")
extendAtlas
SCpubr::do_DimPlot(extendAtlas,group.by = "C7_named",raster = T,pt.size = 8,raster.dpi = 4000,colors.use = mypal)
ggsave("result/1.1_fig6/C7_umap.pdf",width = 6,height = 6)
  