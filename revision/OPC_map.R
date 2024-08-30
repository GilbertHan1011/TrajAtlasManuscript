library(dplyr)
library(Seurat)
wt_integrate <- readRDS("important_processed_data/5.4_wtintegrate_full_seurat.Rds")

#markerGene <- read.table("3.9_wt_integrate/result/4.26_annotation/wtintegrate_harmonized___markers_all_pruned.tsv",header = T)
markerGene2 <- data.table::fread("3.9_wt_integrate/result/4.26_annotation/wtintegrate_harmonized__pruned_annotation_markers_filtered.txt",header = T)
#markerGene3 <- data.table::fread("3.9_wt_integrate/result/4.26_annotation/comparisons_all_update.csv",header = T)
genes <- markerGene %>% filter(p_val_adj<0.05) %>% select(gene) %>% unique %>% unlist

genes <- markerGene3 %>% filter(p_val_adj<0.05& pct.1>0.4) %>% select(gene)  %>% unique %>% unlist
sscTable <- readxl::read_xlsx("revision/data/Supp_table3_ssc_anno.xlsx")
genesTable <- sscTable$Marker
splitGenes <- strsplit(genesTable, ";") %>% unlist


splitGenes %in% genes %>% sum
geneSpecific <- splitGenes[splitGenes %in% genes]




testGene <- markerGene3[markerGene3$gene%in%splitGenes[[5]]]
testGene[!testGene$parent%in%testGene$cluster_id]

filterMarker <- function(marker,geneFilter){
  marker <- marker %>% filter(p_val_adj<0.05& pct.1>0.2)
  slectGene <- marker[marker$gene%in%geneFilter]
  slectGene <- slectGene[!slectGene$parent%in%slectGene$cluster_id]
  return(slectGene)
}

filterMarker(markerGene2,"Foxa2")
filterMarker(markerGene2,"Gli1")
geneList <- lapply(geneSpecific,filterMarker,marker = markerGene2)


geneCluster <- lapply(geneList,function(x) x$clean_names)

names(geneCluster) <- geneSpecific

FeaturePlot(wt_integrate,"Hes1")
filterMarker(markerGene2,"Hoxa11")


FeaturePlot(wt_integrate,"Sox9")


LimbAdult <- wt_integrate[,wt_integrate$Organ=="Limb_adult"]
LimbEmbryo <- wt_integrate[,wt_integrate$Organ=="Limb_embryo"]


FeaturePlot(LimbEmbryo,"Sox9")
DimPlot(LimbEmbryo,group.by = "anno_level_3")
Idents(LimbEmbryo) <- LimbEmbryo$anno_level_2

sc_sox9_cell <- WhichCells(LimbEmbryo,idents=c("LE_Chondrocyte"))

p1 <- DimPlot(LimbEmbryo,cells.highlight = sc_sox9_cell)
p2 <- FeaturePlot(LimbEmbryo,"Sox9")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2

ggsave("3.9_wt_integrate/result/4.22_stem_cell_anno/sc_Sox9_cell.pdf",width = 10,height = 6)
#filterMarker(markerGene2,"Gli1")

p1 <- FeaturePlot(wt_integrate,"Nes",pt.size = 2)+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p2 <- DimPlot(wt_integrate,cells.highlight = pthlhCell)

wt_integrate$Null <- NA
p3 <- FeaturePlot(wt_integrate,"Gli1",pt.size = 2)+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p4 <- DimPlot(wt_integrate,group.by = "Null")

p3|p4
p3
ggsave("3.9_wt_integrate/result/4.22_stem_cell_anno/Gli1.pdf",width = 5,height = 6)


FeaturePlot(LimbEmbryo,"Hoxa11")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))

DimPlot(LimbEmbryo,group.by = "anno_level_3")
Idents(LimbEmbryo) <- LimbEmbryo$anno_level_3
sc_hoxa11_cell <- WhichCells(LimbEmbryo,idents=c("LE_Early.Mesenchyme"))
p1 <- FeaturePlot(LimbEmbryo,"Hoxa11")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p2 <-  DimPlot(wt_integrate,cells.highlight = sc_hoxa11_cell)
p1|p2

FeaturePlot(LimbAdult,"Nes")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("3.9_wt_integrate/result/4.22_stem_cell_anno/Hoxa11.pdf",width = 5,height = 6)


FeaturePlot(LimbAdult,"Nfatc1")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
