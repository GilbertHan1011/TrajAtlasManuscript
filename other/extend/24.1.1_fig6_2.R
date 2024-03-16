setwd("/home/gilberthan/Desktop/disk1/limb/8.25_full_integration/")
dpt <- read.csv("../important_processed_data/24.1.4_dpt.csv")
lineage2 <- read.csv("process_data/lineage/1.3_label.csv")
rownames(dpt) <- lineage$X

#dpt <- dpt%>%column_to_rownames("...1")

dpt <- dpt[colnames(extendAtlas),]
extendAtlas$pseudotime <- dpt[["X0"]]
rownames(lineage) <- lineage$X
lineage <- lineage[2:5]
lineage <- apply(lineage,1,as.logical)
lineage <- as.data.frame(lineage)
lineage <- lineage[colnames(extendAtlas),]
lineageRow <- colSums(lineage)>0
extendAtlas$lineage <- lineageRow

extendAtlas@meta.data[c("lineage_laFibro", "lineage_lepr", "lineage_mesenchyme", 
                       "lineage_chondro")] <- lineage[c("lineage_laFibro", "lineage_lepr", "lineage_mesenchyme", 
                                                       "lineage_chondro")]
lineage <- t(lineage)%>%as.data.frame()
get_column_name <- function(row) {
  true_indices <- which(row)
  if (length(true_indices) == 1) {
    return(names(true_indices))
  } else if (length(true_indices) > 1) {
    selected_index <- sample(true_indices, 1)
    return(names(selected_index))
  } else {
    return(NA)  # No TRUE values in the row
  }
}
# result_column <- apply(lineage2, 1, get_column_name)
# result_column <- result_column%>%unlist()
colnames(lineage) <- c("Fibroblast","Lepr+ BMSC","Mesenchyme","Chondrocyte")
result_column2 <- apply(lineage, 1, get_column_name)
extendAtlas$lineage_sum <-result_column2
DimPlot(extendAtlas,group.by = "lineage_sum")
mypal <- pal_npg("nrc", alpha = 0.7)(5)
names(mypal) <- c("none","Chondrocyte","Fibroblast","Lepr+ BMSC","Mesenchyme")
SCpubr::do_DimPlot(extendAtlas,group.by = "lineage_sum",raster = T,pt.size = 5,na.value = "grey90",colors.use = mypal,shuffle  = T,raster.dpi = 4000)
ggsave("result/1.1_fig6/extend_lineage.pdf",width = 6,height = 6)
extendAtlas$pseudotime[is.na(result_column2)] <- NA
SCpubr::do_FeaturePlot(extendAtlas, features = c("pseudotime"),na.value = "grey90",viridis.palette="H",use_viridis = T,raster = T,pt.size = 4,raster.dpi = 4000)
ggsave("result/1.1_fig6/extend_pseudotime.pdf",width = 6,height = 6)

result_column2%>%as.data.frame()%>%write.csv(.,"process_data/traj_diff/1.8_lineage_df.csv")

#== plot trajMap----------------------------
library(ggsci)
Idents(trajMap2) <- trajMap2$Lineage
mypal <- pal_npg("nrc", alpha = 0.7)(4)
names(mypal) <- unique(Idents(trajMap2))
SCpubr::do_DimPlot(sample = trajMap2,pt.size = 5,reduction = "UMAP",colors.use = mypal)









