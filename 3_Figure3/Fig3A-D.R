rm(wt_intergrate)
library(zenith)
library(ggsci)
wt_intergrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
wtMeta <- read.csv("../important_processed_data/10.26_wt_integrate_meta.csv")
wt_intergrate@meta.data <- wtMeta
rownames(wt_intergrate@meta.data) <- wt_intergrate$X

wt_intergrate$start[wt_intergrate$start=="None"]=NA
mypal <- pal_npg("nrc", alpha = 0.7)(5)
names(mypal) <- unique(wt_intergrate$start)[2:6]
SCpubr::do_DimPlot(wt_intergrate,group.by = "start",raster = T,pt.size = 5,na.value = "grey70",reduction = "X_draw_graph_fa",colors.use = mypal)
ggsave("result/12.24_fig3/start_point.pdf",width = 6,height = 6)

Idents(wt_intergrate) <- wt_intergrate$C90_named

endcell <- WhichCells(wt_intergrate,idents = "Bglap3.Ob")
wt_intergrate$start[endcell]="Mature.Ob"
wt_intergrate$start[wt_intergrate$start=="LimbEmbryo_Mesenchyme"|wt_intergrate$start=="Head_Mesenchyme"]="Mesenchyme"
wt_intergrate$start[wt_intergrate$start=="Limb_Hyperchondrocyte"]="Chondrocyte"
wt_intergrate$start[wt_intergrate$start=="LimbAdult_BMSC"]="Lepr+ BMSC"
wt_intergrate$start[wt_intergrate$start=="LimbAdult_Fibro"]="Fibroblast"


C49Tb <- table(wt_intergrate$start,wt_intergrate$C49_named)%>%as.data.frame()
C36Tb <- table(wt_intergrate$start,wt_intergrate$C36_named)%>%as.data.frame()
C19Tb <- table(wt_intergrate$start,wt_intergrate$C19_named)%>%as.data.frame()
C90Tb <- table(wt_intergrate$start,wt_intergrate$C90_named)%>%as.data.frame()
mypal <- pal_npg("nrc", alpha = 0.7)(5)
names(mypal) <- unique(wt_intergrate$start)[2:6]
SCpubr::do_DimPlot(wt_intergrate,group.by = "start",raster = T,pt.size = 5,na.value = "grey90",reduction = "X_draw_graph_fa",colors.use = mypal)
ggsave("result/12.24_fig3/start_point_withend_1.15.pdf",width = 6,height = 6)

lineage2 <- read.csv("../important_processed_data/10.29_slingshot_lineage_mod.csv",row.names= 1)
lineage <- read.csv("../important_processed_data/11.2_dpt_lineage_infer.csv")
lineage <- lineage%>%column_to_rownames("Unnamed..0")
lineage <- apply(lineage, MARGIN = c(1, 2), as.logical)


lineage[rownames(lineage2)]
lineage2 <- apply(lineage2, MARGIN = c(1, 2), as.logical)

duplicateCell <- rownames(lineage2)[rowSums(lineage2)>1]
fibroLepr <- rownames(lineage2)[lineage2[,2] & lineage2[,3]]
lineage2[,3][sample(fibroLepr,size = length(fibroLepr)/2)] <- FALSE
wt_intergrate$lineage <- NA
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
colnames(lineage) <- c("Chondrocyte","Fibroblast","Lepr+ BMSC","Mesenchyme")
result_column2 <- apply(lineage, 1, get_column_name)

# wt_intergrate$lineage <- result_column
wt_intergrate$lineage <- NA
wt_intergrate$lineage[names(result_column2)] <- result_column2
# mypal2 <- pal_npg("nrc", alpha = 0.7)(5)
# names(mypal2) <-colnames(lineage2)
SCpubr::do_DimPlot(wt_intergrate,group.by = "lineage",raster = T,pt.size = 5,na.value = "grey90",reduction = "X_draw_graph_fa",colors.use = mypal,shuffle  = T)
ggsave("result/12.24_fig3/lineage_1.15.pdf",width = 6,height = 6)

dptTime <- read.csv("../important_processed_data/11.19_lightGBM_pred_and_orig.csv")
lightDpt <- dptTime$light_pred
names(lightDpt) <- dptTime$Unnamed..0
wt_intergrate$pseudotime <- NA
wt_intergrate$pseudotime[names(lightDpt)] <- lightDpt
SCpubr::do_FeaturePlot(wt_intergrate, 
                       features = c("pseudotime"),reduction = "X_draw_graph_fa",na.value = "grey90")+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))

ggsave("result/12.24_fig3/pseudotime.pdf",width = 6,height = 6)

SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 3,
                       features = c("pseudotime"),reduction = "X_draw_graph_fa",na.value = "grey90")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 3,
                       features = c("pseudotime"),reduction = "X_draw_graph_fa",
                       na.value = "grey90",viridis.palette="H",use_viridis = T)
ggsave("result/12.24_fig3/pseudotime_1.15.pdf",width = 6,height = 6)

rm(wt_intergrate)
dpt <- wt_intergrate[,dptTime$Unnamed..0]
mat <- dpt@assays$originalexp$data
C2_pathway <- decoupleR::run_aucell(mat, M2gmt,
                                    .source="geneset",
                                    .target="genesymbol",nproc = 20)


