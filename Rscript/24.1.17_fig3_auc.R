aucSeurat <- zellkonverter::readH5AD("../important_processed_data/12.24_aucAdata.h5ad")
wt_intergrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
selectDf <- read.csv("processed_data/24.1.17_trajAuc/aucSelect.csv",row.names = 1)
selectDf$REACTOME_CELL_CYCLE_CHECKPOINTS
wt_intergrate$REACTOME_CELL_CYCLE_CHECKPOINTS <- NA
wt_intergrate$REACTOME_CELL_CYCLE_CHECKPOINTS[rownames(selectDf)] <- selectDf$REACTOME_CELL_CYCLE_CHECKPOINTS
wt_intergrate$REACTOME_MET_ACTIVATES_PTK2_SIGNALING <- NA
wt_intergrate$REACTOME_MET_ACTIVATES_PTK2_SIGNALING[rownames(selectDf)] <- selectDf$REACTOME_MET_ACTIVATES_PTK2_SIGNALING
wt_intergrate$REACTOME_NCAM1_INTERACTIONS <- NA
wt_intergrate$REACTOME_NCAM1_INTERACTIONS[rownames(selectDf)] <- selectDf$REACTOME_NCAM1_INTERACTIONS
wt_intergrate$REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS <- NA
wt_intergrate$REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS[rownames(selectDf)] <- selectDf$REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS
wt_intergrate$REACTOME_SIGNALING_BY_PDGF <- NA
wt_intergrate$REACTOME_SIGNALING_BY_PDGF[rownames(selectDf)] <- selectDf$REACTOME_SIGNALING_BY_PDGF

SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 4,
                       features = c("REACTOME_NCAM1_INTERACTIONS"),reduction = "X_draw_graph_fa",
                       na.value = "grey90",viridis.palette="H",use_viridis = T,raster.dpi = 2000)
ggsave("result/24.1.17_flexDotPlot/ncam1.pdf",height = 6,width = 6)
SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 4,
                       features = c("REACTOME_CELL_CYCLE_CHECKPOINTS"),reduction = "X_draw_graph_fa",
                       na.value = "grey90",viridis.palette="H",use_viridis = T,raster.dpi = 2000)
ggsave("result/24.1.17_flexDotPlot/REACTOME_CELL_CYCLE_CHECKPOINTS.pdf",height = 6,width = 6)
SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 4,
                       features = c("REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS"),reduction = "X_draw_graph_fa",
                       na.value = "grey90",viridis.palette="H",use_viridis = T,raster.dpi = 2000)
ggsave("result/24.1.17_flexDotPlot/REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS.pdf",height = 6,width = 6)
SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 4,
                       features = c("REACTOME_SIGNALING_BY_PDGF"),reduction = "X_draw_graph_fa",
                       na.value = "grey90",viridis.palette="H",use_viridis = T,raster.dpi = 2000,max.cutoff = 0.25)
ggsave("result/24.1.17_flexDotPlot/REACTOME_SIGNALING_BY_PDGF.pdf",height = 6,width = 6)
SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 4,
                       features = c("REACTOME_SIGNALING_BY_PDGF"),reduction = "X_draw_graph_fa",
                       na.value = "grey90",viridis.palette="H",use_viridis = T,raster.dpi = 2000,max.cutoff = 0.25)
ggsave("result/24.1.17_flexDotPlot/REACTOME_SIGNALING_BY_PDGF.pdf",height = 6,width = 6)
SCpubr::do_FeaturePlot(wt_intergrate, raster = T,pt.size = 4,
                       features = c("REACTOME_MET_ACTIVATES_PTK2_SIGNALING"),reduction = "X_draw_graph_fa",
                       na.value = "grey90",viridis.palette="H",use_viridis = T,raster.dpi = 2000)
ggsave("result/24.1.17_flexDotPlot/REACTOME_MET_ACTIVATES_PTK2_SIGNALING.pdf",height = 6,width = 6)
