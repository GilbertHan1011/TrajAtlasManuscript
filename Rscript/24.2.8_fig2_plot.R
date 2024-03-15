
wt_intergrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
DimPlot(wt_intergrate,group.by = "C7_named",label = TRUE)+	scale_color_npg()
ggsave("result/12.19_figure2_OPC_plot/level2_umap_raster.pdf",width = 10,height = 8)


mypal <- c("#E64B35B2", "#4DBBD5B2", "#3C5488B2","#00A087B2",  "#F39B7FB2","#e377c2ff","#00ff00ff")

names(mypal) <- c("Ob", "Fibroblast","Chondro","MSC", "Lepr+ BMSC", 
                  "Pericyte", "Ly6a+ MSC")

p <- SCpubr::do_DimPlot(wt_intergrate,group.by = "C7_named",raster = T,colors.use = mypal,pt.size = 3)
SCpubr::save_Plot(plot = p,
                  figure_path = "result/24.2.6_fig2_treeplot/",
                  file_name = "C7_raster",output_format = "pdf",
                  create_path = TRUE)
wt_intergrate$C7_named <- full_seurat$C7_named
wt_intergrate$C7_named[wt_intergrate$Organ=="Limb_adult"&wt_intergrate$C7_named=="MSC"] <- "Lepr+ BMSC"
p2 <- SCpubr::do_DimPlot(wt_intergrate,group.by = "C7_named",raster = T,colors.use = mypal,pt.size = 3)
SCpubr::save_Plot(plot = p2,
                  figure_path = "result/24.2.6_fig2_treeplot/",
                  file_name = "C7_mod_raster",output_format = "pdf",
                  create_path = TRUE)
