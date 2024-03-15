#== This script is written to use SCPuber to better plot

library(ggsci)

wt_intergrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
DimPlot(wt_intergrate,group.by = "C7_named",label = TRUE)+	scale_color_npg()
ggsave("result/12.19_figure2_OPC_plot/level2_umap_raster.pdf",width = 10,height = 8)
DimPlot(wt_intergrate,group.by = "C7_named",label = TRUE,raster = F)+	scale_color_npg()
ggsave("result/12.19_figure2_OPC_plot/level2_umap_notRaster.pdf",width = 10,height = 8)
SCpubr::do_DimPlot(wt_intergrate,group.by = "C7_named",raster = T)+	scale_color_npg()
ggsave("result/12.19_figure2_OPC_plot/level2_umap_Raster_scPubr.pdf",width = 10,height = 8)

stemCells <- readRDS("processed_data/4.22_highlightcell.Rds")

sorted_list <- stemCells[order(sapply(stemCells, length), decreasing = TRUE)]
unique_lengths <- unique(sapply(sorted_list, length))
combined_list <- list()

# Loop through unique lengths and combine arrays with the same length
for (length_val in unique_lengths) {
  # Extract arrays with the current length
  arrays_with_length <- sorted_list[sapply(sorted_list, length) == length_val]
  
  # Combine arrays with a comma-separated string as the name
  combined_name <- paste(names(arrays_with_length), collapse = ', ')
  combined_array <- unname(unlist(arrays_with_length))
  
  # Add the combined array to the list
  combined_list[[combined_name]] <- combined_array
}

wt_intergrate$stemAnno <- NA
for (i in seq_len(length(combined_list))){
  wt_intergrate$stemAnno[combined_list[[i]]] <- names(combined_list)[i]
}
p <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                        group.by = "stemAnno",
                        na.value = "grey90",raster = T,pt.size = 1)
SCpubr::save_Plot(plot = p,
                  figure_path = "result/12.19_figure2_OPC_plot/",
                  file_name = "stemCell_highlight_raster",output_format = "pdf",
                  create_path = TRUE)
p <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                        group.by = "Organ",
                        na.value = "grey90",raster = T,pt.size = 5)+	scale_color_npg()
SCpubr::save_Plot(plot = p,
                  figure_path = "result/12.19_figure2_OPC_plot/",
                  file_name = "organ_highlight_raster",output_format = "pdf",
                  create_path = TRUE)
p <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                        group.by = "Organ",
                        na.value = "grey90",raster = T,pt.size = 5)+	scale_color_npg()
SCpubr::save_Plot(plot = p,
                  figure_path = "result/12.19_figure2_OPC_plot/",
                  file_name = "organ_highlight_raster",output_format = "pdf",
                  create_path = TRUE)

wtMeta <- readRDS("../important_processed_data/9.10_wt_meta.Rds")
wt_intergrate$Age <- wtMeta$Age
wtMeta$Age <- unique(wtMeta$Age)
p <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                        group.by = "Age",
                        na.value = "grey90",raster = T,pt.size = 5)+	scale_color_npg()
mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(6)
wt_intergrate$Age <- factor(wt_intergrate$Age,levels = c("Organogenesis stage", "Fetal stage", "Postnatal", 
                                                         "Young Adult","Adult",  "Old"))
p <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                        group.by = "Age",
                        na.value = "grey90",raster = T,pt.size = 5)+	scale_color_manual(values=mycolor)
p

SCpubr::save_Plot(plot = p,
                  figure_path = "result/12.19_figure2_OPC_plot/",
                  file_name = "Age_umap_raster",output_format = "pdf",
                  create_path = TRUE)
p <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                        group.by = "C7_named",
                        na.value = "grey90",raster = T,pt.size = 5)+	scale_color_npg()
SCpubr::save_Plot(plot = p,
                  figure_path = "result/12.19_figure2_OPC_plot/",
                  file_name = "Age_umap_raster",output_format = "pdf",
                  create_path = TRUE)
SCpubr::save_Plot(plot = p,
                  figure_path = "result/12.19_figure2_OPC_plot/",
                  file_name = "C7_umap_raster",output_format = "pdf",
                  create_path = TRUE)
