#== This script is written to use SCPuber to better plot
stemCells <- readRDS("processed_data/4.22_highlightcell.Rds")
stemCells <- stemCells[c("Acta2+ BMSC",  "Axin2+ MSC", 
                         "Clec11a+ BMSC", "Hypertrophic chondrocyte", "Chondrocyte", 
                         "Early perichondrial cell", "Fox2a+ chondrocyte", "Gli1+ MSC", 
                         "Grem1+ BMSC", "Hhip+ MSC", "Cd168+ skeletal stem/progenitor cell", 
                         "Lepr+ BMSC", "Mcam+ BMSC", "Msx2+ MSC", "Pdgfra+ BMSC", "Prrx1+ MSC", 
                         "Pthlh+ chondrocyte", "pvSSC", "ocSSC", "Hes1+ Perichondrium"
)]


newid <- c("AB","AM","CB","HC","Ch","EP","FC","GM","GB","HM","CS","LB","MB","MM","PB","PM","PC","pS","oS","HP")
names(stemCells) <-newid 
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


wt_intergrate$stemAnno2 <- NA
for (i in seq_len(length(combined_list))){
  wt_intergrate$stemAnno2[combined_list[[i]]] <- names(combined_list)[i]
}

wt_intergrate$stemAnno2 <- factor(wt_intergrate$stemAnno2, levels = c("Ch", "CS", 
                                                                    "HC", "FC", "oS", 
                                                                    "EP", "PC", "PM", 
                                                                    "CB", "AM", "MM", "HM", 
                                                                    "HP", "GB", "LB", 
                                                                    "PB", "pS", "AB, MB"))


mycolor2<-colorRampPalette(brewer.pal(8,'Set3'))(18)
names(mycolor2) <- levels(wt_intergrate$stemAnno2)

p3 <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                         group.by = "stemAnno2",
                         na.value = "grey90",raster = T,pt.size = 3,colors.use = mycolor2,label = T,label.size = 8)

SCpubr::save_Plot(plot = p3,
                  figure_path = "revision/results/fig2a/",
                  file_name = "stemCell_highlight_raster3_2",output_format = "pdf",
                  create_path = TRUE)

chondro <-  c("Chondrocyte(Ch)", "Cd168+ skeletal stem/progenitor cell(CS)", 
         "Hypertrophic chondrocyte(HC)", "Fox2a+ chondrocyte(FC)", "ocSSC(oS)", 
         "Early perichondrial cell(EP)", "Pthlh+ chondrocyte(PC)", "Prrx1+ MSC(PM)", 
         "Clec11a+ BMSC(CB)", "Axin2+ MSC(AM)", "Msx2+ MSC(MM)", "Hhip+ MSC(HM)", 
         "Hes1+ Perichondrium(HP)", "Grem1+ BMSC(GB)", "Lepr+ BMSC(LB)", 
         "Pdgfra+ BMSC(PB)", "pvSSC(pS)", "Acta2+ BMSC(AB), Mcam+ BMSC(MB)")

OPC_short <- data.frame(colnames(heatmap_matrix_OPC), c("GM", "PM", "Ch", "PB", "pS", 
                                                        "GB", "LB", "HC", "FC", 
                                                        "CS", "CB", "AM", 
                                                        "MM", "HM", "AB", "MB", "HP", 
                                                        "PC", "oS", "EP"))



colnames(heatmap_matrix_OPC) <- c("GM", "PM", "Ch", "PB", "pS", 
                                  "GB", "LB", "HC", "FC", 
                                  "CS", "CB", "AM", 
                                  "MM", "HM", "AB", "MB", "HP", 
                                  "PC", "oS", "EP")
# 
opcName=OPC_short$colnames.heatmap_matrix_OPC.
opcName <- opcName[2:20]
stemCells <- stemCells[opcName]
OPC_short$name <- paste0(OPC_short[[1]],"(",OPC_short[[2]],")")
names(stemCells) <- OPC_short$name[2:20]

rownames(OPC_short) <- OPC_short$c..GM....PM....Ch....PB....pS....GB....LB....HC....FC....CS...
OPC_short[colnames(heatmap_matrix_OPC),][1]
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


wt_intergrate$stemAnno <- factor(wt_intergrate$stemAnno, levels = c("Chondrocyte(Ch)", "Cd168+ skeletal stem/progenitor cell(CS)", 
                                                                    "Hypertrophic chondrocyte(HC)", "Fox2a+ chondrocyte(FC)", "ocSSC(oS)", 
                                                                    "Early perichondrial cell(EP)", "Pthlh+ chondrocyte(PC)", "Prrx1+ MSC(PM)", 
                                                                    "Clec11a+ BMSC(CB)", "Axin2+ MSC(AM)", "Msx2+ MSC(MM)", "Hhip+ MSC(HM)", 
                                                                    "Hes1+ Perichondrium(HP)", "Grem1+ BMSC(GB)", "Lepr+ BMSC(LB)", 
                                                                    "Pdgfra+ BMSC(PB)", "pvSSC(pS)", "Acta2+ BMSC(AB), Mcam+ BMSC(MB)"))








