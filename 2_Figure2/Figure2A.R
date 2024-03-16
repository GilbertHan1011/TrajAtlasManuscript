library(aplot)
library(rlang)
library(dplyr)
source("../function/fun_plot_figure.R")
source("../function/fun_plotclustertree.R")
results_path_figure2="result/24.2.6_fig2_treeplot/"
annotation_df <- data.table::fread("result/4.26_annotation/wtintegrate_harmonized__pruned_annotation_result.txt")

comparisons_all_update <- read.csv("result/4.26_annotation//comparisons_all_update.csv")
p3 <- readRDS("result/4.26_annotation//5.4_normaltree_full_dotplot.Rds")
full_seurat <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
wtMeta <- readRDS("../important_processed_data/9.10_wt_meta.Rds")
full_seurat$Age <- wtMeta$Age
leaf_level_column = "C49"
leaf_level = 5

load_colors()

#== read edgelist
edgelist <- readRDS("result/4.26_annotation/5.4update_edgelist.Rds")
#== read OPC
stemCells <- readRDS("processed_data/4.22_highlightcell.Rds")

stemCells <- stemCells[c("Acta2+ BMSC",  "Axin2+ MSC", 
                         "Clec11a+ BMSC", "Hypertrophic chondrocyte", "Chondrocyte", 
                         "Early perichondrial cell", "Fox2a+ chondrocyte", "Gli1+ MSC", 
                         "Grem1+ BMSC", "Hhip+ MSC", "Cd168+ skeletal stem/progenitor cell", 
                         "Lepr+ BMSC", "Mcam+ BMSC", "Msx2+ MSC", "Pdgfra+ BMSC", "Prrx1+ MSC", 
                         "Pthlh+ chondrocyte", "pvSSC", "ocSSC", "Hes1+ Perichondrium"
)]
stemDF <- data.frame(full_seurat$cellid,full_seurat$C49)
rownames(stemDF) <- full_seurat$cellid
cellPercent <- stemDF%>%group_by(full_seurat.C49)%>%summarise(count = n())

opcData <- as.data.frame(unique(stemDF$full_seurat.C49))
rownames(opcData) <- opcData[,1]
# i="Acta2+ BMSC"
# cells <- intersect(stemDF$full_seurat.cellid,stemCells[[i]])
# stemDF[cells, i] <- i
# cellCount <- stemDF%>%group_by(full_seurat.C49,!!sym(i))%>%count(name = "presence")
# cellCount <- merge(cellCount,cellPercent,by="full_seurat.C49")%>%mutate(percent=presence/count)%>%
#   .[!is.na(cellCount[[i]]),]
# opcData[[i]] <- 0
# opcData[cellCount$full_seurat.C49,][[i]] <- cellCount$percent
for (i in names(stemCells)){
  cells <- intersect(stemDF$full_seurat.cellid,stemCells[[i]])
  stemDF[cells, i] <- i
  cellCount <- stemDF%>%group_by(full_seurat.C49,!!sym(i))%>%dplyr::count(name = "presence")
  cellCount <- merge(cellCount,cellPercent,by="full_seurat.C49")%>%mutate(percent=presence/count)%>%
    .[!is.na(cellCount[[i]]),]
  opcData[[i]] <- 0
  opcData[cellCount$full_seurat.C49,][[i]] <- cellCount$percent
}
colnames(opcData)[[1]] <- "C49"
heatmap_matrix_OPC = as.matrix(opcData[,2:ncol(opcData)])
# cluster---
hm <- ComplexHeatmap::Heatmap(heatmap_matrix_OPC)
colOrder <- ComplexHeatmap::column_order(hm)
heatmap_matrix_OPC <- heatmap_matrix_OPC[,colOrder]

anno_df = annotation_df %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

#== modify edgelist
edgelist <-edgelist[3:nrow(edgelist),]
edgelist$from[edgelist$from%in%c("C2-1","C2-2")] <- "all"

edgelist <- edgelist%>%dplyr::filter(!(from %in% c("C7-7","C19-19","C36-36","C36-35") | to %in%  c("C7-7","C19-19","C36-36","C36-35")))

# plot cluster tree:

tree_color = "grey70"
circular_tree = plot_cluster_tree(edgelist = edgelist,
                                  leaf_level=6,
                                  anno_df = anno_df ,
                                  metadata=full_seurat@meta.data,
                                  label_size = 4, 
                                  show_genes = TRUE,
                                  vjust_label = -0.5,
                                  edge_color = tree_color, 
                                  node_color = tree_color)
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
heatmap_matrix_OPC <- heatmap_matrix_OPC[,c(
  "Ch","CS","HC","FC","oS","EP","PC",
  "PM","CB", "AM",  "MM", "HM","HP",
  "GB","LB",
  "PB","pS",
  "AB","MB"
)]



heatmap_data3 = full_seurat@meta.data %>% dplyr::select(cellid,Age,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Age) %>% #dplyr::filter(predicted_Campbell!="NA") 
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label")) 
  tidyr::spread(key = Age,value=presence) %>% as.data.frame()

heatmap_matrix_age = as.matrix(heatmap_data3[,2:ncol(heatmap_data3)])
rownames(heatmap_matrix_age) = heatmap_data3[,leaf_level_column]
heatmap_matrix_age[is.na(heatmap_matrix_age)] = 0
heatmap_matrix_age <- heatmap_matrix_age[,c("Organogenesis stage","Fetal stage","Postnatal", "Young Adult",  "Adult","Old"  )]

heatmap_data4 = full_seurat@meta.data %>% dplyr::select(cellid,Organ,!!sym(leaf_level_column)) %>% dplyr::group_by(!!sym(leaf_level_column),Organ) %>% #dplyr::filter(predicted_Campbell!="NA") 
  dplyr::count(name = "presence")  %>% dplyr::group_by(!!sym(leaf_level_column)) %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup() %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label")) 
  tidyr::spread(key = Organ,value=presence) %>% as.data.frame()

heatmap_matrix_organ = as.matrix(heatmap_data4[,2:ncol(heatmap_data4)])
rownames(heatmap_matrix_organ) = heatmap_data4[,leaf_level_column]
heatmap_matrix_organ[is.na(heatmap_matrix_organ)] = 0


normaltree <- ggtree(circular_tree$data)+
  geom_nodelab(aes(x=branch, label=first_cluster_name), color="darkred",vjust=-0.25,size=3,fontface="bold")+
  geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name),color="darkred",vjust=-0.25,size=3, fontface="bold")
normal_tree_heat <- add_heatmap(circular_tree=normaltree,
            heatmap_matrix = heatmap_matrix_OPC,
            heatmap_colors=c(bg_col,"darkred"),
            #scale_limits = c(0,100),
            heatmap_colnames =T,
            legend_title = "OPC Mapping",
            matrix_offset = 0.5,
            matrix_width =0.6,
            colnames_angle=90,
            legend_text_size = 8,
            hjust_colnames=-1,
            na_color = "white",scale_limits = c(0,1),
            heatmap_text_size=5)
normal_tree_heat
#ggsave("result/24.2.6_fig2_treeplot/fig2_test.pdf",width = 8,height = 8)

normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix_age,
                               heatmap_colors=c(bg_col,"darkorange"),
                               #scale_limits = c(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Ag",
                               matrix_offset = 2.9,
                               matrix_width =0.2,
                               colnames_angle=90,
                               legend_text_size = 8,
                               hjust_colnames=-1,
                               na_color = "white",
                               heatmap_text_size=5)
normal_tree_heat
normal_tree_heat = add_heatmap(circular_tree=normal_tree_heat,
                               heatmap_matrix = heatmap_matrix_organ,
                               heatmap_colors=c(bg_col,"dark green"),
                               #scale_limits = orangec(0,100),
                               heatmap_colnames =T,
                               legend_title = "Pct Tissue",
                               matrix_offset = 3.7,
                               matrix_width =0.1,
                               colnames_angle=90,
                               legend_text_size = 8,
                               hjust_colnames=-1,
                               na_color = "white",
                               heatmap_text_size=5)
normal_tree_heat
ggsave("result/24.2.6_fig2_treeplot/normal_tree_text.pdf",width = 11,height = 8)
Idents(full_seurat) <- full_seurat$C49
selectCluster <-levels(full_seurat)
index <- order(as.numeric(sub("C49-", "", selectCluster)))





# Use the index to get the sorted vector
selectCluster <- selectCluster[index]


clusterOrder <- normaltree$data%>%
  dplyr::filter(clusterlevel=="C49")%>%
  arrange(y)%>%
  select(label)%>%
  unlist()

other_genes_remove = rownames(full_seurat@assays$originalexp@counts)[grepl("RP|Gm|Hb[ab]|Rik|-ps|Tomato",rownames(full_seurat@assays$originalexp@counts))]
# make list of all genes that should be removed:
all_exclusion_genes = other_genes_remove

comparisonGene <- comparisons_all_update%>%
  dplyr::filter(cluster_id%in%clusterOrder&!(gene%in%all_exclusion_genes))%>%
  group_by(cluster_id)%>%
  arrange(desc(specificity))%>%
  top_n(1)%>%
  arrange(fct_relevel(cluster_id, clusterOrder))%>%
  ungroup()%>%
  select(gene)%>%
  unlist()
levels(full_seurat) <- clusterOrder
full_seurat <- full_seurat[,full_seurat$C49%in%clusterOrder]
p2 <- DotPlot(full_seurat,
              features = unique(comparisonGene),assay = "originalexp",
              cluster.idents = T,cols = c("lightgrey", "tan1"),dot.scale = 4)+
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90,size = 10,face = "bold",hjust  = 1)
  )
p3 <- p2%>%
  insert_left(normal_tree_heat,width = 1.5)
p3
ggsave("result/24.2.6_fig2_treeplot/normal_tree.pdf",width = 14,height = 8)



#== plot umap-------------------

wt_intergrate <- readRDS("../important_processed_data/5.4_wtintegrate_full_seurat.Rds")
wt_intergrate <- full_seurat
stemCells <- readRDS("processed_data/4.22_highlightcell.Rds")
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

p <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                        group.by = "stemAnno",
                        na.value = "grey90",raster = T,pt.size = 3)
p
SCpubr::save_Plot(plot = p,
                  figure_path = "result/24.2.6_fig2_treeplot/",
                  file_name = "stemCell_highlight_raster",output_format = "pdf",
                  create_path = TRUE)

mycolor<-colorRampPalette(brewer.pal(8,'Spectral'))(18)
names(mycolor) <- levels(wt_intergrate$stemAnno)
p2 <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                   group.by = "stemAnno",
                   na.value = "grey90",raster = T,pt.size = 3,colors.use = mycolor)
SCpubr::save_Plot(plot = p2,
                  figure_path = "result/24.2.6_fig2_treeplot/",
                  file_name = "stemCell_highlight_raster2",output_format = "pdf",
                  create_path = TRUE)
mycolor2<-colorRampPalette(brewer.pal(8,'Set3'))(18)
names(mycolor2) <- levels(wt_intergrate$stemAnno)
p3 <- SCpubr::do_DimPlot(sample = wt_intergrate, 
                         group.by = "stemAnno",
                         na.value = "grey90",raster = T,pt.size = 3,colors.use = mycolor2)
SCpubr::save_Plot(plot = p3,
                  figure_path = "result/24.2.6_fig2_treeplot/",
                  file_name = "stemCell_highlight_raster3",output_format = "pdf",
                  create_path = TRUE)

mycolor2<-colorRampPalette(brewer.pal(8,'Set'))(18)
names(mycolor) <- levels(wt_intergrate$stemAnno)
SCpubr::do_DimPlot(sample = wt_intergrate, 
                   group.by = "stemAnno",
                   na.value = "grey90",raster = T,pt.size = 3,colors.use = mycolor)

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


#== plot overview--------------


df <- data.frame(node=c(47, 48), type=c("A", "B"))
circular_tree+ geom_nodepoint(color="red", alpha=1, size=2)+geom_tippoint(color="#FDAC4F", shape=8, size=1)+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect")

circular_tree2 =ggtree(circular_tree$data,layout = 'circular', branch.length='none',color="grey")+
  geom_nodepoint(color="#FDAC4F", alpha=1, aes(size=nodesize))+geom_tippoint(color="#FDAC4F", aes(size=nodesize))+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect",alpha=0.2)+
  
  circular_tree2

meta <- wt_intergrate@meta.data
treeData <- circular_tree$data
treeDict <- data.frame(treeData$nodesize,treeData$label)
rownames(treeDict) <- treeDict$treeData.label

count1 <- table(meta$C7)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
treeDict[count1$Var1%>%as.character(),][["treeData.nodesize"]] <- count1$nodeSize

count2 <- table(meta$C19)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
count2 <- count2[count2$Var1%in% meta$C19,]
treeDict[count2$Var1%>%as.character(),][["treeData.nodesize"]] <- count2$nodeSize
count3 <- table(meta$C36)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
count3 <- count3[count3$Var1%in% meta$C36,]
treeDict[count3$Var1%>%as.character(),][["treeData.nodesize"]] <- count3$nodeSize

count4 <- table(meta$C49)%>%as.data.frame()%>%mutate(nodeSize = Freq / nrow(meta))
count4 <- count4[count4$Var1%in% meta$C49,]
treeDict[count4$Var1%>%as.character(),][["treeData.nodesize"]] <- count4$nodeSize
treeDict <- is.na()

treeData$nodesize <- as.numeric(treeDict[treeData$label,]$treeData.nodesize)


circular_tree2 =ggtree(treeData,layout = 'circular', branch.length='none',color="grey")+
  geom_nodepoint(color="#FDAC4F", alpha=1, aes(size=nodesize))+geom_tippoint(color="#FDAC4F", aes(size=nodesize))+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect",alpha=0.2)+
  geom_cladelab(node=47, label="another clade",angle=-90,align=TRUE,
                geom='label', fill='lightblue',hjust='right')
circular_tree2

df <- data.frame(node=c(47, 48,49,50,51,52), type=c("MSC", "Ly6a+ MSC","Lepr+ BMSC","Fibroblast","Chondro","Pericyte"))
circular_tree2 =ggtree(treeData,layout = 'circular', branch.length='none',color="grey")+
  geom_nodepoint(color="#FDAC4F", alpha=1, aes(size=nodesize))+geom_tippoint(color="#FDAC4F", aes(size=nodesize))+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect",alpha=0.1,extend=0.5)+
  geom_cladelab(node=47, label="MSC",angle=-100,align=TRUE,
                geom='label',fill=mypal["MSC"],hjust='right',offset=0.2)+
  geom_cladelab(node=48, label="Ly6a+ MSC",angle=-100,align=TRUE,
                geom='label',fill=mypal["Ly6a+ MSC"],hjust='left',offset=0.2)+
  geom_cladelab(node=49, label="Lepr+ BMSC",angle=-100,align=TRUE,
                geom='label',fill=mypal["Lepr+ BMSC"],hjust='right',offset=0.2)+
  geom_cladelab(node=50, label="Fibroblast",angle=-100,align=TRUE,
                geom='label',fill=mypal["Fibroblast"],hjust='left',offset=0.2)+
  geom_cladelab(node=51, label="Chondro",angle=-100,align=TRUE,
                geom='label',fill=mypal["Chondro"],hjust='bottom',offset=0.2)+
  geom_cladelab(node=52, label="Pericyte",angle=-100,align=TRUE,
                geom='label',fill=mypal["Pericyte"],hjust='left',offset=0.2)+
  scale_size(range=c(0,12))+
  scale_fill_manual(values=mypal)
circular_tree2

ggsave("result/24.2.6_fig2_treeplot/tree.pdf")
circular_tree_bk <- circular_tree
