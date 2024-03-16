leaf_level_column = "C137"
leaf_level = 7

stemDF <- data.frame(full_seurat$cellid,full_seurat$C137)
rownames(stemDF) <- full_seurat$cellid
cellPercent <- stemDF%>%group_by(full_seurat.C137)%>%summarise(count = n())

opcData <- as.data.frame(unique(stemDF$full_seurat.C137))
rownames(opcData) <- opcData[,1]
# i="Acta2+ BMSC"
# cells <- intersect(stemDF$full_seurat.cellid,stemCells[[i]])
# stemDF[cells, i] <- i
# cellCount <- stemDF%>%group_by(full_seurat.C137,!!sym(i))%>%count(name = "presence")
# cellCount <- merge(cellCount,cellPercent,by="full_seurat.C137")%>%mutate(percent=presence/count)%>%
#   .[!is.na(cellCount[[i]]),]
# opcData[[i]] <- 0
# opcData[cellCount$full_seurat.C137,][[i]] <- cellCount$percent
for (i in names(stemCells)){
  cells <- intersect(stemDF$full_seurat.cellid,stemCells[[i]])
  stemDF[cells, i] <- i
  cellCount <- stemDF%>%group_by(full_seurat.C137,!!sym(i))%>%dplyr::count(name = "presence")
  cellCount <- merge(cellCount,cellPercent,by="full_seurat.C137")%>%mutate(percent=presence/count)%>%
    .[!is.na(cellCount[[i]]),]
  opcData[[i]] <- 0
  opcData[cellCount$full_seurat.C137,][[i]] <- cellCount$percent
}
colnames(opcData)[[1]] <- "C137"
heatmap_matrix_OPC = as.matrix(opcData[,2:ncol(opcData)])
# cluster---
hm <- ComplexHeatmap::Heatmap(heatmap_matrix_OPC)
colOrder <- ComplexHeatmap::column_order(hm)
heatmap_matrix_OPC <- heatmap_matrix_OPC[,colOrder]

anno_df = annotation_df %>% dplyr::select(cluster_id,clusterlevel,cluster_name = clean_names)
anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})

#== modify edgelist
#edgelist2 <- readRDS("result/4.26_annotation/5.4update_edgelist.Rds")
# edgelist <-edgelist[3:nrow(edgelist),]

edgelist$from[edgelist$from%in%c("C2-1","C2-2")] <- "all"
# 
edgelist <- edgelist%>%dplyr::filter(!(from %in% 
                                         c("C7-7","C19-19","C36-36","C36-35","C49-47","C49-48","C49-49","C90-88","C90-89","C90-90") | to %in%  
                                         c("C7-7","C19-19","C36-36","C36-35","C49-47","C49-48","C49-49","C90-88","C90-89","C90-90")))

# plot cluster tree:

tree_color = "grey70"
circular_tree = plot_cluster_tree(edgelist = edgelist,
                                  leaf_level=8,
                                  anno_df = anno_df ,
                                  metadata=full_seurat@meta.data,
                                  label_size = 4, 
                                  show_genes = TRUE,
                                  vjust_label = -0.5,
                                  edge_color = tree_color, 
                                  node_color = tree_color)
# OPC_short <- data.frame(colnames(heatmap_matrix_OPC), c("GM", "PM", "Ch", "PB", "pS", 
#                                                         "GB", "LB", "HC", "FC", 
#                                                         "CS", "CB", "AM", 
#                                                         "MM", "HM", "AB", "MB", "HP", 
#                                                         "PC", "oS", "EP"))


extract_matches <- function(text){
  words <- str_extract(text,"(?<=\\().*?(?=\\))")
  return(words)
}

matcheslist <- lapply(colnames(heatmap_matrix_OPC),extract_matches)%>%unlist()

colnames(heatmap_matrix_OPC) <-matcheslist
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
  geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name),color="darkred",vjust=-0.25,size=3, fontface="bold")+
  geom_hilight(data=df, aes(node=node, fill=type),  type = "roundrect",alpha=0.1,extend=0.5)+
  scale_fill_manual(values=mypal)
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






ggsave("result/24.2.11_Fig2_supp//normal_tree.pdf",width = 11,height = 15)
Idents(full_seurat) <- full_seurat$C137
selectCluster <-levels(full_seurat)
index <- order(as.numeric(sub("C137-", "", selectCluster)))

#== all heatmap------------------

# Use the index to get the sorted vector
selectCluster <- selectCluster[index]


clusterOrder <- normaltree$data%>%
  dplyr::filter(clusterlevel=="C137")%>%
  arrange(y)%>%
  select(label)%>%
  unlist()

other_genes_remove = rownames(full_seurat@assays$originalexp@counts)[grepl("RP|Gm|Hb[ab]|Rik|-ps|Tomato",rownames(full_seurat@assays$originalexp@counts))]
# make list of all genes that should be removed:
all_exclusion_genes = other_genes_remove
comparisons_all_update <- read.csv("result/4.26_annotation//comparisons_all_update.csv")
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
full_seurat <- full_seurat[,full_seurat$C137%in%clusterOrder]
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
