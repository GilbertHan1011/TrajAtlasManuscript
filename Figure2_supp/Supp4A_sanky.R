library(ggalluvial)
library(ggrepel)
source("../function/9.24_color_fun.R")
sankyMeta=full_seurat@meta.data[c("anno_level_2","mix_level_2","C7_named")]


sanky_tab <- table(sankyMeta)%>%as.data.frame()%>%
  filter(Freq>300)
sanky_tab$anno_level_2 <- factor(sanky_tab$anno_level_2,levels=c("LA_Diaphyseal MSC","C_Sox18+ Dermis", "C_Hypodermis",  "LA_Fibroblast",
                                                                 "C_Chondrocyte",  "LA_Chondrocyte","LE_Chondrocyte",
                                                                 "C_Osteoblast","LE_Osteoblast","LA_Osteoblast", "LA_Pericyte","LE_Smooth muscle cell" ))
sanky_tab$mix_level_2  <- factor(sanky_tab$mix_level_2,c( "LA_Diaphyseal MSC","Mix_Chondrocyte","Mix_Mesenchyme", 
                                                          "LA_Fibroblast","C_Hypodermis", 
                                                          "C_Sox18+ Dermis","Mix_Osteoblast", "Mix_Pericyte"))
sanky_tab$C7_named  <- factor(sanky_tab$C7_named,c( "Lepr+ BMSC","Ly6a+ MSC", "Fibroblast","Chondro", "MSC",  "Ob", "Pericyte"))
library(ggsci)
ggplot(as.data.frame(sanky_tab),
       aes(y = Freq, axis1 = anno_level_2, axis2 = mix_level_2,axis3=C7_named)) +
  scale_fill_adaptive(name = "npg", palette = "nrc")+
  geom_alluvium(aes(fill = C7_named, width = 1/12,na.rm = F)) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum))) +
  ggtitle("sanky")+
  theme_classic()
ggsave("../3.9_wt_integrate/result/9.24_integrate_fundation/sanky_C7.pdf",width = 8,height = 6)


sankyMeta=full_seurat@meta.data[c("anno_level_3","merge_id_level3","C19_named")]


sanky_tab <- table(sankyMeta)%>%as.data.frame()%>%
  filter(Freq>300)
sanky_tab$anno_level_3 <- factor(sanky_tab$anno_level_3,levels=c("LA_Lepr+ MSC","LA_Ppbp+ MSC","LA_A2m.Chondrocyte",
                                                                   "C_Middle Mesenchyme", 
                                                                 "C_Sox18+ Dermis","LA_Angptl7.Fibroblast", 
                                                                 "LA_Aspn+ MSC", "C_Chondrocyte", 
                                                                 "LA_Proliferating.Chondrocyte", "LA_Serpina1a.Chondrocyte", 
                                                                 "LE_Chondrocyte progenitor",  "LE_Mature.Chondrocyte", "LE_Hypertrophic.Chondrocyte",  "LA_Chondrocyte progenitor", "LA_Dpt.Fibroblast", "LA_Mature.Chondrocyte",  
                                                                 "LA_Hypertrophic.Chondrocyte",  "C_Early mesenchyme","LE_Early.Mesenchyme",
                                                                 "LE_Middle.Mesenchyme", 
                                                                 "C_Irx1.Mesenchyme","LE_Irx1.Mesenchyme", "C_Hypodermis",  "LE_Late.Mesenchyme", "C_Late Mesenchyme", 
                                                                 "C_Meninges",
                                                                   "LA_Pericyte", "C_Pre-osteoblast", "LA_Metaphysis MSC",
                                                                 "C_Osteoblast","LA_Osteoblast",
                                                                 "LE_Osteoblast", "LE_Smooth muscle cell"))
sanky_tab$merge_id_level3  <- factor(sanky_tab$merge_id_level3,c("LA_Lepr+ MSC", "LA_Ppbp+ MSC",   "Mix_Mature.Chondrocyte",
                                                                 "Mix_Hypertrophic.Chondrocyte", "LA_A2m.Chondrocyte","LA_Serpina1a.Chondrocyte",
                                                                 "LA_Proliferating.Chondrocyte", "LA_Chondrocyte progenitor",
                                                             "Mix_Early.Mesenchyme", "C_Middle Mesenchyme","Mix_Irx1.Mesenchyme", 
                                                             "LE_Middle.Mesenchyme", "LE_Late.Mesenchyme", "C_Hypodermis",  "C_Meninges", "C_Late Mesenchyme",  "C_Sox18+ Dermis", "Mix_Pericyte", 
                                                             "LE_Chondrocyte progenitor", "LA_Aspn+ MSC", 
                                                             "LA_Dpt.Fibroblast","LA_Angptl7.Fibroblast", 
                                                              "Mix_Pre-osteoblast","Mix_Osteoblast"
                                                             ))
sanky_tab$C19_named  <- factor(sanky_tab$C19_named,c("Lepr+ BMSC","Ly6a+ MSC",
                                                     "HC", "Hmmr+ CPC",   "Col1a1.Chondro","Fibro", 
                                                      "CPC", "Cyp26a1.Chondro", "Mature Chondro",   "Early.MSC",  "Irx1+ MSC","Meninges","Middle.MSC", "Late.MSC", "Pre-ob", 
                                                     "Ob", "Tex14+ HC","Pericyte"))

ggplot(as.data.frame(sanky_tab),
       aes(y = Freq, axis1 = anno_level_3, axis2 = merge_id_level3,axis3=C19_named)) +
  scale_fill_adaptive(name = "npg", palette = "nrc")+
  geom_alluvium(aes(fill = C19_named), width = 1/12,na.rm = F,aes.bind = T) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label_repel(stat = "stratum", label.size = 0.01,aes(label = after_stat(stratum))) +
  ggtitle("sanky")+
  theme_classic()
ggsave("../3.9_wt_integrate/result/9.24_integrate_fundation/10.2_C19_sanky.pdf",width = 10,height = 8)
