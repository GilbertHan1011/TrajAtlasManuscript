#setwd("../../../disk1/limb/3.9_wt_integrate/")
library(Seurat)
library(RColorBrewer)
#wt_integrate <- readRDS("../important_processed_data/5.4_wtintegrate_wt_integrate.Rds")

DimPlot(wt_integrate)
seuratmeta <- read.csv("../important_processed_data/5.22_wtintegrate_metadata.csv",row.names = 1)
wt_integrate@meta.data <- seuratmeta
DimPlot(wt_integrate)

wt_integrate$Age.In.Detail. <- factor(wt_integrate$Age.In.Detail.,levels = c( "E9", "E10.5", "E10","E11",  "E11.5", "E12","E12.5", 
                                                                            "E13", "E13.5",  "E14","E14.5","E15",
                                                                            "E15.5", "E16","E16.5", "E17","E17.5","E18", "E18.5","P2", 
                                                                            "P5", "P7.5", "P10", "P13",  "3W","P19" , "P21","P23",  "P28","4W",  
                                                                            "1M", "1.5M",  "2M", "3M",  "4M", "Adult", "16M",  "18M"))

my_color<-colorRampPalette(brewer.pal(8,'Spectral'))(38)

wtMeta <- wt_integrate@meta.data

wtMeta$Age_num <- as.numeric(wtMeta$Age.In.Detail.)
wtGroup <- wtMeta%>%
  group_by(C49_named) %>%
  summarize(mean_value = mean(Age_num))%>%
  arrange(mean_value)

wtMeta$C49_named <- factor(wtMeta$C49_named,levels = wtGroup$C49_named)
sort(wtGroup$mean_value)
wtMeta <- table(wtMeta$C49_named,wtMeta$Age.In.Detail.)

ggplot(wtMeta, aes(x=C49_named, fill=Age.In.Detail.)) + 
  geom_bar(position = "fill")+theme_classic()+
  scale_fill_manual(values=my_color)+
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1,size = 10,face = "bold"),title = NULL)

ggsave("result/24.2.15_supp3/startpoint.pdf",width = 10,height = 5)



cytotraceMeta <- read.csv("../important_processed_data/5.22_wtintegrate_cytotrace_metadata.csv")
wt_integrate$cytotrace <- cytotraceMeta$ct_pseudotime

Idents(wt_integrate) <- wt_integrate$C36_named
wtGroup2 <- wtMeta%>%
  group_by(C36_named) %>%
  summarize(mean_value = mean(cytotrace))%>%
  arrange(mean_value)

wt_integrate$C36_named <- factor(wt_integrate$C36_named,wtGroup2$C36_named)

VlnPlot(wt_integrate,features = "cytotrace",pt.size = 0)+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size = 15,face = "bold"),
                                                               axis.text.y = element_text(angle = 0, vjust = 1, hjust=1,size = 15,face = "bold"),
                                                               axis.title.y = element_text(hjust = 0.5, face="bold", size=15),
                                                               axis.title.x = element_text(hjust = 0.5, face="bold", size=15),
                                                               legend.text = element_text(face="bold", size=12))
ggsave("result/24.2.15_supp3/cytotrace.pdf",width = 15,height = 6)

FeaturePlot(wt_integrate,"cytotrace",reduction = "X_draw_graph_fa",split.by = "Organ")&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
ggsave("result/5.22_start_point/cytotrace_umap.pdf",width = 8,height = 6)


changecell <- colnames(wt_integrate)[wt_integrate$Organ=="Limb_adult"& wt_integrate$C7_named=="MSC"]
changecell2 <- colnames(wt_integrate)[wt_integrate$Organ=="Limb_adult"& wt_integrate$C7_named=="Fibroblast"]
DimPlot(wt_integrate,cells.highlight = changecell,reduction = "X_draw_graph_fa",repel = T)
DimPlot(wt_integrate,cells.highlight = changecell2,reduction = "X_draw_graph_fa",repel = T)
wt_integrate$C7_named[changecell] <- "Lepr+ BMSC"
write.csv(wt_integrate@meta.data,"../important_processed_data/10.26_wt_integrate_meta.csv")


wt_integrate$start[endcell] <- "Bglap3.Ob"
