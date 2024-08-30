library(RColorBrewer)
wtMeta <- readRDS("important_processed_data/9.10_wt_meta.Rds")
limbBud <- wtMeta[wtMeta$Organ=="Limb_embryo" & wtMeta$Age == "Fetal stage",]
limbBudPlotDf <- table(limbBud$Sample,limbBud$C7_named)%>%as.data.frame()
limbBudPlotDf =limbBudPlotDf[limbBudPlotDf$Freq>10,]

colnames(limbBudPlotDf) <- c("Sample","Anno","Freq")
limbBudPlotDf$Sample <- as.character(limbBudPlotDf$Sample)
ggplot(limbBudPlotDf, aes(x=Sample, y=Freq,fill=Anno)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_manual(values=my_palette)+
  theme_classic()+labs(x = NULL)+
  theme(
    plot.title = element_text(size = 10),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black",angle = 90),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = "white"))

wide_df <- limbBudPlotDf %>%# Select relevant columns
  pivot_wider(
    names_from = Anno,   # Column with names for the new columns
    values_from = Freq     # Column with values for the new columns
  )

wide_df[is.na(wide_df)] <- 0
wide_df <- wide_df %>% column_to_rownames("Sample")
wide_df <- wide_df/rowSums(wide_df)

limbbudMeta <- wtMeta[wtMeta$Sample%in%rownames(wide_df),][c("Tissue","Age.In.Detail.","Sample")] %>% unique
limbbudMeta$Tissue <- as.character(limbbudMeta$Tissue)
limbbudMeta$Tissue[limbbudMeta$Tissue!="Forelimb"] <- "Hindlimb"
rownames(limbbudMeta) <- limbbudMeta$Sample

tissueHA <- limbbudMeta[rownames(wide_df),]$Tissue
ageHA <- limbbudMeta[rownames(wide_df),]$Age.In.Detail

my_color1<-brewer.pal(3,'Set1')[1:2]
names(my_color1) <- unique(tissueHA)

my_color2<-brewer.pal(4,'Spectral')
names(my_color2) <- c("E15","E16.5" ,"E18","E18.5" )


ha = HeatmapAnnotation(
  tissue = factor(tissueHA),
  age = factor(ageHA),
  col = list(
    tissue=my_color1,
    age=my_color2
  )
)

library(ComplexHeatmap)



pdf("revision/results/fig2d/limbbud_anno_hm.pdf",width = 6,height = 5)
Heatmap(t(wide_df),cluster_rows = F,cluster_columns = T,top_annotation = ha,name = "cell composition",col=colorRamp2(c(0, 0.4, 0.8), c("DeepSkyBlue3", "white", "red")))
dev.off()



rowAnno <- ComplexHeatmap::annotat


longbone <- wtMeta[wtMeta$Organ=="Limb_adult" & wtMeta$Age == "Postnatal",]
longbonePlotDf <- table(longbone$Sample,longbone$C7_named)%>%as.data.frame()
longbonePlotDf <- longbonePlotDf[longbonePlotDf$Freq>10,]
colnames(longbonePlotDf) <- c("Sample","Anno","Freq")
my_palette <- brewer.pal(n = 7, name = "Set1")
longbonePlotDf$Sample <- as.character(longbonePlotDf$Sample)
ggplot(longbonePlotDf, aes(x=Sample, y=Freq,fill=Anno)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_manual(values=my_palette)+
  theme_classic()+labs(x = NULL)+
  theme(
    plot.title = element_text(size = 10),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black",angle = 90),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = "white"))

# meta---------------
longboneMeta <- wtMeta[wtMeta$Sample%in%unique(longbonePlotDf$Sample),][c("Tissue","Tissue.Specific.","Age.In.Detail.","Cre","Sample")] %>% unique
# limbbudMeta$Tissue <- as.character(limbbudMeta$Tissue)
# limbbudMeta$Tissue[limbbudMeta$Tissue!="Forelimb"] <- "Hindlimb"
rownames(longboneMeta) <- longboneMeta$Sample
longboneMeta$Tissue.Specific. <- as.character(longboneMeta$Tissue.Specific.)
longboneMeta$Tissue.Specific.[longboneMeta$Tissue.Specific.%in%c("Limb Bud")] <- "Diaphysis"

wide_df2 <- longbonePlotDf %>%# Select relevant columns
  pivot_wider(
    names_from = Anno,   # Column with names for the new columns
    values_from = Freq     # Column with values for the new columns
  )

wide_df2[is.na(wide_df2)] <- 0
wide_df2 <- wide_df2 %>% column_to_rownames("Sample")
wide_df2 <- wide_df2/rowSums(wide_df2)

Heatmap(t(wide_df2),cluster_rows = F)

tissueHA2 <- longboneMeta[rownames(wide_df2),]$Tissue.Specific.
ageHA2 <- longboneMeta[rownames(wide_df2),]$Age.In.Detail
CreHA <- longboneMeta[rownames(wide_df2),]$Cre

my_color_long_1<-brewer.pal(5,'Set1')
names(my_color_long_1) <- unique(tissueHA2)

my_color_long_2<-brewer.pal(7,'Spectral')
names(my_color_long_2) <- c("P2", "P7.5", "P13", "P19","P21","P23", "P28" )

ha2 = HeatmapAnnotation(
  tissue = factor(tissueHA2),
  age = factor(ageHA2),
  col = list(
    tissue=my_color_long_1,
    age=my_color_long_2
  )
)

pdf("revision/results/fig2d/longbone_anno_hm.pdf",width = 7,height = 5.5)
Heatmap(t(wide_df2),cluster_rows = F,cluster_columns = T,top_annotation = ha2,name = "cell composition",col=colorRamp2(c(0, 0.4, 0.8), c("DeepSkyBlue3", "white", "red")))
dev.off()






my_color1<-brewer.pal(3,'Set1')[1:2]
names(my_color1) <- unique(tissueHA)

my_color2<-brewer.pal(4,'Spectral')
names(my_color2) <- c("E15","E16.5" ,"E18","E18.5" )

