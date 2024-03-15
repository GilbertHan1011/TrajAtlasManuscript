library(ggplot2)
library(purrr)
library(cowplot)
library(gridExtra)
library(magick)

hmCluster <- read.csv("result/10.30_lineage_draw/hm_cluster.csv")
hm <- readRDS("result/10.30_lineage_draw/11.6_lineage_hm.Rds")
aucAssay <- readRDS("../important_processed_data/11.6_scenic_aucassay.Rds")
aucDpt <- CreateSeuratObject(aucAssay)
dpt_bin <- read.csv("../important_processed_data/11.4_dpt_bin_df.csv")


leprCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Lepr1.0","Lepr2.0","Lepr3.0","Lepr4.0","Lepr5.0")]%>%unique()
chondroCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("chondro7.0", "chondro6.0", "chondro4.0", "chondro5.0", 
                                                  "chondro3.0", "chondro2.0", "chondro1.0")]%>%unique()
MesCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Mesenchyme2.0", "Mesenchyme1.0","Mesenchyme3.0","Mesenchyme4.0")]%>%unique()
FibroCell <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro2.0", "Fibro3.0", "Fibro1.0" , "Fibro4.0" )]%>%unique()
Fibro_Mes <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro5.0", "Mesenchyme5.0")]%>%unique()
Fibro_Mes_Lepr <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro6.0", "Mesenchyme6.0","Lepr6.0","Fibro7.0", "Mesenchyme7.0","Lepr7.0")]
Osteo <- dpt_bin$X[dpt_bin$dpt_lineage%in%c("Fibro8.0", "Mesenchyme8.0","Lepr8.0","chondro8.0",
                                            "Fibro9.0", "Mesenchyme9.0","Lepr9.0","chondro9.0",
                                            "Fibro10.0", "Mesenchyme10.0","Lepr10.0","chondro10.0",
                                            "Fibro11.0","Lepr11.0","chondro11.0")]%>%unique()

aucDpt$coverage_node <- "None"
aucDpt$coverage_node[chondroCell] <- "Chondro"
aucDpt$coverage_node[leprCell] <- "Lepr"
aucDpt$coverage_node[MesCell] <- "Mesenchyme"
aucDpt$coverage_node[FibroCell] <- "Fibro"
aucDpt$coverage_node[Fibro_Mes] <- "Fibro_Mes"
aucDpt$coverage_node[Fibro_Mes_Lepr] <- "Fibro_Mes_Lepr"
aucDpt$coverage_node[Osteo] <- "Osteo"


dpt_dict <- read.csv("../unimportant_processed_data/11.4_dpt_lineage_temp.csv",row.names = 1)

dpt_bin <- dpt_bin %>%
  left_join(dpt_dict, by = "dpt_lineage")

hmcell <- dpt_bin$hmLineage
names(hmcell) <- dpt_bin$X
aucDpt$hmLineage <- hmcell[colnames(aucDpt)]
Idents(aucDpt) <- aucDpt$hmLineage
avg <- AverageExpression(aucDpt,return.seurat = F)
avg <- avg$RNA
avg <- as.data.frame(avg)
avg_rowname <- avg%>%
  rownames_to_column("Pathway")


position <- read.csv("processed_data/24.2.4_fig3_pic_add/position.csv")

plotFun <- function(pathway){
  pathway <- avg[pathway,]%>%t%>%as.data.frame()
  pathway[is.na(pathway)]
  position$pathway <- pathway[position$Point,] 
  
  g <- ggplot(position,aes(x=X,y=Y,color=pathway))+
    geom_point(size=6.7)+
    theme_bw()+  theme(axis.line = element_line(colour = "black"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(),
                       legend.position = "none"
    ) + lims(x = c(0, 4), y = c(0, 4)) +  scale_color_distiller(palette = "Spectral")
  ggdraw() +
    draw_plot(g)+
    draw_image(stock_image3, x = 0, y = 0, scale=1)
}
stock_image3 <- image_read("processed_data/24.2.4_fig3_pic_add/Fig3_overview_4.png")
plotFun("Sox9(+)")
ggsave("result/24.2.15_supp3/regulon_Sox9.pdf",width = 5,height = 3)
plotFun("Sox8(+)")
ggsave("result/24.2.15_supp3/regulon_Sox8.pdf",width = 5,height = 3)
plotFun("Tead1(+)")
ggsave("result/24.2.15_supp3/regulon_Tead1.pdf",width = 5,height = 3)
plotFun("Tbx1(+)")
ggsave("result/24.2.15_supp3/regulon_Tbx1.pdf",width = 5,height = 3)
plotFun("Sox7(+)")
ggsave("result/24.2.15_supp3/regulon_Sox7.pdf",width = 5,height = 3)
plotFun("Scx(+)")
ggsave("result/24.2.15_supp3/regulon_Scx.pdf",width = 5,height = 3)
plotFun("Hoxa11(+)")
ggsave("result/24.2.15_supp3/regulon_Hoxa11.pdf",width = 5,height = 3)
plotFun("Twist1(+)")
ggsave("result/24.2.15_supp3/regulon_Twist1.pdf",width = 5,height = 3)
plotFun("Xbp1(+)")
ggsave("result/24.2.15_supp3/regulon_Xbp1.pdf",width = 5,height = 3)
plotFun("Twist(+)")
ggsave("result/24.2.15_supp3/regulon_Twist.pdf",width = 5,height = 3)
plotFun("Runx2(+)")
ggsave("result/24.2.15_supp3/regulon_Sp7.pdf",width = 5,height = 3)
plotFun("Scx(+)")
ggsave("result/24.2.15_supp3/regulon_Scx.pdf",width = 5,height = 3)
plotFun("Twist(+)")
ggsave("result/24.2.15_supp3/regulon_Twist.pdf",width = 5,height = 3)
plotFun("Sp7(+)")
ggsave("result/24.2.15_supp3/regulon_Sp7.pdf",width = 5,height = 3)
plotFun("Stat5a(+)")
ggsave("result/24.2.15_supp3/regulon_Stat5a.pdf",width = 5,height = 3)

plotFun("Hoxa11(+)")
ggsave("result/24.2.15_supp3/regulon_Twist.pdf",width = 5,height = 3)
plotFun("Sp7(+)")
ggsave("result/24.2.15_supp3/regulon_Sp7.pdf",width = 5,height = 3)