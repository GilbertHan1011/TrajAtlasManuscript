library(ggplot2)
library(purrr)
library(cowplot)
library(gridExtra)
library(magick)
library(ggplot2)

stock_image <- image_read_svg("processed_data/24.2.4_fig3_pic_add/FIg2_overview.svg")
stock_image2 <- image_read("processed_data/24.2.4_fig3_pic_add/Fig3_overview.png")
any_ggplot1 <- qplot(mtcars$hp, mtcars$mpg)+ theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


test <- ggdraw() +
  draw_image(stock_image, x = 0.25, y = 0.25, scale=10) +
  draw_plot(any_ggplot1)
grid.arrange(grobs=list(test))
ggdraw() +
  draw_image(stock_image, x = 0.25, y = 0.25, scale=0.5,)
ggdraw() +
  draw_plot(any_ggplot1)+
  draw_image(stock_image2, x = 0, y = 0, scale=1)

ggdraw() +
  draw_image(stock_image2, x = 0, y = 0, scale=1)

position$X <- 0.3
position$Y <- 0.87
position$X[3] <- 0.687
position$Y[2] <- 1.54
position <- read.csv("processed_data/24.2.4_fig3_pic_add/position.csv")

g <- ggplot(position,aes(x=X,y=Y))+
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
                     axis.line.y = element_blank()
                     ) + lims(x = c(0, 4), y = c(0, 4)) 

ggdraw() +
  draw_plot(g)+
  draw_image(stock_image2, x = 0, y = 0, scale=1)
ggsave("processed_data/24.2.4_fig3_pic_add/middlefile.pdf",width = 5,height = 3)

dpt=zellkonverter::readH5AD("../important_processed_data/12.24_aucAdata.h5ad")
dptSeurat <- as.Seurat(dpt,counts = "X",data = "X")
dpt_bin=read.csv("../unimportant_processed_data/11.4_dpt_lineage_temp.csv")
dpt_lineage <- read.csv("../important_processed_data/11.4_dpt_bin_df.csv")

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

dptSeurat$coverage_node <- "None"
dptSeurat$coverage_node[chondroCell] <- "Chondro"
dptSeurat$coverage_node[leprCell] <- "Lepr"
dptSeurat$coverage_node[MesCell] <- "Mesenchyme"
dptSeurat$coverage_node[FibroCell] <- "Fibro"
dptSeurat$coverage_node[Fibro_Mes] <- "Fibro_Mes"
dptSeurat$coverage_node[Fibro_Mes_Lepr] <- "Fibro_Mes_Lepr"
dptSeurat$coverage_node[Osteo] <- "Osteo"

dpt_dict <- read.csv("../unimportant_processed_data/11.4_dpt_lineage_temp.csv",row.names = 1)

dpt_bin <- dpt_bin %>%
  left_join(dpt_dict, by = "dpt_lineage")

hmcell <- dpt_bin$hmLineage
names(hmcell) <- dpt_bin$X
dptSeurat$hmLineage <- hmcell[colnames(dptSeurat)]
Idents(dptSeurat) <- dptSeurat$hmLineage
avg <- AverageExpression(dptSeurat,return.seurat = F)
avg <- avg$originalexp
avg <- as.data.frame(avg)
avg_rowname <- avg%>%
  rownames_to_column("Pathway")
pathwayPdgf <- avg["REACTOME-SIGNALING-BY-PDGF",]%>%t%>%as.data.frame()
pathwayPdgf[is.na(pathwayPdgf)]
position$pathwayPdgf <- pathwayPdgf[position$Point,] 

g <- ggplot(position,aes(x=X,y=Y,color=pathwayPdgf))+
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
stock_image3 <- image_read("processed_data/24.2.4_fig3_pic_add/Fig3_overview_4.png")
ggdraw() +
  draw_plot(g)+
  draw_image(stock_image3, x = 0, y = 0, scale=1)
ggsave("result/24.2.4_fig3_pathway/pdgf.pdf",width = 5,height = 3)
ggsave("processed_data/24.2.4_fig3_pic_add/middlefile.pdf",width = 5,height = 3)


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
plotFun("REACTOME-SIGNALING-BY-WNT")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-WNT.pdf",width = 5,height = 3)

plotFun("REACTOME-SIGNALING-BY-BMP")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-BMP.pdf",width = 5,height = 3)

plotFun("REACTOME-SIGNALING-BY-TGFB-FAMILY-MEMBERS")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-TGFB-FAMILY-MEMBERS.pdf",width = 5,height = 3)

plotFun("REACTOME-SIGNALING-BY-TGFB-FAMILY-MEMBERS")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-TGFB-FAMILY-MEMBERS.pdf",width = 5,height = 3)

plotFun("REACTOME-SIGNALING-BY-VEGF")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-VEGF.pdf",width = 5,height = 3)


plotFun("REACTOME-SIGNALING-BY-FGFR")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-FGFR.pdf",width = 5,height = 3)

plotFun("REACTOME-SIGNALING-BY-NOTCH")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-NOTCH.pdf",width = 5,height = 3)

plotFun("REACTOME-EPHRIN-SIGNALING")
ggsave("result/24.2.4_fig3_pathway/REACTOME-EPHRIN-SIGNALING.pdf",width = 5,height = 3)


plotFun("REACTOME-MET-ACTIVATES-PTK2-SIGNALING")
ggsave("result/24.2.4_fig3_pathway/REACTOME-MET-ACTIVATES-PTK2-SIGNALING.pdf",width = 5,height = 3)
plotFun("REACTOME-WNT-MEDIATED-ACTIVATION-OF-DVL")
ggsave("result/24.2.4_fig3_pathway/REACTOME-WNT-MEDIATED-ACTIVATION-OF-DVL.pdf",width = 5,height = 3)
plotFun("REACTOME-SIGNALING-BY-ALK")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-ALK.pdf",width = 5,height = 3)

plotFun("REACTOME-VEGF-LIGAND-RECEPTOR-INTERACTIONS")
ggsave("result/24.2.4_fig3_pathway/REACTOME-VEGF-LIGAND-RECEPTOR-INTERACTIONS.pdf",width = 5,height = 3)
plotFun("REACTOME-SIGNALING-BY-ACTIVIN")
ggsave("result/24.2.4_fig3_pathway/EACTOME-SIGNALING-BY-ACTIVIN.pdf",width = 5,height = 3)

plotFun("REACTOME-INTERLEUKIN-2-SIGNALING")
ggsave("result/24.2.4_fig3_pathway/REACTOME-INTERLEUKIN-2-SIGNALING.pdf",width = 5,height = 3)

plotFun("REACTOME-GROWTH-HORMONE-RECEPTOR-SIGNALING")
ggsave("result/24.2.4_fig3_pathway/REACTOME-GROWTH-HORMONE-RECEPTOR-SIGNALING.pdf",width = 5,height = 3)

plotFun("REACTOME-FGFR2C-LIGAND-BINDING-AND-ACTIVATION")
ggsave("result/24.2.4_fig3_pathway/REACTOME-FGFR2C-LIGAND-BINDING-AND-ACTIVATION.pdf",width = 5,height = 3)

plotFun("REACTOME-ERBB2-ACTIVATES-PTK6-SIGNALING")
ggsave("result/24.2.4_fig3_pathway/REACTOME-ERBB2-ACTIVATES-PTK6-SIGNALING.pdf",width = 5,height = 3)
plotFun("REACTOME-SIGNALING-BY-HEDGEHOG")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-HEDGEHOG.pdf",width = 5,height = 3)
plotFun("REACTOME-SIGNALING-BY-GPCR")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALING-BY-GPCR.pdf",width = 5,height = 3)

plotFun("REACTOME-SIGNALLING-TO-ERKS")
ggsave("result/24.2.4_fig3_pathway/REACTOME-SIGNALLING-TO-ERKS.pdf",width = 5,height = 3)


plotFun("REACTOME-REGULATION-OF-INSULIN-LIKE-GROWTH-FACTOR-IGF-TRANSPORT-AND-UPTAKE-BY-INSULIN-LIKE-GROWTH-FACTOR-BINDING-PROTEINS-IGFBPS")
ggsave("result/24.2.4_fig3_pathway/REACTOME-REGULATION-OF-INSULIN-LIKE-GROWTH-FACTOR-IGF-TRANSPORT-AND-UPTAKE-BY-INSULIN-LIKE-GROWTH-FACTOR-BINDING-PROTEINS-IGFBPS.pdf",width = 5,height = 3)


Idents(dptSeurat) <- dptSeurat$C19_named
Idents(dptSeurat) <- dptSeurat$C7_named
markers <- FindAllMarkers(dptSeurat)
mesMarker <- FindMarkers(dptSeurat,"Early.MSC")
osteoMarker <- FindMarkers(dptSeurat,"Ob")
leprMarker <- FindMarkers(dptSeurat,"Lepr+ BMSC")
fibroMarker <- FindMarkers(dptSeurat,"Fibro")
chondroMarker <- FindMarkers(dptSeurat,"Chondro")
leprMarker <- leprMarker[grep("REACTOME",rownames(leprMarker)),]
fibroMarker <- fibroMarker[grep("REACTOME",rownames(fibroMarker)),]
fibroMarker <- fibroMarker[grep("REACTOME",rownames(fibroMarker)),]
chondroMarker <- chondroMarker[grep("REACTOME",rownames(chondroMarker)),]
mesMarker <- mesMarker[grep("REACTOME",rownames(mesMarker)),]


rownames(leprMarker)


ggplot(position,aes(x=X,y=Y,color=pathwayPdgf))+
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
  ) + lims(x = c(0, 4), y = c(0, 4)) +  scale_color_distiller(palette = "Spectral")
ggsave("result/24.2.4_fig3_pathway/legend.pdf",width = 5,height = 3)

library(zellkonverter)
seurat <- readH5AD("../important_processed_data/11.16_dpt.h5ad")

full_seurat <- as.Seurat(seurat,counts = "counts",data = "X")
mesenchymal_gam <- readRDS("../important_processed_data/12.26_figgam.Rds")
colors <- c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2")

# Create a data frame with a single column containing the colors
data <- data.frame(colors)

# Create a bar plot with the colors
ggplot(data, aes(x = 1, fill = colors)) +
  geom_bar(stat = "count", show.legend = FALSE) +
  scale_fill_identity() +
  labs(title = "Color Palette Visualization", x = NULL, y = NULL) +
  theme_void()
