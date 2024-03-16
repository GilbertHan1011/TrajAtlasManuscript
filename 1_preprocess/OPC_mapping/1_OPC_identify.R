#==1---------------
FeaturePlot(full_seurat,"Pdgfra1")
limbAdult <- limgAdult
limbEmbryo <- full_seurat[,basic_boolean_embro]
rm(limgAdult)
limbAdult <- full_seurat[,basic_boolean_limbadult]
Idents(limbAdult) <- limbAdult$anno_level_2
DimPlot(limbAdult,group.by = "anno_level_2")
sc_pdgfra_cell <- WhichCells(limbAdult,idents=c("LA_Fibroblast", "LA_Diaphyseal MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_pdgfra_cell)
p2 <- FeaturePlot(limbAdult,"Pdgfra")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_pdgfra_cell.pdf",width = 10,height = 6)

#==2---------------
FeaturePlot(limbAdult,"Acta2")
sc_acta2_cell <- WhichCells(limbAdult,idents=c("LA_Pericyte"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_acta2_cell)
p2 <- FeaturePlot(limbAdult,"Acta2")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_acta2_cell.pdf",width = 10,height = 6)

#==3----------------
FeaturePlot(limbAdult,"Lepr")
sc_lepr_cell <- WhichCells(limbAdult,idents=c("LA_Diaphyseal MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_lepr_cell)
p2 <- FeaturePlot(limbAdult,"Lepr")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_Lepr_cell.pdf",width = 10,height = 6)


#==4-----------------------
FeaturePlot(limbAdult,"Cxcl12")
sc_cxcl12_cell <- WhichCells(limbAdult,idents=c("LA_Diaphyseal MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_lepr_cell)
p2 <- FeaturePlot(limbAdult,"Cxcl12")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_cxcl12_cell.pdf",width = 10,height = 6)


#==5-------------------------------------------
Idents(limbAdult) <- limbAdult$anno_level_3
DimPlot(limbAdult,label = T)
FeaturePlot(limbAdult,"Clec11a")
sc_clec11a_cell <- WhichCells(limbAdult,idents=c("LA_Metaphysis MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_clec11a_cell)
p2 <- FeaturePlot(limbAdult,"Clec11a")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_clec11a_cell.pdf",width = 10,height = 6)

#==6---------------------------------

Idents(limbAdult) <- limbAdult$anno_level_2
DimPlot(limbAdult,label = T)
FeaturePlot(limbAdult,"Clec11a")
sc_adipoq_cell <- WhichCells(limbAdult,idents=c("LA_Diaphyseal MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_adipoq_cell)
p2 <- FeaturePlot(limbAdult,"Adipoq")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_Adipoq_cell.pdf",width = 10,height = 6)


#==7---------------------------------

Idents(limbAdult) <- limbAdult$anno_level_3
DimPlot(limbAdult,label = T)
FeaturePlot(limbAdult,"Nes")
sc_adipoq_cell <- WhichCells(limbAdult,idents=c("LA_Diaphyseal MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_adipoq_cell)
p2 <- FeaturePlot(limbAdult,"Nes")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p2
ggsave("result/4.22_stem_cell_anno/sc_Nes_cell.pdf",width = 8,height = 6)

#== 8--------------------------------
Idents(limbAdult) <- limbAdult$anno_level_3
DimPlot(limbAdult,label = T)
FeaturePlot(limbAdult,"Grem1")
sc_grem1_cell <- WhichCells(limbAdult,idents=c("LA_A2m.Chondrocyte","LA_Lepr+ MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_grem1_cell)
p2 <- FeaturePlot(limbAdult,"Grem1")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_grem1_cell.pdf",width = 10,height = 6)

#== 9--------
Idents(limbAdult) <- limbAdult$anno_level_2
FeaturePlot(limbAdult,"Ly6a")
sc_pvssc_cell <- WhichCells(limbAdult,idents=c("LA_Fibroblast"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_pvssc_cell)
p2 <- FeaturePlot(limbAdult,c("Ly6a","Pdgfra","Cd24a"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_pvssc_cell.pdf",width = 12,height = 6)


#== 9--------
Idents(limbAdult) <- limbAdult$anno_level_2
FeaturePlot(limbAdult,"Cdh2")
sc_cdh2_cell <- WhichCells(limbAdult,idents=c("LA_Diaphyseal MSC"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_cdh2_cell)
p2 <- FeaturePlot(limbAdult,c("Cdh2"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_cdh2_cell.pdf",width = 12,height = 6)


#== 10--------
Idents(limbAdult) <- limbAdult$anno_level_2
FeaturePlot(limbAdult,"Mcam")
sc_mcam_cell <- WhichCells(limbAdult,idents=c("LA_Pericyte"))
p1 <- DimPlot(limbAdult,cells.highlight = sc_mcam_cell)
p2 <- FeaturePlot(limbAdult,c("Mcam"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_mcam_cell.pdf",width = 12,height = 6)

#== 11---------------------------

Idents(limbAdult) <- limbAdult$anno_level_3
Idents(limbEmbryo) <- limbEmbryo$anno_level_3
DimPlot(limbEmbryo)
FeaturePlot(limbAdult,"Col10a1")
sc_col10a1_cell <- WhichCells(limbAdult,idents=c("LA_Hypertrophic.Chondrocyte"))
sc_col10a1_cell <- c(sc_col10a1_cell,WhichCells(limbEmbryo,idents=c("LE_Hypertrophic.Chondrocyte")))


p1 <- DimPlot(full_seurat,cells.highlight = sc_col10a1_cell)
p2 <- FeaturePlot(full_seurat,c("Col10a1"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_col10a1_cell.pdf",width = 10,height = 6)

#== 12---------------------------------------

Idents(limbAdult) <- limbAdult$anno_level_2
Idents(limbEmbryo) <- limbEmbryo$anno_level_2
DimPlot(limbEmbryo)
FeaturePlot(limbAdult,"Col2a1")
FeaturePlot(limbEmbryo,"Col2a1")
sc_col2a1_cell <- WhichCells(limbAdult,idents=c("LA_Chondrocyte"))
sc_col2a1_cell <- c(sc_col2a1_cell,WhichCells(limbEmbryo,idents=c("LE_Chondrocyte")))


p1 <- DimPlot(full_seurat,cells.highlight = sc_col2a1_cell)
p2 <- FeaturePlot(full_seurat,c("Col2a1"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_col2a1_cell.pdf",width = 10,height = 6)

#== 13-------------------------------


Idents(limbAdult) <- limbAdult$anno_level_2
Idents(limbEmbryo) <- limbEmbryo$anno_level_3
DimPlot(limbEmbryo)
FeaturePlot(limbEmbryo,"Dlx5")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
sc_dlx5_cell <- WhichCells(limbEmbryo,idents=c("LE_Trim47.Middle.Mesenchyme"))

p1 <- DimPlot(limbEmbryo,cells.highlight = sc_dlx5_cell)
p2 <- FeaturePlot(limbEmbryo,c("Dlx5"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_dlx5_cell.pdf",width = 10,height = 6)


#15--------------------------
FeaturePlot(limbEmbryo,"Cd200")
DimPlot(limbAdult)
FeaturePlot(limbAdult,"Cd200")

FeaturePlot(limbAdult,"Itgav")

sc_ssc_cell <- WhichCells(limbAdult,idents=c("LA_Serpina1c.Hypertrophic.Chondrocyte","LA_Cds1.Mature.Chondrocyte"))

p1 <- DimPlot(limbAdult,cells.highlight = sc_ssc_cell)
p2 <- FeaturePlot(limbAdult,c("Cd200"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_ssc_cell.pdf",width = 10,height = 6)


#16---------------------------------
Idents(limbAdult) <- limbAdult$anno_level_3
Idents(limbEmbryo) <- limbEmbryo$anno_level_4
Idents(head) <- head$anno_level_2
FeaturePlot(full_seurat,"Gli1")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
FeaturePlot(limbEmbryo,"Gli1")&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
VlnPlot(limbAdult,"Gli1")
VlnPlot(head,"Gli1")
VlnPlot(limbEmbryo,"Gli1")
dput(levels(limbEmbryo))
sc_gli_cell <- WhichCells(limbEmbryo,idents=c("LE_Trim47.Middle.Mesenchyme","LE_Pmaip1.Early.Mesenchyme","LE_Cytl1.Mature.Chondrocyte","LE_Nr2f1.Middle.Mesenchyme"))
sc_gli_cell <- c(sc_gli_cell,WhichCells(limbAdult,idents=c("LA_Metaphysis MSC","LA_Proliferating.Chondrocyte", 
  "LA_Chondrocyte progenitor", "LA_Mature.Chondrocyte")))
sc_gli_cell <- c(sc_gli_cell,WhichCells(head,idents=c("C_Mesenchyme","C_Sox18+ Dermis")))
p1 <- DimPlot(full_seurat,cells.highlight = sc_gli_cell)
p2 <- FeaturePlot(full_seurat,c("Gli1"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_gli1_cell.pdf",width = 10,height = 6)





#==17------------------------------------

FeaturePlot(limbAdult,"Pthlh")

Idents(limbEmbryo) <- limbEmbryo$anno_level_5
Idents(limbAdult) <- limbAdult$anno_level_4
VlnPlot(limbAdult,"Pthlh")
sc_pthlh_cell <- WhichCells(limbEmbryo,idents=c("LE_Irx1.Chondrocyte progenitor"))
sc_pthlh_cell <- c(sc_pthlh_cell,WhichCells(limbAdult,idents=c("LA_Pappa2.Proliferating.Chondrocyte")))

DimPlot(limbAdult)
p1 <- DimPlot(full_seurat,cells.highlight = sc_pthlh_cell)
p2 <- FeaturePlot(full_seurat,c("Pthlh"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_pthlh_cell.pdf",width = 10,height = 6)


#== 18-----------------------------------------

FeaturePlot(limbAdult,"Foxa2")
FeaturePlot(limbEmbryo,"Foxa2")
Idents(limbEmbryo) <- limbEmbryo$anno_level_5
Idents(limbAdult) <- limbAdult$anno_level_4
VlnPlot(limbAdult,"Foxa2")
VlnPlot(limbEmbryo,"Foxa2")
sc_foxa2_cell <- WhichCells(limbEmbryo,idents=c("LE_Aldh1a3.Hypertrophic.Chondrocyte","LE_Smpd3.Mature.Chondrocyte"))
sc_foxa2_cell <- c(sc_foxa2_cell,WhichCells(limbAdult,idents=c("LA_Serpina1c.Hypertrophic.Chondrocyte", 
                                                               "LA_Ihh.Mature.Chondrocyte", "LA_Cds1.Mature.Chondrocyte")))

DimPlot(limbAdult)
p1 <- DimPlot(full_seurat,cells.highlight = sc_foxa2_cell)
p2 <- FeaturePlot(full_seurat,c("Foxa2"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_foxa2_cell.pdf",width = 10,height = 6)

#== 19 -----------------------------
FeaturePlot(limbAdult,"Hmmr")
FeaturePlot(limbEmbryo,"Hmmr")
VlnPlot(limbAdult,"Hmmr")
Idents(limbEmbryo) <- limbEmbryo$anno_level_4
VlnPlot(limbEmbryo,"Hmmr")
DimPlot(limbAdult,group.by = "anno_level_4")
sc_hmmr_cell <- WhichCells(limbAdult,idents=c("LA_Pbk.Proliferating.Chondrocyte"))
sc_hmmr_cell <- c(sc_hmmr_cell,WhichCells(limbEmbryo,idents=c("LE_Trim47.Middle.Mesenchyme","LE_Nr2f1.Middle.Mesenchyme",
                                                              "LE_Pmaip1.Early.Mesenchyme"
                                                              )))
p1 <- DimPlot(full_seurat,cells.highlight = sc_hmmr_cell)
p2 <- FeaturePlot(full_seurat,c("Hmmr"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_hmmr_cell.pdf",width = 10,height = 6)

#==20---------------------------------
head <- full_seurat[,basic_boolean_head]
Idents(head) <- head$anno_level_4
FeaturePlot(head,"Axin2")


VlnPlot(head,"Axin2")
DimPlot(limbAdult,group.by = "anno_level_4")
sc_axin2_cell <- WhichCells(head,idents=c("C_Nusap1.Early mesenchyme","C_Lef1.Irx1.Mesenchyme"))
p1 <- DimPlot(head,cells.highlight = sc_axin2_cell)
p2 <- FeaturePlot(head,c("Axin2"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_axin2_cell.pdf",width = 10,height = 6)


#== 21-----------------

Idents(head) <- head$anno_level_3
FeaturePlot(head,"Hhip")


VlnPlot(head,"Hhip")
DimPlot(limbAdult,group.by = "anno_level_4")
sc_hhip_cell <- WhichCells(head,idents=c("C_Pre-osteoblast"))
p1 <- DimPlot(head,cells.highlight = sc_hhip_cell)
p2 <- FeaturePlot(head,c("Hhip"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_hhip_cell.pdf",width = 10,height = 6)

#== 22-------------------------

FeaturePlot(head,"Msx2")
DimPlot(head,label = T)

VlnPlot(head,"Msx2")
Idents(head) <- head$anno_level_4
sc_msx2_cell <- WhichCells(head,idents=c("C_Tnmd.Middle Mesenchyme","C_Matn4.Middle Mesenchyme","C_Tnn.Pre-osteoblast"))
Idents(head) <- head$anno_level_5
sc_msx2_cell <- c(sc_msx2_cell,WhichCells(head,idents=c("C_Tnn.Smpd3.Pre-osteoblast")))
p1 <- DimPlot(head,cells.highlight = sc_msx2_cell)
p2 <- FeaturePlot(head,c("Msx2"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_msx2_cell.pdf",width = 10,height = 6)


#== prrx1---------------------------
Idents(head) <- head$anno_level_2
sc_prrx1_cell <- WhichCells(head,idents=c( "C_Mesenchyme"))


#==23--------------------
limbEmbryo <- full_seurat[,full_seurat$Organ=="Limb_embryo"]
FeaturePlot(limbEmbryo,"Hes1")
FeaturePlot(limbEmbryo,"Sox9")
Idents(limbEmbryo) <- limbEmbryo$anno_level_3
DimPlot(limbEmbryo,label = T)
sc_hes1_cell <- WhichCells(limbEmbryo,idents=c("LE_Early.Mesenchyme"))
p1 <- DimPlot(limbEmbryo,cells.highlight = sc_hes1_cell)
p2 <- FeaturePlot(limbEmbryo,c("Hes1"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_hes1_cell.pdf",width = 10,height = 6)

#==24----------------------------------------------
unique(full_seurat$C49_named)
Idents(full_seurat) <- full_seurat$C49_named
sstr2_cell <-  WhichCells(full_seurat,idents=c( "Vcan+ MP BMSC"))
sstr2_cell <-  colnames(full_seurat)[full_seurat$C49_named=="Vcan+ MP BMSC"& full_seurat$Organ=="Limb_adult" ]
p1 <- DimPlot(full_seurat,cells.highlight = sstr2_cell)
p2 <- FeaturePlot(full_seurat,c("Sstr2"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))
p1+p2
ggsave("result/4.22_stem_cell_anno/sc_sstr2_cell.pdf",width = 10,height = 6)
#== 25----------------------------------

stemCell$`Sstr2+ MSC` <- sstr2_cell

totalCell <- unique(stemCell%>%unlist)
table(full_seurat$C2_named)["Non-osteo"]
length(totalCell)/table(full_seurat$C2_named)["Non-osteo"]

cellTable <- table(stemCell%>%unlist)%>%as.data.frame()
full_seurat$stemAnno <- 0

full_seurat$stemAnno[as.character(cellTable$Var1)] <- cellTable$Freq
FeaturePlot(full_seurat,"stemAnno")+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")),values = c(0,0.4,0.55,0.65,1.0))+
  ggtitle("Annotation times")
ggsave("result/4.22_stem_cell_anno/10.4_anno_times.pdf",width = 10,height = 8)
DimPlot(full_seurat,group.by = "stemAnno")

write.csv(full_seurat$stemAnno,"10.4_meta_stemAnno.csv")
colnames(full_seurat)[full_seurat$stemAnno==0]
table(full_seurat$mix_level_2[full_seurat$stemAnno==0])
sum(stemCell$`Lepr+ BMSC`%in%colnames(full_seurat)[full_seurat$stemAnno==0])
sscFreq <- table(full_seurat$stemAnno)%>%as.data.frame()
VlnPlot(full_seurat,features = "stemAnno",group.by = "C7_named")
ggplot(sscFreq, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "SSC Freq", x = "Freq", y = "Cell number")+
  theme_minimal() +
  scale_fill_manual(values  =  rev(brewer.pal(n = 8, name = "RdBu")))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), # Make title bold
        axis.text =  element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"), # Make x-axis label bold
        axis.title.y = element_text(face = "bold"))
