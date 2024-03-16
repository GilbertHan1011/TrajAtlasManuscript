library(zellkonverter)
library(tidyverse)
seurat <- readH5AD("../important_processed_data/4.17_wt_integrate_SCRAN_log1p_scANVI.h5ad")
full_seurat <- as.Seurat(seurat,counts = "counts",data = "X")

seurat_meta <- readRDS("../important_processed_data/4.20_wt_metadata.Rds")
full_seurat@meta.data <- seurat_meta
rm(seurat)

DimPlot(full_seurat,group.by = "anno_level_1")
Idents(full_seurat) <- full_seurat$anno_level_1
sample_level_1 <- downsample_and_predict(full_seurat,scANVI = "X_scANVI",subsample_size = 1000)
sample_level_1_NA <- sample_level_1 
sample_level_1_NA[sample_level_1_NA<0.2]=0
pheatmap::pheatmap(sample_level_1_NA)
merge_id_level1 <- c("Mix-Non-osteo cells", "Mix-Osteoblast", "Mix-Non-osteo cells", 
                     "Mix-Osteoblast", "Mix-Osteoblast", "Mix-Non-osteo cells")
names(merge_id_level1) <- levels(full_seurat)  
full_seurat <- RenameIdents(full_seurat,merge_id_level1)
full_seurat$mix_level_1 <- Idents(full_seurat)

#== level2----------------------
Idents(full_seurat) <- full_seurat$anno_level_2
sample_level_2 <- downsample_and_predict(full_seurat,scANVI = "X_scANVI",subsample_size = 1000,threshold = 0.18)
pheatmap::pheatmap(sample_level_2)
merge_id_level2 <- c("Mix_Mesenchyme", "C_Hypodermis", "Mix_Osteoblast", "C_Sox18+ Dermis", 
                     "Mix_Chondrocyte", "Mix_Chondrocyte", "Mix_Mesenchyme", "Mix_Osteoblast", 
                     "Mix_Pericyte", "Mix_Osteoblast", "LA_Fibroblast", "LA_Diaphyseal MSC", 
                     "Mix_Pericyte", "Mix_Chondrocyte")
names(merge_id_level2) <- levels(full_seurat)  
full_seurat <- RenameIdents(full_seurat,merge_id_level2)
full_seurat$mix_level_2 <- Idents(full_seurat)
DimPlot(full_seurat)
#== level3----------------------------------
Idents(full_seurat) <- full_seurat$anno_level_3
sample_level_3 <- downsample_and_predict(full_seurat,scANVI = "X_scANVI",subsample_size = 700,threshold = 0.18)
pheatmap::pheatmap(sample_level_3)
merge_id_level3 <-c("Mix_Pre-osteoblast", "C_Hypodermis", "Mix_Osteoblast", "C_Meninges", 
                    "C_Late Mesenchyme", "C_Middle Mesenchyme", "Mix_Early.Mesenchyme", 
                    "Mix_Irx1.Mesenchyme", "C_Sox18+ Dermis", "Mix_Mature.Chondrocyte", "Mix_Hypertrophic.Chondrocyte", 
                    "LE_Middle.Mesenchyme", "Mix_Osteoblast", "LE_Late.Mesenchyme", 
                    "Mix_Early.Mesenchyme", "Mix_Pericyte", "Mix_Mature.Chondrocyte", 
                    "LE_Chondrocyte progenitor", "Mix_Irx1.Mesenchyme", "Mix_Osteoblast", 
                    "LA_Aspn+ MSC", "LA_Lepr+ MSC", "Mix_Pericyte", "Mix_Pre-osteoblast", 
                    "LA_Dpt.Fibroblast", "LA_Ppbp+ MSC", "LA_A2m.Chondrocyte", "LA_Angptl7.Fibroblast", 
                    "Mix_Hypertrophic.Chondrocyte", "LA_Proliferating.Chondrocyte", 
                    "LA_Chondrocyte progenitor", "Mix_Mature.Chondrocyte", "LA_Serpina1a.Chondrocyte"
)
names(merge_id_level3) <- levels(full_seurat)  
full_seurat <- RenameIdents(full_seurat,merge_id_level3)
full_seurat$merge_id_level3 <- Idents(full_seurat)
DimPlot(full_seurat)

test1 <- downsample_and_predict(full_seurat,scANVI = "X_scANVI",subsample_size = 700,threshold = 0.18)
pheatmap::pheatmap(test1)


#== level4----------------------------------
Idents(full_seurat) <- full_seurat$anno_level_4
sample_level_4 <- downsample_and_predict(full_seurat,scANVI = "X_scANVI",subsample_size = 500,threshold = 0.18)
pheatmap::pheatmap(sample_level_4)
full_seurat$Organ <- full_seurat$batch_atlas



merge_id_level4 <-c("C_Tnn.Pre-osteoblast", "Mix_Ly6a_Mesenchyme", "C_Smpd3.Pre-osteoblast", 
                    "Mix_Ifitm5.Osteoblast", "C_Meninges", "Mix_Mmp13.Pre-osteoblast", "C_Late Mesenchyme", 
                    "C_Matn4.Middle Mesenchyme", "C_Nusap1.Early mesenchyme", "C_Barx1.Middle Mesenchyme", 
                    "C_Gjb6.Middle Mesenchyme", "C_Tnmd.Middle Mesenchyme", "Mix_Irx1.Mesenchyme", 
                    "C_Barx1.Early mesenchyme", "Mix_Ogn.Late.Mesenchyme", "C_Sox18+ Dermis", 
                    "C_Mmp9.Chondrocyte", "Mix_Matn3.Chondrocyte", "Mix_Hmmr.Chondrocyte progenitor", 
                    "LE_Hypertrophic.Chondrocyte", "LE_Nr2f1.Middle.Mesenchyme", 
                    "Mix_Ifitm5.Osteoblast", "Mix_Ogn.Late.Mesenchyme", "LE_Tubb5.Early.Mesenchyme", 
                    "LE_Gsc.Late.Mesenchyme", "LE_Col1a1.Chondrocyte", "Mix_Matn3.Chondrocyte", 
                    "LE_Pmaip1.Early.Mesenchyme", "LE_Chondrocyte progenitor", "LE_Smpd3.Mature.Chondrocyte", 
                    "Mix_Ly6a_Mesenchyme", "LE_Cenpa.Late.Mesenchyme", "LE_Cilp2.Pre-osteoblast", 
                    "LE_Foxp2.Late.Mesenchyme", "Mix_Pericyte", "Mix_Hmmr.Chondrocyte progenitor", 
                    "LE_Tnmd.Mature.Chondrocyte", "LE_Cxcl14.Late.Mesenchyme", "Mix_Irx1.Mesenchyme", 
                    "LE_Rpgrip1.Early.Mesenchyme", "LE_Cyp26a1.Mature.Chondrocyte", 
                    "LA_Lipc.Osteoblast", "LA_Tnn.Aspn+ MSC", "LA_Esm1.Lepr+ MSC", 
                    "LA_Ccl2.Lepr+ MSC", "Mix_Pericyte", "LA_Angpt4.Metaphysis MSC", 
                    "Mix_Mmp13.Pre-osteoblast", "Mix_Ifitm5.Osteoblast", "LA_Tspan11.Dpt.Fibroblast", 
                    "LA_Kcnk2.Lepr+ MSC", "LA_Ppbp+ MSC", "LA_Kit.Lepr+ MSC", "LA_A2m.Chondrocyte", 
                    "Mix_Ly6a_Mesenchyme", "LA_Tnxb.Angptl7.Fibroblast", "LA_Chodl.Angptl7.Fibroblast", 
                    "LA_Ffar4.Hypertrophic.Chondrocyte", "LA_Tspan15.Dpt.Fibroblast", 
                    "Mix_Hmmr.Chondrocyte progenitor", "LA_Cytl1.Metaphysis MSC", 
                    "LA_Mmp3.Aspn+ MSC", "LA_Chondrocyte progenitor", "LA_Pappa2.Proliferating.Chondrocyte", 
                    "LA_Hs6st2.Proliferating.Chondrocyte", "LA_R3hdml.Proliferating.Chondrocyte", 
                    "LA_Il17b.Mature.Chondrocyte", "LA_Serpina1c.Hypertrophic.Chondrocyte", 
                    "LA_Ihh.Mature.Chondrocyte", "LA_Cds1.Mature.Chondrocyte", "LA_Serpina1a.Chondrocyte"
)
cranio <- full_seurat[,full_seurat$Organ=="Head"]
la <- full_seurat[,full_seurat$Organ=="Limb_adult"]
le <- full_seurat[,full_seurat$Organ=="Limb_embryo"]
DimPlot(le)
FeaturePlot(la,"Satb2")
FeaturePlot(cranio,"Ly6a")
FeaturePlot(le,"Spp1")
test <- subset(full_seurat,idents=c("LE_Col4a1.Smooth muscle cell","LA_Pericyte","LE_Spp1.Smooth muscle cell"))
DimPlot(test)
VlnPlot(test,c("Col4a1","Acta2"))
merge_submit(sample_level_4,0.25)
pheatmap(sample_level_4)
names(merge_id_level4) <- levels(full_seurat)  
full_seurat <- RenameIdents(full_seurat,merge_id_level4)
full_seurat$merge_id_level4 <- Idents(full_seurat)
le$anno_level_3
DimPlot(full_seurat)
DimPlot(le,group.by = "C22")
le$C21_named
test1 <- downsample_and_predict(full_seurat,scANVI = "X_scANVI",subsample_size = 700,threshold = 0.18)
pheatmap::pheatmap(test1)


#== level 5--------------------------
Idents(full_seurat) <- full_seurat$anno_level_5
sample_level_5 <- downsample_and_predict(full_seurat,scANVI = "X_scANVI",subsample_size = 500,threshold = 0.25)
merge_submit(sample_level_5,0.4)
merge_id_level5 <- c("C_Tnn.Pre-osteoblast", "Mix_Ly6a_Mesenchyme", "C_Ifitm5.Smpd3.Pre-osteoblast", 
                     "Mix_Ifitm5.Osteoblast", "C_Meninges", "C_Cp.Mmp13.Pre-osteoblast", 
                     "C_C1qtnf3.Late Mesenchyme", "C_Gpha2.Matn4.Middle Mesenchyme", 
                     "C_Pax3.Nusap1.Early mesenchyme", "C_Tnn.Smpd3.Pre-osteoblast", 
                     "C_Cenpf.Nusap1.Early mesenchyme", "C_Ostn.Mmp13.Pre-osteoblast", 
                     "C_Mecom.Barx1.Middle Mesenchyme", "C_Tac1.Late Mesenchyme", 
                     "C_Slc6a13.Gjb6.Middle Mesenchyme", "Mix_Hmmr.Chondrocyte progenitor", 
                     "Mix_Lipc.Osteoblast", "C_Hoxd1.Lef1.Irx1.Mesenchyme", "C_Fxyd5.Gjb6.Middle Mesenchyme", 
                     "C_Aldh1a2.Gjb6.Middle Mesenchyme", "C_Epyc.Matn4.Middle Mesenchyme", 
                     "C_Dusp9.Barx1.Early mesenchyme", "C_Agtr2.Irx1.Mesenchyme", 
                     "C_Tbx1.Nusap1.Early mesenchyme", "C_Vtn.Barx1.Early mesenchyme", 
                     "C_Sox18+ Dermis", "C_Nusap1.Lef1.Irx1.Mesenchyme", "C_Dlx1.Barx1.Early mesenchyme", 
                     "C_Sox18.Lef1.Irx1.Mesenchyme", "C_Hells.Barx1.Early mesenchyme", 
                     "C_Mmp9.Chondrocyte", "C_Pax3.Lef1.Irx1.Mesenchyme", "C_Heyl.Barx1.Early mesenchyme", 
                     "Mix_Scrg1.Mature.Chondrocyte", "C_Cytl1.Matn3.Chondrocyte", "Mix_Hmmr.Chondrocyte progenitor", 
                     "C_Mzb1.Lef1.Irx1.Mesenchyme", "C_Pax9.Barx1.Middle Mesenchyme", 
                     "Mix_Smpd3.Mature.Chondrocyte", "LE_Steap4.Hypertrophic.Chondrocyte", 
                     "LE_Cenpf.Nr2f1.Middle.Mesenchyme", "LE_Aldh1a3.Hypertrophic.Chondrocyte", 
                     "LE_Hp.Satb2.Osteoblast", "LE_Crtac1.Ogn.Late.Mesenchyme", "LE_Tubb5.Early.Mesenchyme", 
                     "LE_Gsc.Late.Mesenchyme", "LE_Col1a1.Chondrocyte", "Mix_Scrg1.Mature.Chondrocyte", 
                     "LE_Nusap1.Nr2f1.Middle.Mesenchyme", "LE_Smpd3.Cytl1.Mature.Chondrocyte", 
                     "LE_Lhx9.Pmaip1.Early.Mesenchyme", "LE_Irx1.Chondrocyte progenitor", 
                     "LE_Gdf5.Ogn.Late.Mesenchyme", "LE_Ebf1.Chondrocyte progenitor", 
                     "Mix_Smpd3.Mature.Chondrocyte", "LE_Hist1h1b.Nr2f1.Middle.Mesenchyme", 
                     "Mix_Ly6a_Mesenchyme", "LE_Ccnb2.Cenpa.Late.Mesenchyme", 
                     "Mix_Ifitm5.Osteoblast", "LE_Gdf5.Nr2f1.Middle.Mesenchyme", 
                     "LE_Cilp2.Pre-osteoblast", "LE_Dlx5.Foxp2.Late.Mesenchyme", "Mix_Pericyte", 
                     "Mix_Hmmr.Chondrocyte progenitor", "LE_Tbx18.Ogn.Late.Mesenchyme", 
                     "LE_Tnmd.Mature.Chondrocyte", "LE_Coch.Ogn.Late.Mesenchyme", 
                     "LE_Ung.Pmaip1.Early.Mesenchyme", "LE_Hoxd12.Cxcl14.Late.Mesenchyme", 
                     "LE_Nrgn.Irx1.Mesenchyme", "LE_Osr1.Pmaip1.Early.Mesenchyme", 
                     "LE_Thy1.Cxcl14.Late.Mesenchyme", "LE_Hoxc5.Ogn.Late.Mesenchyme", 
                     "LE_Dkk2.Ogn.Late.Mesenchyme", "LE_Rpgrip1.Early.Mesenchyme", 
                     "LE_Gria2.Foxp2.Late.Mesenchyme", "LE_Tbx18.Foxp2.Late.Mesenchyme", 
                     "LE_Msx1.Cenpa.Late.Mesenchyme", "LE_Cyp26a1.Mature.Chondrocyte", 
                     "LE_Dlx5.Trim47.Middle.Mesenchyme", "LE_Il33.Irx1.Mesenchyme", 
                     "Mix_Lipc.Osteoblast", "LA_Tnn.Aspn+ MSC", "LA_Esm1.Lepr+ MSC", 
                     "LA_Ccl2.Lepr+ MSC", "Mix_Pericyte", "LA_Angpt4.Metaphysis MSC", 
                     "LA_Vcan.Metaphysis MSC", "Mix_Ifitm5.Osteoblast", "LA_Tspan11.Dpt.Fibroblast", 
                     "LA_Kcnk2.Lepr+ MSC", "LA_Ppbp+ MSC", "LA_Kit.Lepr+ MSC", "LA_A2m.Chondrocyte", 
                     "Mix_Ly6a_Mesenchyme", "LA_Tnxb.Angptl7.Fibroblast", "LA_Chodl.Angptl7.Fibroblast", 
                     "LA_Ffar4.Hypertrophic.Chondrocyte", "LA_Tspan15.Dpt.Fibroblast", 
                     "LA_Pbk.Proliferating.Chondrocyte", "LA_Cytl1.Metaphysis MSC", 
                     "LA_Mmp3.Aspn+ MSC", "LA_Chondrocyte progenitor", "LA_Pappa2.Proliferating.Chondrocyte", 
                     "LA_Hs6st2.Proliferating.Chondrocyte", "LA_R3hdml.Proliferating.Chondrocyte", 
                     "LA_Il17b.Mature.Chondrocyte", "LA_Serpina1c.Hypertrophic.Chondrocyte", 
                     "LA_Ihh.Mature.Chondrocyte", "Mix_Smpd3.Mature.Chondrocyte", "LA_Serpina1a.Chondrocyte"
)
test <- subset(full_seurat,idents=c("LE_Thy1.Cxcl14.Late.Mesenchyme","C_Hoxd1.Lef1.Irx1.Mesenchyme"))
DimPlot(test)
VlnPlot(test,c("Thy1","Hoxd1"))
names(merge_id_level5) <- levels(full_seurat)  
full_seurat <- RenameIdents(full_seurat,merge_id_level5)
full_seurat$merge_id_level5 <- Idents(full_seurat)

df <- full_seurat@meta.data[c("anno_level_2","mix_level_2")]
df_test <- df%>%
  group_by(anno_level_2)%>%
  unique()%>%column_to_rownames("anno_level_2")
for (i in rownames(df_test)){
  if(df_test[i,]!=i){
    print(paste0("Change ",i, " to ",df_test[i,]))
  }
  else{
    df_test <- subset(df_test, row.names(df_test) != i)
  }
}
anno_update_1 <- show_ident_change(full_seurat,previous_anno = "anno_level_1",latter_anno = "mix_level_1")
anno_update_2 <- show_ident_change(full_seurat,previous_anno = "anno_level_2",latter_anno = "mix_level_2")
anno_update_3 <-  show_ident_change(full_seurat,previous_anno = "anno_level_3",latter_anno = "merge_id_level3")
anno_update_4 <-  show_ident_change(full_seurat,previous_anno = "anno_level_4",latter_anno = "merge_id_level4")
anno_update_5 <-  show_ident_change(full_seurat,previous_anno = "anno_level_5",latter_anno = "merge_id_level5")
write.csv(anno_update_1,"result/4.20_similarity/level1_update.csv")
write.csv(anno_update_2,"result/4.20_similarity/level2_update.csv")
write.csv(anno_update_3,"result/4.20_similarity/level3_update.csv")
write.csv(anno_update_4,"result/4.20_similarity/level4_update.csv")
write.csv(anno_update_5,"result/4.20_similarity/level5_update.csv")

saveRDS(full_seurat@meta.data,"../important_processed_data/4.28_wt_metadata.Rds")
previous_anno="anno_level_2"
latter_anno = "mix_level_2"
anno_update_2 <- full_seurat@meta.data[c("anno_level_2","mix_level_2")]
anno_update_3 <- full_seurat@meta.data[c("anno_level_3","mix_level_3")]
