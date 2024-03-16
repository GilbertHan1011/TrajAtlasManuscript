## Use similiarity function to prove organ atlas can be integrated

#=== data import and function--------------------
library(ggsci)
rm(cranio)
varGene=read.csv("../3.9_wt_integrate/7.6_software/trajMap/src/trajMap/datasets/variable_2000.csv")
varGene=varGene$X0%>%unlist()
full_sub=full_seurat[varGene,]
library(randomForest)
downsample_and_predict <- function(reflatent_seurat, scANVI, subsample_size) {
  require(randomForest)
  cell_ranks <- levels(reflatent_seurat)
  downsampled_cells <- c()
  for (i in 1:length(cell_ranks)) {
    current_cells <- WhichCells(reflatent_seurat, idents = cell_ranks[i])
    if (min(length(current_cells), subsample_size) < subsample_size) {
      warning(paste0("The current ident ", cell_ranks[i], " is smaller than sample size, adjust to the ident size"))
    }
    downsampled_cells <- c(downsampled_cells, sample(current_cells, size = min(length(current_cells), subsample_size), replace = FALSE))
  }
  downsample_seurat <- reflatent_seurat[, downsampled_cells]
  if (scANVI!="X"){
    rf_model <- randomForest(downsample_seurat@reductions[[scANVI]]@cell.embeddings, y = Idents(downsample_seurat))
  } else{
    rf_model <- randomForest(t(as.matrix(downsample_seurat@assays$originalexp@data)), y = Idents(downsample_seurat))
  }

  votes <- rf_model$votes
  identmatrix <- as.data.frame(Idents(downsample_seurat))
  names(identmatrix) <- "ident"
  votes <- cbind(votes, identmatrix)
  summary <- votes %>%
    group_by(ident) %>%
    summarize(across( everything(),.fns = mean))%>%
    column_to_rownames("ident")%>%
    as.matrix()
  summary <- apply(summary, 2, function(x) (x)/max(x))
  diag(summary) <- 0
  summary <- t(summary)
  return(summary)
}
write.csv(sum,"../3.9_wt_integrate/result/9.24_integrate_fundation/9.24_rf_predict_similiarity.csv")
Idents(full_sub) <- "anno_level_2"
sum=downsample_and_predict(full_sub,scANVI = "X",subsample_size = 1000)
pheatmap(sum,cluster_cols = F,cluster_rows = F)
cellCranio=c("C_Mesenchyme","C_Chondrocyte", "C_Osteoblast", "C_Sox18+ Dermis",  "C_Hypodermis")
cellLA=c(  "LA_Fibroblast",   "LA_Chondrocyte","LA_Osteoblast",
           "LA_Pericyte","LA_Diaphyseal MSC")
cellLE=c( "LE_Mesenchyme", "LE_Chondrocyte", "LE_Osteoblast", 
           "LE_Smooth muscle cell")
max_color_value=0.15
pdf("../3.9_wt_integrate/result/9.24_integrate_fundation/9.24_la_C.pdf")
pheatmap(sum[cellCranio,cellLA],cluster_cols = F,cluster_rows = F,
         breaks = seq(0, max_color_value, length.out = 101))
dev.off()
pdf("../3.9_wt_integrate/result/9.24_integrate_fundation/9.24_c_le.pdf")
pheatmap(sum[cellCranio,cellLE],cluster_cols = F,cluster_rows = F,
         breaks = seq(0, max_color_value, length.out = 101))  
dev.off()
pdf("../3.9_wt_integrate/result/9.24_integrate_fundation/9.24_la_le.pdf")
pheatmap(sum[cellLA,cellLE],cluster_cols = F,cluster_rows = F,
         breaks = seq(0, max_color_value, length.out = 101))
dev.off()

sankyMeta=full_sub@meta.data[c("Organ","anno_level_2")]
sankyMetaC=sankyMeta$anno_level_2[sankyMeta$Organ=="Head"]
sankyMetaLA=sankyMeta$anno_level_2[sankyMeta$Organ=="Limb_adult"]
sankyMetaLE=sankyMeta$anno_level_2[sankyMeta$Organ=="Limb_embryo"]

sankyTableC=table(sankyMetaC)%>%as.data.frame()
sankyTableC$Freq=(sankyTableC$Freq/sum(sankyTableC$Freq))
sankyTableLA=table(sankyMetaLA)%>%as.data.frame()
sankyTableLA$Freq=(sankyTableLA$Freq/sum(sankyTableLA$Freq))
sankyTableLE=table(sankyMetaLE)%>%as.data.frame()
sankyTableLE$Freq=(sankyTableLE$Freq/sum(sankyTableLE$Freq))

sankyTableC$common=c("Chondro","Hypodermis","Mes","Ob","Dermis")
sankyTableLA$common=c("Chondro","Lepr","Fibro","Ob","smc")
sankyTableLE$common=c("Chondro","Mes","Ob","smc")

sankyTableC$organ="Head"
sankyTableLE$organ="LE"
sankyTableLA$organ="LA"
names(sankyTableC) <- c("sankyMeta", "Freq", "common", "organ")
names(sankyTableLA) <- c("sankyMeta", "Freq", "common", "organ")
names(sankyTableLE) <- c("sankyMeta", "Freq", "common", "organ")
plotData=rbind(sankyTableC,sankyTableLE,sankyTableLA)

plotData$organ <- factor(plotData$organ,c("Head", "LE", "LA"))

plotData$sankyMeta <- gsub("C_","",plotData$sankyMeta)%>%gsub("LE_","",.)
plotData$sankyMeta <- gsub("LA_","",plotData$sankyMeta)
plotData$sankyMeta <- factor(plotData$sankyMeta,c("Chondrocyte",  "Osteoblast", "Mesenchyme",
  "Smooth muscle cell","Pericyte","Sox18+ Dermis",  "Hypodermis", "Diaphyseal MSC", "Fibroblast"
))

plotData$common <- gsub("Chondro","Chondrocyte",plotData$common)
plotData$common <- gsub("Ob","Osteoblast",plotData$common)
plotData$common <- gsub("Mes","Mesenchyme",plotData$common)
plotData$common <- gsub("smc","Pericyte",plotData$common)

plotData$common <- factor(plotData$common,c("Chondrocyte","Osteoblast", "Pericyte","Mesenchyme",   "Lepr", "Hypodermis", "Dermis",
                                     "Fibro"))
gg <- ggplot(data = plotData,
             aes(x = organ, y = Freq,stratum=common, alluvium = common,fill=sankyMeta,label=sankyMeta)) +
  theme_bw() +
  scale_fill_adaptive(name = "npg", palette = "nrc")

# proportional knot positioning (default)
gg +
  geom_alluvium(aes(fill = common),
                alpha = .75,  width = 1/2) +
  geom_stratum(aes(stratum = common),  width = 1/2)+
  ggfittext::geom_fit_text(stat = "stratum",min.size = 0.1) +theme(legend.position = "none")+
  scale_x_discrete(labels=c('Head', 'Limb bud', 'Long bone'))+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16)) 
ggsave("../3.9_wt_integrate/result/9.24_integrate_fundation/sanky.pdf",width = 10,height = 6)



sankyTableC=table(sankyMetaC)%>%as.data.frame()

sankyTableLA=table(sankyMetaLA)%>%as.data.frame()

sankyTableLE=table(sankyMetaLE)%>%as.data.frame()
sankyTableC$common=c("Chondro","Hypodermis","Mes","Ob","Dermis")
sankyTableLA$common=c("Chondro","Lepr","Fibro","Ob","smc")
sankyTableLE$common=c("Chondro","Mes","Ob","smc")

sankyTableC$organ="Head"
sankyTableLE$organ="LE"
sankyTableLA$organ="LA"
names(sankyTableC) <- c("sankyMeta", "Freq", "common", "organ")
names(sankyTableLA) <- c("sankyMeta", "Freq", "common", "organ")
names(sankyTableLE) <- c("sankyMeta", "Freq", "common", "organ")
sankydata=rbind(sankyTableC,sankyTableLE,sankyTableLA)
plotData2=plotData
plotData2$Freq <- sankydata$Freq
gg <- ggplot(data = plotData2,
             aes(x = organ, y = Freq,stratum=common, alluvium = common,fill=sankyMeta,label=sankyMeta)) +
  theme_bw() +
  scale_fill_adaptive(name = "npg", palette = "nrc")

# proportional knot positioning (default)
gg +
  geom_alluvium(aes(fill = common),
                alpha = .75,  width = 1/2) +
  geom_stratum(aes(stratum = common),  width = 1/2)+
  ggfittext::geom_fit_text(stat = "stratum",min.size = 0.1) +theme(legend.position = "none")+
  scale_x_discrete(labels=c('Head', 'Limb bud', 'Long bone'))+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16)) 
ggsave("../3.9_wt_integrate/result/9.24_integrate_fundation/sanky_orig.pdf",width = 10,height = 6)

sum(plotData$Freq)