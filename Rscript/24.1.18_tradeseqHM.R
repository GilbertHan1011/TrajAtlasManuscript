library(Seurat)
myColor=mypal
names(myColor) <- c("Ob","Fibroblast", "MSC",  "Chondro","Lepr+ BMSC")
scMeta=read.csv("../important_processed_data/5.22_wt_integrate_meta.csv",row.names = 1)
wt_intergrate$C7_named2 <- scMeta$C7_named
wt_intergrate$C19_named2 <- scMeta$C19_named
DimPlot(wt_intergrate,group.by = c("C7_named","C7_named2","C19_named"))

mes_meta3 <- wt_intergrate@meta.data[c("C7_named2")]%>%as.data.frame()
mes_meta3 <- mes_meta3[rownames(pseduo),]
mes_meta3 <- merge(pseduo,mes_meta3,by=0,all=TRUE)%>%column_to_rownames("Row.names")
mes_meta <- as.data.frame(pseduo)
mes_meta$C7_named <- mes_meta3
line1_meta <- na.exclude(mes_meta[c(1,5)])%>%filter(.[[1]]>0)
line2_meta <- na.exclude(mes_meta[c(2,5)])%>%filter(.[[1]]>0)
line3_meta <- na.exclude(mes_meta[c(3,5)])%>%filter(.[[1]]>0)
line4_meta <- na.exclude(mes_meta[c(4,5)])%>%filter(.[[1]]>0)

#update function----------------
library(Hmisc) 
makeHa <- function(meta){
  meta_plot <-  meta%>%
    mutate(rnk=floor(.[[1]]/(max(.[[1]])+0.001) * 50)+1)
  tabMetaIdent<-table(meta_plot[c("C7_named","rnk")])%>%as.matrix()
  maxVal <- apply(tabMetaIdent,2,function(x) names(which.max(x)))
  return(list(maxVal))
}
line1_max <- makeHa(line1_meta)
line2_max <- makeHa(line2_meta)
line3_max <- makeHa(line3_meta)
line4_max <- makeHa(line4_meta)

line2_max[[1]] <- c(rep("Fibroblast",17),line2_max[[1]][17:49])
identHa <- c(line1_max[[1]],line2_max[[1]],line3_max[[1]],line4_max[[1]])

haSingle = HeatmapAnnotation(
  level_1 = factor(identHa,levels = unique(identHa)),
  col = list(
    level_1=myColor
  )
)
hmAllOrder2 <- ComplexHeatmap::Heatmap(plotDataAllOrder,cluster_columns = F,show_column_names = F,cluster_rows = F,
                                      show_row_names = F,col=colorRamp2(c(-2, 0, 3), c("DeepSkyBlue3", "white", "red")),
                                      row_split = rowSplit,column_split = columnSplit,top_annotation=haSingle,use_raster = T,border = T)
pdf("result/6.21_traj_hm/1.18_traj_hm_raster2.pdf",width = 10,height = 12)
hmAllOrder2 <- draw(hmAllOrder2)
dev.off()
saveRDS(plotDataAllOrder,"result/6.18_traj_map/24.1.19hmData.Rds")
