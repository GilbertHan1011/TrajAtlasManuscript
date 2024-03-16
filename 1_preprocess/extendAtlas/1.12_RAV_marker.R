setwd("../../8.25_full_integration/")
library(MuDataSeurat)
trajMap <- ReadH5MU("process_data/trajMap/trajMap_1.h5mu")
DefaultAssay(trajMap) <- "TRAV"
trajMes <- trajMap[,trajMap$Lineage=="Mesenchyme"]
InjuryMarker <- FindMarkers(trajMes,group.by = "MesInjury",ident.1 = "MesInjury",logfc.threshold = 0.01)
InjuryMarker3 <- FindMarkers(trajMes,group.by = "Stage",ident.1 = "Injury(Regeneration)",logfc.threshold = 0.01)
InjuryMarker2 <- FindMarkers(trajMap,group.by = "Stage",ident.1 = "Injury(Regeneration)",logfc.threshold = 0.01)
View(avg_load)


trajMap2 <- ReadH5MU("process_data/trajMap/1.12_trajmap_2nd.h5mu")
DefaultAssay(trajMap2) <- "TRAV"
trajMes <- trajMap2[,trajMap2$Lineage=="Mesenchyme"]
InjuryMarker <- FindMarkers(trajMes,group.by = "MesInjury",ident.1 = "MesInjury",logfc.threshold = 0.01)
InjuryMarker3 <- FindMarkers(trajMes,group.by = "Stage",ident.1 = "Injury(Regeneration)",logfc.threshold = 0.01)
InjuryMarker2 <- FindMarkers(trajMap,group.by = "Stage",ident.1 = "Injury(Regeneration)",logfc.threshold = 0.01)
View(avg_load)


plotData <- InjuryMarker3%>%
  rownames_to_column("gene")%>%
  mutate(logP=-log10(p_val))%>%
  mutate(type = ifelse((logP>2&avg_log2FC > 5), "UP",
                       ifelse((logP>2&avg_log2FC < -5),"DOWN","NONE")))


g <- ggplot(plotData, aes(x = avg_log2FC, y = logP, color = type))+ 
  geom_point(size = 2)+
  scale_color_manual(values=c('#0000FF', 'grey','red'))+ 
  theme_classic()+scale_x_continuous(name="log2FoldChange")+
  scale_y_continuous(name="-log10Pvalue") + 
  labs(title = "Injury_vs_Non-injury")+
  theme(plot.title = element_text(hjust=0.5,size=16, face="bold"))+
  geom_vline(xintercept = 5, linetype="dashed")+ geom_vline(xintercept = -5, linetype="dashed")+
  geom_hline(yintercept = 2, linetype="dashed")+
  theme(plot.title = element_text(hjust=0.5,size=16, face="bold"))+
  theme(text = element_text(size= 20))+
  theme(axis.line = element_line(size=1, colour = "black"))+
  theme(axis.text.x = element_text(size = 15,face = "bold"), axis.text.y = element_text(size = 15,face = "bold"))+
  theme(legend.position="none")

g <- g+geom_label_repel(data=dplyr::filter(plotData, logP>4&avg_log2FC > 5), aes(label=gene),box.padding = unit(0.9, "lines"),
                  segment.color = "black",segment.size = 0.5,size=5)
g
ggsave("result/1.13_volcano/1.13_volcano.pdf",width = 6,height = 6)
