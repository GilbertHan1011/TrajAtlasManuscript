library(clusterProfiler)
library(ggstance)
source("../function/enrichFun.R")

library(ReactomePA)


GC_cluster <- read.csv("result/6.21_traj_hm/12.26_symbolOrderCluster.csv",row.names = 1)
geneSplit <-  split(GC_cluster$symbolOrdered, GC_cluster$X.1)

m2GMT <- read.gmt("data/m2.all.v2023.2.Mm.symbols (1).gmt")
term <- unique(m2GMT$term)
reactomeTerm <- grep("REACTOME",term)
reactomeTerm <- term[reactomeTerm]
reactomeGMT <- m2GMT[m2GMT$term%in%reactomeTerm,]

goEnrichReactome <- lapply(geneSplit,FUN = function(x){ enricher(x,TERM2GENE=reactomeGMT)})

Res1 <- makeEnrichTable(goEnrichReactome,2,c(1:5))
ggplot(Res1, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,5)),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")
ggsave("result/24.1.20_trajdiff_dotplot/Res1_bar.pdf",width = 5,height = 3)

for (i in 1:12){
  name=paste0("GC_",i)
  saveName=paste0("result/24.2.19_suppGeneHMEnrich/",name,".pdf")
  Res1 <- makeEnrichTable(goEnrichReactome,name,c(1:5))
  ggplot(Res1, aes(x = logP, y = ID,fill=logP)) +
    geom_barh(stat = "identity", color = "black") +
    geom_text(aes(label = ID,y=ID,x=rep(0,5)),hjust=0) +theme_minimal()+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
    scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")+ggtitle(name)
  ggsave(saveName,width = 5,height = 3)
}
View(goEnrichReactome[["GC_11"]]@result)

Res1 <- makeEnrichTable(goEnrichReactome,"GC_11",c(2,3,5,6,12))
ggplot(Res1, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,5)),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")+
  ggtitle("test")
ggsave("result/24.2.19_suppGeneHMEnrich/GC11.pdf",width = 5,height = 3)

#== preOb---------------------------------
mscUp <- preobMarker%>%
  rownames_to_column("gene")%>%
  arrange(desc(avg_log2FC))%>%top_n(500,wt = avg_log2FC)%>%select(gene)%>%unlist()
mscUpEnrich <- enricher(mscUp,TERM2GENE=reactomeGMT)
mscDown <- preobMarker%>%
  rownames_to_column("gene")%>%
  arrange(desc(avg_log2FC))%>%top_n(500,wt = desc(avg_log2FC))%>%select(gene)%>%unlist()
mscDownEnrich <- enricher(mscDown,TERM2GENE=reactomeGMT)



mesEnrich <- list(mscUpEnrich,mscDownEnrich)
Res1 <- makeEnrichTable(mesEnrich,1,c(1,2,4,6,21))
ggplot(Res1, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,5)),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")+
  ggtitle("MSC-OPCST Preosteoblast")
ggsave("result/24.2.19_suppGeneHMEnrich/msc_opcst_preob.pdf",width = 5,height = 3)
Res2 <- makeEnrichTable(mesEnrich,2,c(1,2,4,6,15))
ggplot(Res2, aes(x = logP, y = ID,fill=logP)) +
  geom_barh(stat = "identity", color = "black") +
  geom_text(aes(label = ID,y=ID,x=rep(0,5)),hjust=0) +theme_minimal()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth =1)) +
  scale_fill_gradient(low = "LightSkyBlue", high = "DeepSkyBlue3")+
  ggtitle("Lepr+ BMSC-OPCST Preosteoblast")
ggsave("result/24.2.19_suppGeneHMEnrich/lepr_opcst_preob.pdf",width = 5,height = 3)



wt_integrate$OPCST <- NA
wt_integrate$OPCST[rownames(lineage)[lineage$lineage_mesenchyme=="True"]] <- "MSC-OPCST"
wt_integrate$OPCST[rownames(lineage)[lineage$lineage_lepr=="True"]] <- "lepr+ BMSC-OPCST"
wt_integrate$OPCST[wt_integrate$C19_named!="Pre-ob"] <- NA
DimPlot(wt_integrate,group.by = "OPCST")
SCpubr::do_DimPlot(wt_integrate,group.by = "OPCST",reduction = "X_draw_graph_fa",
                   label = T,raster = T,raster.dpi = 600,pt.size = 3,na.value = "grey90")
ggsave("result/24.2.19_suppGeneHMEnrich/preOB_umap.pdf",width = 6,height = 6)
