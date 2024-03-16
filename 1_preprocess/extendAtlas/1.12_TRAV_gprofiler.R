library(clusterProfiler)
colnames(RAVextend@assays@data$RAVindex) <- rownames(valAll)

genelist_540 <-RAVindex(RAVextend)[,"RAV_540"]%>%sort(decreasing = TRUE)%>%.[1:100]%>%names
genelist_766 <-RAVindex(RAVextend)[,"RAV_766"]%>%sort(decreasing = TRUE)%>%.[1:100]%>%names
genelist_767 <-RAVindex(RAVextend)[,"RAV_767"]%>%sort(decreasing = TRUE)%>%.[1:100]%>%names
genelist <- list(genelist_540,genelist_766,genelist_767)
names(genelist) <- c("RAV_540","RAV_766","RAV_767")
enrichGene <- clusterProfiler::compareCluster(genelist,fun = "enrichGO",keyType="SYMBOL",ont="BP",OrgDb="org.Mm.eg.db")
dotplot(enrichGene,showCategory = 4)+ theme(axis.text.x = element_text(size = 15,face = "bold"), axis.text.y = element_text(size = 15,face = "bold"))
ggsave("result/1.13_volcano/dotplot.pdf",width = 8,height = 7)
