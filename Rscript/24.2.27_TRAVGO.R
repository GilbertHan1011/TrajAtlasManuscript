RAVmodel <- readRDS("process_data/trajMap/")
valAll <- read.csv("process_data/trajMap/1.12_valAll_extend.csv",row.names = 1)
ravIndex=RAVindex(RAVextend)
colnames(ravIndex) <- rownames(valAll)

genelist_540 <- ravIndex[,"RAV_540"]%>%sort(decreasing = TRUE)%>%.[1:100]%>%names
genelist_766 <-ravIndex[,"RAV_766"]%>%sort(decreasing = TRUE)%>%.[1:100]%>%names
genelist_767 <-ravIndex[,"RAV_767"]%>%sort(decreasing = TRUE)%>%.[1:100]%>%names
genelist <- list(genelist_540,genelist_766,genelist_767)
names(genelist) <- c("TRAV540","TRAV766","TRAV777")
enrichGene <- clusterProfiler::compareCluster(genelist,fun = "enrichGO",keyType="SYMBOL",ont="BP",OrgDb="org.Mm.eg.db")
dotplot(enrichGene,showCategory = 4)+ theme(axis.text.x = element_text(size = 15,face = "bold"), axis.text.y = element_text(size = 15,face = "bold"))

saveRDS(enrichGene,"process_data/trajMap/2.27_TRAVGO.Rds")
