makeEnrichTable <- function(enrichRes,index,goindex){
  Res1 <- enrichRes[[index]]@result%>%
    mutate(logP=-log10(p.adjust))%>%
    dplyr::select(ID,logP)%>%
    mutate(ID = gsub("REACTOME_", "", ID))
  Res1Sub <- Res1[goindex,]
  Res1Sub$ID <- gsub("REACTOME_","",Res1Sub$ID)
  Res1Sub <- Res1Sub %>%
    mutate(ID = reorder(ID, logP))
  return(Res1Sub)
}