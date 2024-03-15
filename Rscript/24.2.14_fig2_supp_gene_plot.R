tab_bk <- readRDS("processed_data/12.23_genetab.Rds")
pb <- readRDS("processed_data/12.23_fig2_supp/pseudobulk.Rds")
library(dreamlet)
res.proc <- processAssays(pb, ~Age_num, min.count = 5)
plotVoom(res.proc)
vp.lst <- fitVarPart(res.proc, ~Age_num)
genes <- vp.lst$gene[2:4]
plotPercentBars(vp.lst[vp.lst$gene %in% genes, ])
plotVarPart(vp.lst, label.angle = 60)
res.dl <- dreamlet(res.proc, ~Age_num)
genes <- c("Apoe", "Bglap","Sp7","Fbn2","Cd200","Aspn","Postn","Cxcl12","Mfap4","Ovol2")
p <- plotGeneHeatmap(res.dl, coef = "Age_num", genes = genes)

tab <- tab_bk
tab$z.std[tab$adj.P.Val>0.05]=NA
tab <- tab[,-4]
df_wide <- pivot_wider(tab, names_from = assay, values_from = z.std)
df_wide <- df_wide[rowSums(is.na(df_wide))<4,]
df_wide[is.na(df_wide)] <- 0
same_rows <- apply(df_wide[2:6], 1, function(row) all(row >= 0,na.rm = T)| all(row <= 0,na.rm = T))
gene_up <- df_wide$ID[apply(df_wide[2:6], 1, function(row) all(row >= 0,na.rm = T))]
gene_down <- df_wide$ID[apply(df_wide[2:6], 1, function(row) all(row <= 0,na.rm = T))]
geneOther <- setdiff(df_wide$ID, c(gene_up,gene_down))



genes=c(as.character(gene_up[1:10]),as.character(gene_down[1:10]),as.character(geneOther[1:10]))
plotGeneHeatmap(res.dl, coef = "Age_num", genes = genes)
plotGeneHeatmap(res.dl, coef = "Age_num", genes = c("Foxq1","Suv39h1"))




res_zenith <- readRDS("processed_data/12.23_fig2_supp/MsigDb_zenit_C2.Rds")

library(zenith)
C2 <- get_MSigDB("C2",org = "mmu",to = "SYMBOL")
library(decoupleR)
gmt_to_decoupler <- function(pth) {
  pathways <- list()
  
  con <- file(pth, "r")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    fields <- strsplit(line, "\t")[[1]]
    name <- fields[[1]]
    genes <- fields[-1]
    pathways[[name]] <- genes
  }
  close(con)
  
  data <- do.call(rbind.data.frame, lapply(names(pathways), function(name) {
    data.frame(geneset = rep(name, length(pathways[[name]])), genesymbol = pathways[[name]], stringsAsFactors = FALSE)
  }))
  
  return(data)
}
test <- gmt_to_decoupler("data/m2.all.v2023.2.Mm.symbols (1).gmt")

geneSymbol <- test$genesymbol[test$geneset=="REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS"]
gene_up <- intersect(gene_up,geneSymbol)
gene_down <- intersect(gene_down,geneSymbol)
gene_other <- intersect(geneOther,geneSymbol)

genes <- c(gene_up[1:10],gene_down,gene_other)
plotGeneHeatmap(res.dl, coef = "Age_num", genes = genes)+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size = 15,face = "bold"),
                                                               axis.text.y = element_text(angle = 0, vjust = 1, hjust=1,size = 15,face = "bold"),
                                                               axis.title.y = element_text(hjust = 0.5, face="bold", size=15),
                                                               axis.title.x = element_text(hjust = 0.5, face="bold", size=15),
                                                               legend.text = element_text(face="bold", size=12))
ggsave("result/24.2.11_Fig2_supp/geneHM.pdf",width = 5,height = 9)
