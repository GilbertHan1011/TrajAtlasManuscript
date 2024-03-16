
selectGO <- function(RAVmodel, gsea.dir) {
  
  ## If you want to select only a subset of RAVs with the specific cluster size
  # ind <- which(metadata(RAVmodel)$size > 3)
  # gsea_all <- vector(mode = "list", length = length(ind))
  # names(gsea_all) <- colnames(RAVmodel)[ind]
  
  gsea_all <- vector(mode = "list", length = ncol(RAVmodel))
  names(gsea_all) <- colnames(RAVmodel)
  gsea.dir <- gsea.dir
  
  for (i in seq_len(ncol(RAVmodel))) {
    pathToRes <- file.path(gsea.dir, paste0("gsea_", i, ".rds"))
    res <- readRDS(pathToRes)
    res=res@result
    # If there is no enriched pathways
    if (nrow(res) == 0) {
      resName <- paste0("RAV", i)
      gsea_all[[resName]] <- res[, c("Description", "GeneRatio", "pvalue", "qvalue","Count"), drop = FALSE]
      print(paste("RAV", i, "has no enriched pathways."))
      next
    }
    
    res <- res[, c("Description", "GeneRatio", "pvalue", "qvalue","Count"), drop = FALSE]
    #resName <- paste0("RAV", i)
    resName=colnames(RAVmodel)[i]
    gsea_all[[resName]] <- res
    print(paste("RAV", i, "is added."))
  }
  
  return(gsea_all)
}


#### 5. RAVmodel from 536 studies C2 ###########################################

# if (!dir.exists(out_dir)) {dir.create(out_dir)}

library(GenomicSuperSignature)
library(clusterProfiler)

#== GSEA-----------------------------

# MSigDB C2
term2gene <- clusterProfiler::read.gmt("data/m2.all.v2023.2.Mm.symbols (1).gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(RAVmodel))) {
  fname <- paste0("gsea_", i, ".rds")
  fpath <- file.path("processed_data/12.12_TRAV/geneset/", "gsea_var_nmf", fname)
  
  geneList <- RAVindex(RAVmodel)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  res <- clusterProfiler::GSEA(geneList, TERM2GENE = term2gene,
                               pvalueCutoff = 0.05, seed = TRUE)
  saveRDS(res, fpath)
}


# test1 <- readRDS("processed_data/12.12_TRAV/geneset/gsea_var_nmf/gsea_1.rds")
# test2 <- enricher(names(geneList)[1:100], TERM2GENE=term2gene)




# MSigDB C2
term2gene <- clusterProfiler::read.gmt("data/m2.all.v2023.2.Mm.symbols (1).gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(RAVmodel))) {
  fname <- paste0("gsea_", i, ".rds")
  fpath <- file.path("processed_data/12.12_TRAV/geneset/", "gsea_var_nmf", fname)
  
  geneList <- RAVindex(RAVmodel)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- names(geneList)[geneList>0]
  k = min(length(geneList),100)
  res <- clusterProfiler::enricher(geneList[1:k], TERM2GENE = term2gene,
                               pvalueCutoff = 0.05)
  saveRDS(res, fpath)
}


# test1 <- readRDS("processed_data/12.12_TRAV/geneset/gsea_var_nmf/gsea_1.rds")
# test2 <- enricher(names(geneList)[1:100], TERM2GENE=term2gene)


#== GO-----------------------------------
for (i in seq_len(ncol(RAVmodel))) {
  fname <- paste0("go_", i, ".rds")
  fpath <- file.path("processed_data/12.12_TRAV/geneset/", "gsea_go_var", fname)
  
  geneList <- RAVindex(RAVmodel)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- names(geneList)[geneList>0]
  k = min(length(geneList),100)
  res <- clusterProfiler::enrichGO(geneList[1:k], keyType = "SYMBOL",ont="BP",OrgDb = "org.Mm.eg.db",
                                   pvalueCutoff = 0.05)
  saveRDS(res, fpath)
}


#=============================
# MSigDB C3
term2gene <- clusterProfiler::read.gmt("data/m3.gtrd.v2023.2.Mm.symbols.gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(RAVmodel))) {
  fname <- paste0("gsea_", i, ".rds")
  fpath <- file.path("processed_data/12.12_TRAV/geneset/go_C3_tf/", fname)
  
  geneList <- RAVindex(RAVmodel)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- names(geneList)[geneList>0]
  k = min(length(geneList),100)
  res <- clusterProfiler::enricher(geneList[1:k], TERM2GENE = term2gene,
                                   pvalueCutoff = 0.05)
  saveRDS(res, fpath)
}
