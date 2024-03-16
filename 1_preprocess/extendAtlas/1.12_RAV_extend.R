#== function------------------------------

selectGO <- function(RAVmodel, gsea.dir) {
  
  ## If you want to select only a subset of RAVs with the specific cluster size
  # ind <- which(metadata(RAVmodel)$size > 3)
  # gsea_all <- vector(mode = "list", length = length(ind))
  # names(gsea_all) <- colnames(RAVmodel)[ind]
  
  gsea_all <- vector(mode = "list", length = ncol(RAVmodel))
  names(gsea_all) <- colnames(RAVmodel)
  gsea.dir <- gsea.dir
  
  for (i in seq_len(ncol(RAVmodel))) {
    pathToRes <- file.path(gsea.dir, paste0("go_", i, ".rds"))
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

#==read data--------------


RAVNormal <- readRDS("../3.9_wt_integrate/processed_data/12.12_TRAV/12.26_ravmodel_nmf2.Rds")
trainingPC <- readRDS("../3.9_wt_integrate/processed_data/12.12_TRAV/trainingData_PCclusters_60.Rds")
normalNMF <- readRDS("../important_processed_data/24.1.trainingdataNMF.Rds")
normalLoad <- readRDS("../important_processed_data/24.1.core_avgload.Rds")
varGene <- readRDS("../3.9_wt_integrate/processed_data/12.12_TRAV/topgene_var2000.Rds")
annot_database <- "MSigDB C5" 
#== merge-----------------------
allZ_list_normal <- lapply(normalNMF, function(x) x$rotation)
allZ_list_normal <- lapply(allZ_list_normal, function(x) {
  x <- as.data.frame(x)
  x <- x[varGene,]
  rownames(x) <- varGene
  x[is.na(x)] <- 0
  return(x)
})

for (i in 1:length(allZ_list_normal)){
  colnames(allZ_list_normal[[i]]) <- paste0(names(allZ_list_normal)[i],"_",colnames(allZ_list_normal[[i]]))
}
allZNomal <- Reduce(cbind, allZ_list_normal)
allNomal <- t(allZNomal) 
RAVassay1=RAVNormal@assays@data$RAVindex
RAVassay2=RAVmodel@assays@data$RAVindex
training1 <- DataFrame(allNomal)
training2 <- RAVmodel@trainingData



RAVextend <-  PCAGenomicSignatures(assays = list(RAVindex =cbind(RAVassay1,RAVassay2)),trainingData =rbind(training1,training2))
RAVextend@metadata$cluster <- c(RAVNormal@metadata$cluster,RAVmodel@metadata$cluster)
RAVextend@metadata$size <- c(RAVNormal@metadata$size,RAVmodel@metadata$size)
RAVextend@metadata$k <- RAVNormal@metadata$k+RAVmodel@metadata$k
RAVextend@metadata$n <- 15
geneSets(RAVextend) <- annot_database

studyCombine <- c(trainingPC$studies,studies(RAVmodel))
studies(RAVextend) <- studyCombine

saveRDS(RAVextend,"../important_processed_data/24.2.27_RAV_extend.Rds")


#== go-----------------

term2gene <- clusterProfiler::read.gmt("../8.25_full_integration/data/m5.all.v2023.2.Mm.symbols.gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(RAVextend))) {
  fname <- paste0("go_", i, ".rds")
  fpath <- file.path("process_data/trajMap/TRAV_GO/", "M5", fname)
  
  geneList <- RAVindex(RAVextend)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- names(geneList)[geneList>0]
  k = min(length(geneList),100)
  res <- clusterProfiler::enricher(geneList[1:k], TERM2GENE = term2gene,
                                   pvalueCutoff = 0.05)
  saveRDS(res, fpath)
}



gsea_m5_var <- selectGO(RAVextend, "process_data/trajMap/TRAV_GO/M5/")
gsea(RAVextend) <- gsea_m5_var

#== RAV extend------------------


val_list <- list()
# rav_test <- readRDS("processed_data/12.12_TRAV/Ravmodel_30_scale_3000gene.Rds")


for (index in 1:length(filterFileNameSub)) {
  x <- data.table::fread(filterFileNameSub[index], 
                         showProgress = FALSE,header = T)
  study=shortNameSub[index]
  # x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
  x <- data.frame(x[,-1], row.names = x$V1)
  topGeneInstert <- intersect(varGene,rownames(x))
  x <- x[topGeneInstert,,drop=FALSE]
  #rownames(x) <- topGene48
  #x[is.na(x)] <- 0
  x <- t(scale(t(x)))
  #x[is.na(x)] <- 0
  # Normalization
  #x <- sweep(x, 1, m)
  #x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
  tryCatch({
    valide_df <- validate_local(x, RAVextend)
    val_list[[study]] <- valide_val
    valide_val <- valide_df$score
  },error = function(err) {
    cat("An error occurred in iteration", index, ":", conditionMessage(err), "\n")
    errorIndex <<- c(errorIndex,index)
    # Skip to the next iteration
    return(invisible(NULL))
  })
  
  # PCA
  rm(x)
  
}
valAll <- data.frame(val_list)
RAVName <- paste0("RAV_",1:dim(valAll)[1])
rownames(valAll) <- RAVName

write.csv(valAll,"process_data/trajMap/1.12_valAll_extend.csv")

#== 3 -----------------


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
