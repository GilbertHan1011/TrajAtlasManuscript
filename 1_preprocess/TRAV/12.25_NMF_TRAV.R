
setwd('~/limb/3.9_wt_integrate/')

#== run with NMF----------------------
trainingData_NMF=list()
shortName <- c(paste0("mesenchyme_",mesLineageName),paste0("lepr_",leprLineageName),
               paste0("chondro_",chondroLineageName),paste0("fibro_",fibroLineageName))
shortName2 <- c(paste0(mesLineageName,"_sep_Mesenchyme"),paste0(leprLineageName,"_sep_Lepr_BMSC"),
               paste0(chondroLineageName,"_sep_Chondro"),paste0(fibroLineageName,"_sep_Fibroblast"))



n <- 15 #== base on cross-validation
index <- 1
for (index in 1:length(filterName)) {
  x <- data.table::fread(filterName[index], 
                         showProgress = FALSE,header = T)
  study=shortName[index]
  # x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
  x <- data.frame(x[,-1], row.names = x$V1)
  
  topGeneInstert <- intersect(varGene,rownames(x))
  x <- x[topGeneInstert,,drop=FALSE]
  #rownames(x) <- topGene48
  #x[is.na(x)] <- 0
  x <- t(scale(t(x)))
  x[is.na(x)] <- 0
  # Normalization
  #x <- sweep(x, 1, m)
  #x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
  
  model_test <- RcppML::nmf(t(x), n, verbose = F, seed = 1234)
  w <- model_test$w
  d <- model_test$d
  h <- model_test$h
  colnames(h) <- rownames(x)
  trainingData_NMF[[study]]$rotation <- t(h)
  # colnames(trainingData_PCA[[study]]$rotation) <- paste0(study, ".PC", c(1:n))
  # eigs <-d
  # pca_summary <- rbind(SD = sqrt(eigs),
  #                      Variance = eigs/sum(eigs),
  #                      Cumulative = cumsum(eigs)/sum(eigs))
  # trainingData_PCA[[study]]$variance <- pca_summary[,1:n]
  #colnames(trainingData_NMF[[study]]$variance) <- paste0(study, ".PC", c(1:n))
  
  rm(x)
}
#str(trainingData_PCA$CalvariaP4_Ayturk)
#saveRDS(trainingData_PCA,"processed_data/12.12_TRAV/whole_PCA.Rds")

# cluster------------
# d <- 4
# x <- test
# test <- allZ_list$mesenchyme_CalvariaP4_Ayturk%>%as.data.frame()
allZ_list <- lapply(trainingData_NMF, function(x) x$rotation)
allZ_list <- lapply(allZ_list, function(x) {
  x <- as.data.frame(x)
  x <- x[varGene,]
  rownames(x) <- varGene
  x[is.na(x)] <- 0
  return(x)
})
for (i in 1:length(allZ_list)){
  colnames(allZ_list[[i]]) <- paste0(names(allZ_list)[i],"_",colnames(allZ_list[[i]]))
}
allZ <- Reduce(cbind, allZ_list)
all <- t(allZ) 
all <- all[rowSums(all)!=0,]
#rownames(all) <- 1:dim(all)[1]
print(paste("Dimension of allZ is", dim(allZ)))
res.dist <- factoextra::get_dist(all, method = "spearman")
ClusterKnot=8
clusterNum=floor(dim(allZ)[2]/ClusterKnot)
res.hcut <- factoextra::hcut(res.dist, k = clusterNum, hc_func = "hclust", 
                             hc_method = "ward.D", hc_metric = "spearman")


trainingData_PCclusters_60 <- buildAvgLoading(allZ, clusterNum, cluster = res.hcut$cluster)
avg_load <- trainingData_PCclusters_60$avgLoading

#== 30-----------------

cl <- trainingData_PCclusters_60$cluster
silh_res <- cluster::silhouette(cl, res.dist)
cl_silh_width <- summary(silh_res)$clus.avg.widths
trainingData_PCclusters_60$sw <- cl_silh_width  # add silhouette width to the result

# pca_summary <- list()
# for (i in seq_along(trainingData_NMF)) {
#   pca_summary[[i]] <- trainingData_NMF[[i]]$variance
#   names(pca_summary)[i] <- names(trainingData_NMF)[i]
# }

annot_database <- "MSigDB C2"   # will be added into the RAVmodel metadata
train_data=trainingData_PCclusters_60$avgLoading
colnames(train_data) <- names(trainingData_PCclusters_60$studies)
RAVmodel <- PCAGenomicSignatures(assays = list(RAVindex = as.matrix(trainingData_PCclusters_60$avgLoading)),trainingData =DataFrame(train_data))
RAVmodel <- PCAGenomicSignatures(trainingData = as.data.frame(trainingData_PCclusters_60$avgLoading))
DataFrame(trainingData_PCclusters_60[c("cluster")])
metadata(RAVmodel) <- trainingData_PCclusters_60[c("cluster", "size", "k", "n")]
names(metadata(RAVmodel)$size) <- paste0("RAV", seq_len(ncol(RAVmodel)))
geneSets(RAVmodel) <- annot_database
studies(RAVmodel) <- trainingData_PCclusters_60$studies
silhouetteWidth(RAVmodel) <- trainingData_PCclusters_60$sw
#metadata(RAVmodel)$MeSH_freq <- MeSH_freq
trainingData(RAVmodel)$PCAsummary <- pca_summary[rownames(trainingData(RAVmodel))]
#mesh(RAVmodel) <- trainingData_MeSH[rownames(trainingData(RAVmodel))]
updateNote(RAVmodel) <- note

#== GSEA------------------
dir2 <- system.file("scripts", package = "GenomicSuperSignature")
searchPathways_func <- file.path(dir2, "searchPathways.R")
source(searchPathways_func)  # load the function
# gsea_all <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_h/")
# gsea_m5 <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_m5/")
gsea_nmf_3000 <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_nmf_3000/")


debug(selectGO)
gsea_m5_var <- selectGO(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_var_nmf/")
gsea(RAVmodel) <- gsea_m5_var
descript <- lapply(gsea_m5_var, function(x) x["Description"])
descriptDf <-  plyr::ldply(descript, rbind)

# 
# lapply(gsea_nmf_3000[1:20], function(x) x["Description"])
saveRDS(RAVmodel,"processed_data/12.12_TRAV/12.26_ravmodel_nmf2.Rds")
RAVmodel
saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_30_scale.Rds")
testModel <- readRDS("processed_data/12.12_TRAV/12.26_ravmodel_nmf.Rds")

# 
# 
# #== top------------------------------------
# # Function to select top 100 values' indices
# select_top_100_indices <- function(column) {
#   top_100_indices <- order(column, decreasing = TRUE)[1:200]
#   return(top_100_indices)
# }
# 
# # Apply the function to each column
# top_100_indices_by_column <- apply(avg_load, 2, select_top_100_indices)
# test <- apply(top_100_indices_by_column, 2, function(x) rownames(avg_load)[x])
# #geneList <-  lapply(test, unlist)
# geneList <-  apply(test, 2,as.list)
# geneList <- lapply(geneList,unlist)
# 
# write.csv(test,"processed_data/12.12_TRAV/nmf_60_top100.csv")
# 
# test1 <- clusterProfiler::enrichGO(test[,1],OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
# ck <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
# saveRDS(ck,"processed_data/12.12_TRAV/NMF_GO.Rds")
# 
# 
# 
# 
# 
# #=== retrive validate score--------------------
# val_list <- list()
# rav_test <- readRDS("processed_data/12.12_TRAV/Ravmodel_30_scale_3000gene.Rds")
# 
# 
# for (index in 1:length(filterName)) {
#   x <- data.table::fread(filterName[index], 
#                          showProgress = FALSE,header = T)
#   study=shortName[index]
#   # x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
#   x <- data.frame(x[,-1], row.names = x$V1)
#   topGeneInstert <- intersect(topGene48,rownames(x))
#   x <- x[topGeneInstert,,drop=FALSE]
#   #rownames(x) <- topGene48
#   #x[is.na(x)] <- 0
#   x <- t(scale(t(x)))
#   #x[is.na(x)] <- 0
#   # Normalization
#   #x <- sweep(x, 1, m)
#   #x <- sweep(x, 1, s, "/")
#   # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
#   valide_df <- validate(x, rav_test)
#   valide_val <- valide_df$score
#   # PCA
#   rm(x)
#   
#   val_list[[study]] <- valide_val
# }
# valAll <- data.frame(val_list)
# 
# scale_val <- t(scale(t(valAll)))
# 
# pcaVal = princomp(t(scale_val))
# fviz_pca_biplot(pcaVal) +
#   ggtitle("")
# 
# 
# rownames(scale_val) <- rownames(valide_df)
# val_Seurat <- CreateSeuratObject(scale_val)
# valReduction <- CreateDimReducObject(embeddings = as.matrix(pcaVal$scores), key = "Val_", assay = DefaultAssay(val_Seurat))
# val_Seurat@reductions$pca <- valReduction
# DimPlot(val_Seurat)
# FeaturePlot(val_Seurat,"RAV1")
# 

#== readGO-----------------


val_list <- list()
# rav_test <- readRDS("processed_data/12.12_TRAV/Ravmodel_30_scale_3000gene.Rds")


for (index in 1:length(filterName)) {
  x <- data.table::fread(filterName[index], 
                         showProgress = FALSE,header = T)
  study=shortName[index]
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
  valide_df <- validate(x, RAVmodel)
  valide_val <- valide_df$score
  # PCA
  rm(x)
  
  val_list[[study]] <- valide_val
}
valAll <- data.frame(val_list)

scale_val <- t(scale(t(valAll)))

pcaVal = princomp(scale_val)
fviz_pca_biplot(pcaVal) +
  ggtitle("")


rownames(scale_val) <- rownames(valide_df)
val_Seurat <- CreateSeuratObject(scale_val)
#valReduction <- CreateDimReducObject(embeddings = as.matrix(pcaVal$), key = "Val_", assay = DefaultAssay(val_Seurat))
val_Seurat@reductions$pca <- valReduction
DimPlot(val_Seurat)
FeaturePlot(val_Seurat,"RAV1")
val_Seurat <- ScaleData(val_Seurat)
val_Seurat <- RunPCA(val_Seurat,features = rownames(val_Seurat))

DimPlot(val_Seurat)

FeaturePlot(val_Seurat,"Cl226-108 (30/30)")

colnames(valAll) <- shortName
colnames(valAll) <- shortName2
rownames(valAll) <- paste0("RAV", seq_len(ncol(RAVmodel)))
write.csv(valAll,"processed_data/12.12_TRAV/12.26_valide_origin_nmf.csv")


val_Seurat@assays$RNA@data <- as.matrix(valAll)

test3 <- val_Seurat@assays$RNA@scale.data


saveRDS(trainingData_NMF,"../important_processed_data/24.1.trainingdataNMF.Rds")
