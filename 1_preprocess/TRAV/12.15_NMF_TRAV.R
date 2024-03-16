cutoff <- 0.3  # Adjust this if you want different cutoff

for (i in seq_len(length(yList))){
  k <- length(yList[[i]])
  k <- round(k*cutoff)
  topGenes <- names(yList[[i]])[1:k]
  topGeneList[[i]] <- topGenes
}

topGene <- topGeneList%>%unlist%>%table%>%as.data.frame()
topGene84 <- topGene%>%filter(Freq>60)%>%dplyr::select(".")%>%unlist()%>%as.character()


sigGene <- read.csv("../important_processed_data/10.29_rf_importance_gene.csv")
sigGene <- sigGene%>%arrange("importance")%>%
  top_n(2000)%>%
  dplyr::select(gene)%>%unlist()

varGene <- rownames(wt_integrate)[wt_integrate@assays$originalexp@meta.features$highly_variable]


rbGene <- varGene[grepl("^Rp[sl]",varGene)]
mtGene <-  varGene[grepl("^mt-",varGene)]
HbbGene <-  varGene[grepl("^Hb",varGene)]
varGene <- setdiff(varGene,rbGene)%>%setdiff(.,mtGene)%>%setdiff(.,HbbGene)
# 
# topGenes <- c()
# 
# k <- length(yList[[1]])
# k <- round(k*cutoff)
# topGenes <- names(yList[[1]])[1:k]
# 
# 
# 
# for (i in 2:length(yList)) {
#   l <- length(yList[[i]])
#   l <- round(l*cutoff)
#   ls <- list(topGenes, names(yList[[i]])[1:l])
#   topGenes <- Reduce(intersect, ls)
# }
# 
# topGenes <- topGenes[nzchar(topGenes)] 
saveRDS(topGenes,"processed_data/12.12_TRAV/mes_topgene.Rds")
saveRDS(topGene48,"processed_data/12.12_TRAV/topgene_48.Rds")
saveRDS(varGene,"processed_data/12.12_TRAV/topgene_var2000.Rds")
trainingData_PCA=list()
#== Run PCA-------------------------------
n <- 20
index <- 1
for (index in 1:length(filterName)) {
  x <- data.table::fread(filterName[index], 
                         showProgress = FALSE,header = T)
  study=shortName[index]
  # x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
  x <- data.frame(x[,-1], row.names = x$V1)
  
  topGeneInstert <- intersect(sigGene,rownames(x))
  x <- x[topGeneInstert,,drop=FALSE]
  #rownames(x) <- topGene48
  #x[is.na(x)] <- 0
  x <- t(scale(t(x)))
  x[is.na(x)] <- 0
  # Normalization
  #x <- sweep(x, 1, m)
  #x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
  
  # PCA
  pca_res <- prcomp(t(x))   # x is a matrix with genes(row) x samples(column)
  trainingData_PCA[[study]]$rotation <- pca_res$rotation[,1:20]
  colnames(trainingData_PCA[[study]]$rotation) <- paste0(study, ".PC", c(1:n))
  eigs <- pca_res$sdev^2
  pca_summary <- rbind(SD = sqrt(eigs),
                       Variance = eigs/sum(eigs),
                       Cumulative = cumsum(eigs)/sum(eigs))
  trainingData_PCA[[study]]$variance <- pca_summary[,1:n]
  colnames(trainingData_PCA[[study]]$variance) <- paste0(study, ".PC", c(1:n))
  
  rm(x)
}
#str(trainingData_PCA$CalvariaP4_Ayturk)
#saveRDS(trainingData_PCA,"processed_data/12.12_TRAV/whole_PCA.Rds")

# cluster------------
d <- 4
x <- test
test <- allZ_list$mesenchyme_CalvariaP4_Ayturk%>%as.data.frame()
allZ_list <- lapply(trainingData_PCA, function(x) x$rotation)
allZ_list <- lapply(allZ_list, function(x) {
  x <- as.data.frame(x)
  x <- x[sigGene,]
  rownames(x) <- sigGene
  x[is.na(x)] <- 0
  return(x)
})
allZ <- Reduce(cbind, allZ_list)
all <- t(allZ) 
print(paste("Dimension of allZ is", dim(allZ)))
res.dist <- factoextra::get_dist(all, method = "spearman")
# 
# nbMax <- fviz_nbclust(as.matrix(res.dist), kmeans, method = "wss",k.max = 100)
# #fviz_nbclust(as.matrix(res.dist), kmeans, method = "gap_stat",k.max = 200)
# nbMax+ geom_vline(xintercept = 35, linetype = 2)
# k <- round(nrow(all)/d, 0)
# start <- Sys.time()
# res.hcut <- factoextra::hcut(res.dist, k = k, hc_func = "hclust", 
#                              hc_method = "ward.D", hc_metric = "spearman")
# end <- Sys.time()
# d <- 4
# k <- round(nrow(all)/d, 0)
# start <- Sys.time()
# res.hcut_605 <- factoextra::hcut(res.dist, k = k, hc_func = "hclust", 
#                                  hc_method = "ward.D", hc_metric = "spearman")
res.hcut_30 <- factoextra::hcut(res.dist, k = 30, hc_func = "hclust", 
                                hc_method = "ward.D", hc_metric = "spearman")
# fviz_dend(res.hcut_30, rect = TRUE, k_colors = 30, cex = 0.7, main = "Hierarchical Clustering Dendrogram")
# 
# 
# res.hcut_10 <- factoextra::hcut(res.dist, k = 10, hc_func = "hclust", 
#                                 hc_method = "ward.D", hc_metric = "spearman")
# fviz_dend(res.hcut_10, rect = TRUE, k_colors = 10,show_labels = F, cex = 0.7, main = "Hierarchical Clustering Dendrogram")

# determine optimized cluster number
#fviz_nbclust(my_data, kmeans, method = "gap_stat")

# fviz_nbclust(as.matrix(res.dist), kmeans, method = "gap_stat",k.max = 200)
# 
# trainingData_PCclusters_605 <- buildAvgLoading(allZ, k, cluster = res.hcut_605$cluster)

trainingData_PCclusters_30 <- buildAvgLoading(allZ, 30, cluster = res.hcut_30$cluster)
avg_load <- trainingData_PCclusters_30$avgLoading
#== 30-----------------

cl <- trainingData_PCclusters_30$cluster
silh_res <- cluster::silhouette(cl, res.dist)
cl_silh_width <- summary(silh_res)$clus.avg.widths
trainingData_PCclusters_30$sw <- cl_silh_width  # add silhouette width to the result




pca_summary <- list()
for (i in seq_along(trainingData_PCA)) {
  pca_summary[[i]] <- trainingData_PCA[[i]]$variance
  names(pca_summary)[i] <- names(trainingData_PCA)[i]
}
# 
# installed_pkgs <- installed.packages()
# installed_pkgs["GenomicSuperSignaturePaper", ]

annot_database <- "MSigDB C5"   # will be added into the RAVmodel metadata
RAVmodel <- PCAGenomicSignatures(assays = list(RAVindex = as.matrix(trainingData_PCclusters_30$avgLoading)))
metadata(RAVmodel) <- trainingData_PCclusters_30[c("cluster", "size", "k", "n")]
names(metadata(RAVmodel)$size) <- paste0("RAV", seq_len(ncol(RAVmodel)))
geneSets(RAVmodel) <- annot_database
studies(RAVmodel) <- trainingData_PCclusters_30$studies
silhouetteWidth(RAVmodel) <- trainingData_PCclusters_30$sw
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
gsea_m5_all_30 <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_m5_3000/")
lapply(gsea_m5_all_30[1], function(x) x["Description"])
saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_30.Rds")
RAVmodel
saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_30_scale_3000gene.Rds")




#== run with NMF----------------------
trainingData_NMF=list()

n <- 20
index <- 1
for (index in 1:length(filterName)) {
  x <- data.table::fread(filterName[index], 
                         showProgress = FALSE,header = T)
  study=shortName[index]
  # x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
  x <- data.frame(x[,-1], row.names = x$V1)
  
  topGeneInstert <- intersect(sigGene,rownames(x))
  x <- x[topGeneInstert,,drop=FALSE]
  #rownames(x) <- topGene48
  #x[is.na(x)] <- 0
  x <- t(scale(t(x)))
  x[is.na(x)] <- 0
  # Normalization
  #x <- sweep(x, 1, m)
  #x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
  
  model_test <- RcppML::nmf(t(x), 30, verbose = F, seed = 1234)
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
  x <- x[sigGene,]
  rownames(x) <- sigGene
  x[is.na(x)] <- 0
  return(x)
})
allZ <- Reduce(cbind, allZ_list)
all <- t(allZ) 
all <- all[rowSums(all)!=0,]
rownames(all) <- 1:dim(all)[1]
print(paste("Dimension of allZ is", dim(allZ)))
res.dist <- factoextra::get_dist(all, method = "spearman")

res.hcut <- factoextra::hcut(res.dist, k = 60, hc_func = "hclust", 
                                hc_method = "ward.D", hc_metric = "spearman")


trainingData_PCclusters_60 <- buildAvgLoading(allZ, 60, cluster = res.hcut$cluster)
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

annot_database <- "MSigDB C5"   # will be added into the RAVmodel metadata
RAVmodel <- PCAGenomicSignatures(assays = list(RAVindex = as.matrix(trainingData_PCclusters_60$avgLoading)))
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

lapply(gsea_nmf_3000[1:20], function(x) x["Description"])
saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_30.Rds")
RAVmodel
saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_30_scale.Rds")



#== top------------------------------------
# Function to select top 100 values' indices
select_top_100_indices <- function(column) {
  top_100_indices <- order(column, decreasing = TRUE)[1:200]
  return(top_100_indices)
}

# Apply the function to each column
top_100_indices_by_column <- apply(avg_load, 2, select_top_100_indices)
test <- apply(top_100_indices_by_column, 2, function(x) rownames(avg_load)[x])
#geneList <-  lapply(test, unlist)
geneList <-  apply(test, 2,as.list)
geneList <- lapply(geneList,unlist)

write.csv(test,"processed_data/12.12_TRAV/nmf_60_top100.csv")

test1 <- clusterProfiler::enrichGO(test[,1],OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
ck <- compareCluster(geneCluster = geneList, fun = enrichGO,OrgDb = "org.Mm.eg.db",ont="ALL",keyType="SYMBOL")
saveRDS(ck,"processed_data/12.12_TRAV/NMF_GO.Rds")





#=== retrive validate score--------------------
val_list <- list()
rav_test <- readRDS("processed_data/12.12_TRAV/Ravmodel_30_scale_3000gene.Rds")


for (index in 1:length(filterName)) {
  x <- data.table::fread(filterName[index], 
                         showProgress = FALSE,header = T)
  study=shortName[index]
  # x <- data.frame(x[,-1], row.names = x$V1) %>% as.matrix
  x <- data.frame(x[,-1], row.names = x$V1)
  topGeneInstert <- intersect(topGene48,rownames(x))
  x <- x[topGeneInstert,,drop=FALSE]
  #rownames(x) <- topGene48
  #x[is.na(x)] <- 0
  x <- t(scale(t(x)))
  #x[is.na(x)] <- 0
  # Normalization
  #x <- sweep(x, 1, m)
  #x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
  valide_df <- validate(x, rav_test)
  valide_val <- valide_df$score
  # PCA
  rm(x)
  
  val_list[[study]] <- valide_val
}
valAll <- data.frame(val_list)

scale_val <- t(scale(t(valAll)))

pcaVal = princomp(t(scale_val))
fviz_pca_biplot(pcaVal) +
  ggtitle("")


rownames(scale_val) <- rownames(valide_df)
val_Seurat <- CreateSeuratObject(scale_val)
valReduction <- CreateDimReducObject(embeddings = as.matrix(pcaVal$scores), key = "Val_", assay = DefaultAssay(val_Seurat))
val_Seurat@reductions$pca <- valReduction
DimPlot(val_Seurat)
FeaturePlot(val_Seurat,"RAV1")



#== valide NMF-------------------------
valide(x,RAVmodel)