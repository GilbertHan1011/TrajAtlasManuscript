#== This script consider injury as individual datasets--------------------


library(GenomicSuperSignature)
library(factoextra)
library(ComplexHeatmap)
#rm(list=ls())

injuryFile <- list.files("process_data/traj_diff/cpm/injury/",full.names = T)
shortFile <- list.files("process_data/traj_diff/cpm/injury/",full.names = F)%>%gsub(".csv","",.)
filterFileInjury <- paste0("process_data/traj_diff/cpm/filter/mes_sep_",shortFile,".csv")


test=data.table::fread(injuryFile[10],header = T)
test=data.table::fread("process_data/traj_diff/cpm/mes/RibRegeneraton_Serowoky_Control.csv",header = T)


filterFile <- function(index,fileName,shortFileName,filterFileName){
  countFile <- data.table::fread(fileName[index],header = T)
  countFile[is.na(countFile)] <- 0
  geneName=countFile$V1
  countFile <- countFile%>%as.data.frame()
  rownames(countFile) <-geneName
  countFile <- countFile[,-1]
  countFile=countFile[!colSums(countFile)==0]
  countFile=countFile[!rowSums(countFile)==0,]
  stopifnot(dim(countFile)[2]>30)
  data.table::fwrite(countFile,filterFileName[index],row.names=T)
}
errorIndex=list()
#dderrorIndex
for (i in 1:length(shortFile)){
  tryCatch({
    filterFile(i,injuryFile,shortFile,filterFileInjury)
  },error = function(err) {
    cat("An error occurred in iteration", i, ":", conditionMessage(err), "\n")
    errorIndex <<- c(errorIndex,i)
    cat("Expanding list")
    # Skip to the next iteration
    return(invisible(NULL))
  })
}
errorIndexNum <- unlist(errorIndex)
filterFileNameSub=filterFileName[!(seq_along(filterFileName) %in% errorIndexNum)]
filterFileNameSub <- unique(c(filterFileNameSub,filterFileInjury[!(seq_along(filterFileInjury) %in% errorIndexNum)]))
# 
shortNameSub=filterFileNameSub%>%gsub("process_data/traj_diff/cpm/filter/","",.)%>%gsub(".csv","",.)
#filterFileNameSub=filterFileName[!(seq_along(filterFileName) %in% errorIndexNum)]
#shortNameSub=shortFileName[!(seq_along(filterFileName) %in% errorIndexNum)]
findHighlyGene <- function(index,fullFileName){
  countFile <- data.table::fread(fullFileName[index])
  y <- apply(countFile[,-1], 1, sd) 
  names(y) <- countFile$V1
  y <- y[order(y, decreasing = TRUE)]
  return(y)
}
yList=list()
for (i in 1:length(filterFileNameSub)){
  y=findHighlyGene(i,filterFileNameSub)
  yList[[filterFileNameSub[i]]]=y
}

topGeneList <- list()
cutoff <- 0.6  # Adjust this if you want different cutoff

for (i in seq_len(length(yList))){
  k <- length(yList[[i]])
  k <- round(k*cutoff)
  topGenes <- names(yList[[i]])[1:k]
  topGeneList[[i]] <- topGenes
}

topGene <- topGeneList%>%unlist%>%table%>%as.data.frame()
topGeneSelect <- topGene%>%filter(Freq>260)%>%dplyr::select(".")%>%unlist()%>%as.character()

saveRDS(topGene,"process_data/trajMap/topGeneList.Rds")
saveRDS(topGeneSelect,"process_data/trajMap/topGeneSelect.Rds")

varGene <- readRDS("../3.9_wt_integrate/processed_data/12.12_TRAV/topgene_var2000.Rds")
trainingData_NMF=list()

n <- 15 #== base on cross-validation
index <- 1
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
  x[is.na(x)] <- 0
  model_test <- RcppML::nmf(t(x), n, verbose = F, seed = 1234)
  w <- model_test$w
  d <- model_test$d
  h <- model_test$h
  colnames(h) <- rownames(x)
  filterH <- rowSums(h)!=0
  h <- h[filterH,]
  trainingData_NMF[[study]]$rotation <- t(h)
  rm(x)
}

allZ_list <- lapply(trainingData_NMF, function(x) x$rotation)
allZ_list <- lapply(allZ_list, function(x) {
  x <- as.data.frame(x)
  x <- x[varGene,]
  rownames(x) <- varGene
  x[is.na(x)] <- 0
  return(x)
})
listLen <- lapply(allZ_list,function(x) dim(x)[2])
listLen <- listLen%>%unlist()
allZ_list=allZ_list[listLen>8]
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

normalNMF <- readRDS("../important_processed_data/24.1.trainingdataNMF.Rds")
normalLoad <- readRDS("../important_processed_data/24.1.core_avgload.Rds")

conbineLoad <- cbind(normalLoad,avg_load)
core_or_extend <- c(rep("core",ncol(normalLoad)),rep("extend",ncol(avg_load)))

my_color2<-c("red","blue")
names(my_color2) <- c("core","extend")
top_anno <- HeatmapAnnotation(
  ident =core_or_extend,
  col = list(
    ident=my_color2
  )
)

hm <- Heatmap(conbineLoad,top_annotation = top_anno,show_row_names = F,show_column_names = F,column_km = 20)
hm <- draw(hm)

clusterlist = column_order(hm)


countData <- hm@ht_list$matrix_2@matrix
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = colnames(countData)[clusterlist[[i]]],
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)
rownames(clu_df) <- clu_df$GeneID



#== 30-----------------

cl <- trainingData_PCclusters_60$cluster
silh_res <- cluster::silhouette(cl, res.dist)
cl_silh_width <- summary(silh_res)$clus.avg.widths
trainingData_PCclusters_60$sw <- cl_silh_width  # add silhouette width to the result

RAVmodel <- PCAGenomicSignatures(assays = list(RAVindex = as.matrix(trainingData_PCclusters_60$avgLoading)),trainingData =DataFrame(all))
names(trainingData_PCclusters_60$studies) <- RAVmodel@colData@rownames
colnames(train_data) <-  RAVmodel@colData@rownames
annot_database <- "MSigDB C5"   # will be added into the RAVmodel metadata
train_data=trainingData_PCclusters_60$avgLoading
colnames(train_data) <- names(trainingData_PCclusters_60$studies)

#RAVmodel <- PCAGenomicSignatures(trainingData = as.data.frame(trainingData_PCclusters_60$avgLoading))
DataFrame(trainingData_PCclusters_60[c("cluster")])
metadata(RAVmodel) <- trainingData_PCclusters_60[c("cluster", "size", "k", "n")]
names(metadata(RAVmodel)$size) <- paste0("RAV", seq_len(ncol(RAVmodel)))
geneSets(RAVmodel) <- annot_database
studies(RAVmodel) <- trainingData_PCclusters_60$studies
silhouetteWidth(RAVmodel) <- trainingData_PCclusters_60$sw
updateNote(RAVmodel) <- note


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
  tryCatch({
    valide_df <- validate_local(x, RAVmodel)
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

write.csv(valAll,"process_data/trajMap/1.11_validation_df_2nd.csv")
saveRDS(RAVmodel,"process_data/trajMap/1.11_2nd_ravmodel.Rds")
scale_val <- t(scale(t(valAll)))


# cluster------------
d <- 4

allZ_list <- lapply(trainingData_PCA, function(x) x$rotation)
allZ <- Reduce(cbind, allZ_list)
all <- t(allZ) 
print(paste("Dimension of allZ is", dim(allZ)))
res.dist <- factoextra::get_dist(all, method = "spearman")

nbMax <- fviz_nbclust(as.matrix(res.dist), kmeans, method = "wss",k.max = 100)
#fviz_nbclust(as.matrix(res.dist), kmeans, method = "gap_stat",k.max = 200)
nbMax+ geom_vline(xintercept = 35, linetype = 2)
k <- round(nrow(all)/d, 0)
start <- Sys.time()
res.hcut <- factoextra::hcut(res.dist, k = k, hc_func = "hclust", 
                             hc_method = "ward.D", hc_metric = "spearman")
end <- Sys.time()
d <- 4
k <- round(nrow(all)/d, 0)
start <- Sys.time()
res.hcut_605 <- factoextra::hcut(res.dist, k = k, hc_func = "hclust", 
                                 hc_method = "ward.D", hc_metric = "spearman")
res.hcut_30 <- factoextra::hcut(res.dist, k = 30, hc_func = "hclust", 
                                hc_method = "ward.D", hc_metric = "spearman")
fviz_dend(res.hcut_30, rect = TRUE, k_colors = 30, cex = 0.7, main = "Hierarchical Clustering Dendrogram")


res.hcut_10 <- factoextra::hcut(res.dist, k = 10, hc_func = "hclust", 
                                hc_method = "ward.D", hc_metric = "spearman")
fviz_dend(res.hcut_10, rect = TRUE, k_colors = 10,show_labels = F, cex = 0.7, main = "Hierarchical Clustering Dendrogram")

# determine optimized cluster number
#fviz_nbclust(my_data, kmeans, method = "gap_stat")

fviz_nbclust(as.matrix(res.dist), kmeans, method = "gap_stat",k.max = 200)

trainingData_PCclusters_605 <- buildAvgLoading(allZ, k, cluster = res.hcut_605$cluster)

trainingData_PCclusters_30 <- buildAvgLoading(allZ, 30, cluster = res.hcut_30$cluster)

## Silhouette Width
cl <- trainingData_PCclusters_605$cluster
silh_res <- cluster::silhouette(cl, res.dist)
cl_silh_width <- summary(silh_res)$clus.avg.widths
trainingData_PCclusters_605$sw <- cl_silh_width  # add silhouette width to the result



pca_summary <- list()
for (i in seq_along(trainingData_PCA)) {
  pca_summary[[i]] <- trainingData_PCA[[i]]$variance
  names(pca_summary)[i] <- names(trainingData_PCA)[i]
}


annot_database <- "MSigDB C5"   # will be added into the RAVmodel metadata
RAVmodel <- PCAGenomicSignatures(assays = list(RAVindex = as.matrix(trainingData_PCclusters_605$avgLoading)))
metadata(RAVmodel) <- trainingData_PCclusters_605[c("cluster", "size", "k", "n")]
names(metadata(RAVmodel)$size) <- paste0("RAV", seq_len(ncol(RAVmodel)))
geneSets(RAVmodel) <- annot_database
studies(RAVmodel) <- trainingData_PCclusters_605$studies
silhouetteWidth(RAVmodel) <- trainingData_PCclusters_605$sw
#metadata(RAVmodel)$MeSH_freq <- MeSH_freq
trainingData(RAVmodel)$PCAsummary <- pca_summary[rownames(trainingData(RAVmodel))]
#mesh(RAVmodel) <- trainingData_MeSH[rownames(trainingData(RAVmodel))]
updateNote(RAVmodel) <- note

#== GSEA------------------
dir2 <- system.file("scripts", package = "GenomicSuperSignature")
searchPathways_func <- file.path(dir2, "searchPathways.R")
source(searchPathways_func)  # load the function
gsea_all <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_h/")
gsea_m5 <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_m5/")
gsea_m5_all <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_all_m5/")


saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_605.Rds")

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

