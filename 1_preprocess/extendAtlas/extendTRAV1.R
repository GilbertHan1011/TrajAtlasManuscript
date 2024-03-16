library(GenomicSuperSignature)
library(factoextra)
library(ComplexHeatmap)
#rm(list=ls())
#gene <- read.csv("../3.9_wt_integrate/processed_data/12.12_TRAV/12.12_gene.csv")
lineageTable=read.csv("process_data/paga/paga_long_dataframe.csv",row.names = 1)
lineageTable$dir="mes"
lineageTable$dir[lineageTable$Lineage=="Lepr+ BMSC"]="lepr"
lineageTable$dir[lineageTable$Lineage=="Fibroblast"]="fibro"
lineageTable$dir[lineageTable$Lineage=="Chondrocyte"]="chondro"

fileName <- paste0("process_data/traj_diff/cpm/",lineageTable$dir,"/",lineageTable$Sample,".csv")
shortFileName <- paste0(lineageTable$dir,"_sep_",lineageTable$Sample)
filterFileName <- paste0("process_data/traj_diff/cpm/filter/",shortFileName,".csv")

test=data.table::fread(fileName[10],header = T)
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
for (i in 1:length(shortFileName)){
  tryCatch({
    filterFile(i,fileName,shortFileName,filterFileName)
  },error = function(err) {
    cat("An error occurred in iteration", i, ":", conditionMessage(err), "\n")
    errorIndex <<- c(errorIndex,i)
    cat("Expanding list")
    # Skip to the next iteration
    return(invisible(NULL))
  })
}
errorIndexNum <- unlist(errorIndex)
# 

filterFileNameSub=filterFileName[!(seq_along(filterFileName) %in% errorIndexNum)]
shortNameSub=shortFileName[!(seq_along(filterFileName) %in% errorIndexNum)]
findHighlyGene <- function(index,fullFileName){
  countFile <- data.table::fread(fullFileName[index])
  y <- apply(countFile[,-1], 1, sd) 
  names(y) <- countFile$V1
  y <- y[order(y, decreasing = TRUE)]
  return(y)
}
yList=list()
#filterName <- c(filterFileFull,filterFileLepr,filterFilechondro,filterFileFibro)
#shortName <- c(paste0("mesenchyme_",mesLineageName),paste0("lepr_",leprLineageName),
#               paste0("chondro_",chondroLineageName),paste0("fibro_",fibroLineageName))
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
# shortName <- c(paste0("mesenchyme_",mesLineageName),paste0("lepr_",leprLineageName),
#                paste0("chondro_",chondroLineageName),paste0("fibro_",fibroLineageName))
# shortName2 <- c(paste0(mesLineageName,"_sep_Mesenchyme"),paste0(leprLineageName,"_sep_Lepr_BMSC"),
#                 paste0(chondroLineageName,"_sep_Chondro"),paste0(fibroLineageName,"_sep_Fibroblast"))



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
  # Normalization
  #x <- sweep(x, 1, m)
  #x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
  
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

# breaksList = seq(-2, 2, by = 0.1)
# 
# pheatmap(avg_load,show_rownames = F,scale = "row",breaks =breaksList, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))
# 
# breaksList = seq(0, 0.001, by = 0.000001)
# pheatmap(avg_load,show_rownames = F,breaks =breaksList, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))

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

hm <- Heatmap(conbineLoad,top_annotation = top_anno,show_row_names = F,show_column_names = F,column_km = 12)
hm <- draw(hm)

clusterlist = column_order(hm)


countData <- hm@ht_list$matrix_8@matrix
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = colnames(countData)[clusterlist[[i]]],
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)
rownames(clu_df) <- clu_df$GeneID



hm@matrix

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
#metadata(RAVmodel)$MeSH_freq <- MeSH_freq
#trainingData(RAVmodel)$PCAsummary <- pca_summary[rownames(trainingData(RAVmodel))]
#mesh(RAVmodel) <- trainingData_MeSH[rownames(trainingData(RAVmodel))]
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
  #x[is.na(x)] <- 0
  # Normalization
  #x <- sweep(x, 1, m)
  #x <- sweep(x, 1, s, "/")
  # x <- preprocessCore::normalize.quantiles(x)   # quantile normalization
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

write.csv(valAll,"process_data/trajMap/1.10_validation_df.csv")
saveRDS(valAll,"process_data/trajMap/1.10_orig_ravmodel.Rds")
scale_val <- t(scale(t(valAll)))


pcaVal = princomp(scale_val)
fviz_pca_biplot(pcaVal) +
  ggtitle("")












# 
# 
# 
# trainingData_PCA=list()
# #== Run PCA-------------------------------
# n <- 20
# index <- 1
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
#   
#   # PCA
#   pca_res <- prcomp(t(x))   # x is a matrix with genes(row) x samples(column)
#   trainingData_PCA[[study]]$rotation <- pca_res$rotation[,1:20]
#   colnames(trainingData_PCA[[study]]$rotation) <- paste0(study, ".PC", c(1:n))
#   eigs <- pca_res$sdev^2
#   pca_summary <- rbind(SD = sqrt(eigs),
#                        Variance = eigs/sum(eigs),
#                        Cumulative = cumsum(eigs)/sum(eigs))
#   trainingData_PCA[[study]]$variance <- pca_summary[,1:n]
#   colnames(trainingData_PCA[[study]]$variance) <- paste0(study, ".PC", c(1:n))
#   
#   rm(x)
# }
# str(trainingData_PCA$CalvariaP4_Ayturk)
# saveRDS(trainingData_PCA,"processed_data/12.12_TRAV/whole_PCA.Rds")



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


#== final model--------------------
# 
# dir <- system.file("extdata", package = "GenomicSuperSignaturePaper")
# load(file.path(dir, "MeSH_terms_1399refinebio.rda"))
# x <- mesh_table
# # Update bagOfWords and MeSH_freq
# MeSH_freq <- table(x$name) %>% sort(., decreasing = TRUE)  # freq. of each term
# for (i in 1:nrow(x)) {x$bagOfWords[i] <- MeSH_freq[x$name[i]]}  # add freq. to the main table
# 
# # Update bagOfWords and MeSH_freq
# MeSH_freq <- table(x$name) %>% sort(., decreasing = TRUE)  # freq. of each term
# for (i in 1:nrow(x)) {x$bagOfWords[i] <- MeSH_freq[x$name[i]]}  # add freq. to the main table
# unique_id <- unique(x$identifier)
# 
# all_MeSH <- vector("list", length(unique_id))
# names(all_MeSH) <- unique_id
# for (study in unique_id) {
#   ind <- grepl(study, x$identifier)
#   all_MeSH[[study]] <- x[ind, c("score", "identifier", "concept", "name", 
#                                 "explanation", "bagOfWords")]
# }
# 
# # Subset only to the studies used in the model 
# trainingData_MeSH <- all_MeSH[allStudies]


#source("../function/12.12_TRAV_build_GSEA.R")

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

#== GSEA------------------
dir2 <- system.file("scripts", package = "GenomicSuperSignature")
searchPathways_func <- file.path(dir2, "searchPathways.R")
source(searchPathways_func)  # load the function
# gsea_all <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_h/")
# gsea_m5 <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_m5/")
gsea_m5_all_30 <- searchPathways(RAVmodel, "processed_data/12.12_TRAV/geneset/gsea_all_m5_30/")
lapply(gsea_m5_all_30[25:30], function(x) x["Description"])
saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_30.Rds")
RAVmodel
saveRDS(RAVmodel,"processed_data/12.12_TRAV/Ravmodel_30_scale.Rds")
