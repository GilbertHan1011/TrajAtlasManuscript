library(GenomicSuperSignature)
library(factoextra)
rm(list=ls())
lineageName=read.csv("processed_data/12.9_trajMap/12.9_coorDf.csv",row.names = 1)%>%colnames()

split_names <- strsplit(lineageName, "_sep_")%>%as.data.frame()%>%t()%>%as.data.frame()


mesFile <- list.files("processed_data/12.12_TRAV/cpmList/mesenchyme/",full.names = TRUE)
mesName <- list.files("processed_data/12.12_TRAV/cpmList/mesenchyme/")%>%gsub(".csv","",.)
mesLineageName <- split_names$V1[split_names$V2=="Mesenchyme"]
leprLineageName <- split_names$V1[split_names$V2=="Lepr_BMSC"]
fibroLineageName <- split_names$V1[split_names$V2=="Fibroblast"]
chondroLineageName <- split_names$V1[split_names$V2=="Chondro"]
mesFile <- mesFile[mesName%in%mesLineageName]
leprFile <- paste0("processed_data/12.12_TRAV/cpmList/lepr/",leprLineageName,".csv")
fibroFile <- paste0("processed_data/12.12_TRAV/cpmList/fibro/",fibroLineageName,".csv")
chondroFile <- paste0("processed_data/12.12_TRAV/cpmList/chondro/",chondroLineageName,".csv")

mes_test <- read.csv("processed_data/12.12_TRAV/cpmList/mesenchyme/CalvariaP4_Ayturk.csv")



gene <- read.csv("processed_data/12.12_TRAV/12.12_gene.csv")
rownames(mes_test) <- rownames(gene)
mes_test[is.na(mes_test)] <- 0
mes_test=mes_test[!colSums(mes_test)==0]
mes_test=mes_test[!rowSums(mes_test)==0,]
y <- apply(mes_test, 1, sd) 
y <- y[order(y, decreasing = TRUE)]
countFile1 <- read.csv(mesFile[1])

countFile2 <- data.table::fread(mesFile[1],header = T)
rownames(countFile2) <- rownames(gene)


countFile2[is.na(countFile2)] <- 0
countFile2 <- countFile2%>%as.data.frame()
countFile2=countFile2[!colSums(countFile2)==0]
countFile2=countFile2[!rowSums(countFile2)==0,]


countFile2[is.na(countFile2)] <- 0
countFile2 <- countFile2%>%as.data.frame()
countFile2[!colSums(countFile2) == 0]
filterFile <- function(index,fullFileName,shortFileName,lineageName){
  countFile <- data.table::fread(fullFileName[index],header = T)
  countFile[is.na(countFile)] <- 0
  countFile <- countFile%>%as.data.frame()
  rownames(countFile) <- rownames(gene)
  countFile <- countFile[,-1]
  countFile=countFile[!colSums(countFile)==0]
  countFile=countFile[!rowSums(countFile)==0,]
  data.table::fwrite(countFile,paste0("processed_data/12.12_TRAV/filter/",lineageName,"/",shortFileName[index],".csv"),row.names=T)
}

for (i in 1:length(mesFile)){
  filterFile(i,mesFile,mesLineageName)
}
for (i in 1:length(leprFile)){
  filterFile(i,leprFile,leprLineageName,"lepr")
}
for (i in 1:length(chondroFile)){
  filterFile(i,chondroFile,chondroLineageName,"chondro")
}
for (i in 1:length(fibroFile)){
  filterFile(i,fibroFile,fibroLineageName,"fibro")
}

filterFileFull <- paste0("processed_data/12.12_TRAV/filter/mesenchyme/",mesLineageName,".csv")
filterFileLepr <- paste0("processed_data/12.12_TRAV/filter/lepr/",leprLineageName,".csv")
filterFilechondro <- paste0("processed_data/12.12_TRAV/filter/chondro/",chondroLineageName,".csv")
filterFileFibro <- paste0("processed_data/12.12_TRAV/filter/fibro/",fibroLineageName,".csv")
test <- data.table::fread(filterFileFull[1])

findHighlyGene <- function(index,fullFileName){
  countFile <- data.table::fread(fullFileName[index])
  y <- apply(countFile[,-1], 1, sd) 
  names(y) <- countFile$V1
  y <- y[order(y, decreasing = TRUE)]
  return(y)
}
yList=list()
filterName <- c(filterFileFull,filterFileLepr,filterFilechondro,filterFileFibro)
shortName <- c(paste0("mesenchyme_",mesLineageName),paste0("lepr_",leprLineageName),
               paste0("chondro_",chondroLineageName),paste0("fibro_",fibroLineageName))
for (i in 1:length(filterName)){
  y=findHighlyGene(i,filterName)
  yList[[filterName[i]]]=y
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
topGene48 <- topGene%>%filter(Freq>48)%>%dplyr::select(".")%>%unlist()%>%as.character()
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
str(trainingData_PCA$CalvariaP4_Ayturk)
saveRDS(trainingData_PCA,"processed_data/12.12_TRAV/whole_PCA.Rds")



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


source("../function/12.12_TRAV_build_GSEA.R")

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
