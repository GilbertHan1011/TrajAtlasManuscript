library(scPOP)
source("function/seurat_utils.R")
seuratSubset <- RunPCA()
test <- runSeurat(test)
testSce <- as.SingleCellExperiment(test)

metrics <-run_all_metrics(reduction = reducedDim(testSce, 'PCA'), 
                          metadata = colData(testSce),
                          batch_key = 'seurat_clusters',
                          label1_key = 'ident',
                          label2_key = 'label', 
                          run_name = 'example')
ari(test$seurat_clusters,test$predicted.id)

seuratSubset <- wt_integrate[,wt_integrate$Sample==sampleName[[1]]]


seuratSubset <-  runSeurat(seuratSubset)

ari(seuratSubset$seurat_clusters,seuratSubset$C19_named)
nmi(seuratSubset$seurat_clusters,seuratSubset$C19_named)


hvgNum <- c(50,200,500,1000,2000,4000,8000)


hvgSelect <- function(x,hvgNum,dim=30){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvgNum)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- FindNeighbors(x, dims = 1:dim)
  x <- FindClusters(x, resolution = 0.5)
  score_ari <- ari(Idents(x),x$C19_named)
  score_nmi <- nmi(Idents(x),x$C19_named)
  return(c(score_ari,score_nmi))
}

metric1 <- hvgSelect(seuratSubset,50)

metricList <- lapply(hvgNum,hvgSelect,x = seuratSubset,dim = 30)
metricDf <- do.call(rbind,metricList)
colnames(metricDf) <- c("ARI","NMI")
rownames(metricDf) <- hvgNum

for (i in sampleName){
  print(paste0("working on ",i))
  seuratSubset <- wt_integrate[,wt_integrate$Sample==i]
  metricList <- lapply(hvgNum,hvgSelect,x = seuratSubset,dim = 30)
  metricDf <- do.call(rbind,metricList)
  colnames(metricDf) <- c("ARI","NMI")
  rownames(metricDf) <- hvgNum
  write.csv(metricDf,paste0("revision/results/hvg/metric_",i,".csv"))
}



sampleName <- unique(wt_integrate$Sample)
for (i in sampleName){
  print(paste0("working on ",i))
  seuratSubset <- wt_integrate[,wt_integrate$Sample==i]
  hvgFun(seuratSubset,geneSet = geneSet1,name = i,outDir = "revision/results/hvg/",assay = "originalexp")
}



metric1 <- read.csv("revision/results/hvg/metric_BmscSpecification_Kishor_1.csv",row.names = 1)
metric1 <- metric1 %>% rownames_to_column("hvg_number")
metric1_long <- metric1 %>%
  pivot_longer(cols = c(ARI, NMI), names_to = "Metric", values_to = "Value") %>% 
  mutate(hvg_number = as.numeric(hvg_number))
# ggplot(metric1_long, aes(x = factor(hvg_number), y = Value, fill = Metric)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(
#     x = "hvg_number",
#     y = "Value",
#     title = "ARI and NMI Metrics by Size",
#     fill = "Metric"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 10, face = "bold"),
#     legend.position = "top"
#   )

ari_data <- metric1_long %>% filter(Metric == "ARI")

# Plot ARI
plot_ari <- ggplot(ari_data, aes(x = factor(hvg_number), y = Value)) +
  geom_bar(stat = "identity") +
  labs(
    x = "hvg_number",
    y = "Value",
    title = "ARI Metric"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 10, face = "bold"),
    legend.position = "none" # No need for legend in individual plots
  )
nmi_data <- metric1_long %>% filter(Metric == "NMI")
plot_nmi <- ggplot(nmi_data, aes(x = factor(hvg_number), y = Value)) +
  geom_bar(stat = "identity") +
  labs(
    x = "hvg_number",
    y = "Value",
    title = "NMI Metric"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 10, face = "bold"),
    legend.position = "none" # No need for legend in individual plots
  )
# Combine the two plots
combined <- plot_ari | plot_nmi 
combined+
  plot_annotation(
    title = "BmscSpecification_Kishor_1",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave("revision/results/hvg/combine_metric_BmscSpecification_Kishor_1.pdf")

plotMetric <- function(metric,name){
  metric <- metric %>% rownames_to_column("hvg_number")
  metric1_long <- metric %>%
    pivot_longer(cols = c(ARI, NMI), names_to = "Metric", values_to = "Value") %>% 
    mutate(hvg_number = as.numeric(hvg_number))

  ari_data <- metric1_long %>% filter(Metric == "ARI")
  
  # Plot ARI
  plot_ari <- ggplot(ari_data, aes(x = factor(hvg_number), y = Value)) +
    geom_bar(stat = "identity") +
    labs(
      x = "hvg_number",
      y = "Value",
      title = "ARI Metric"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 10, face = "bold"),
      legend.position = "none" # No need for legend in individual plots
    )
  nmi_data <- metric1_long %>% filter(Metric == "NMI")
  plot_nmi <- ggplot(nmi_data, aes(x = factor(hvg_number), y = Value)) +
    geom_bar(stat = "identity") +
    labs(
      x = "hvg_number",
      y = "Value",
      title = "NMI Metric"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 10, face = "bold"),
      legend.position = "none" # No need for legend in individual plots
    )
  # Combine the two plots
  combined <- plot_ari | plot_nmi 
  combined+
    plot_annotation(
      title = name,
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
}

metric2 <- read.csv("revision/results/hvg/metric_CranioSoxc_Angelozzi_WTE15.5.csv",row.names = 1)
plotMetric(metric2,name = "CranioSoxc_Angelozzi_WTE15.5")
ggsave("revision/results/hvg/combine_metric_CranioSoxc_Angelozzi_WTE15.5.pdf")


metricFile <- list.files(path = "revision/results/hvg/",pattern = "metric_*",full.names = T)
metricFile <- grep(".csv",metricFile,value = T)
files <- lapply(metricFile,read.csv,row.name = 1)
ariMetric <- lapply(files,function(x) hvgNum[which.max(x[,1])])
ariMetric <- unlist(ariMetric)
#hvgNum[which.max(files[[1]][,1])]
nmiMetric <- lapply(files,function(x) hvgNum[which.max(x[,2])])
nmiMetric <- unlist(nmiMetric)

test@assays$RNA@meta.features
