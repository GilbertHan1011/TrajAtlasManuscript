library(Seurat)
library(scPOP)
library(ggpubr)
test <- readRDS("data/finalData/Suture2021_Farmer_E15.Rds")
test
test <- dbFinderFun(test)

test <- runSeurat(test)
label1 <- test$seurat_clusters
test <- runSeurat(test,seeds = 1223)
label3 <- test$seurat_clusters

testSub <- test[,test$scDblFinder_class=="singlet"]
testSub <- runSeurat(testSub)
label2 <- testSub$seurat_clusters

scPOP::ari(label1,label3)

ari(label1[colnames(testSub)],label3[colnames(testSub)])
ari(label1[colnames(testSub)],label2[colnames(testSub)])

testSub$label1 <- label1[colnames(testSub)]
testSub$label3 <- label3[colnames(testSub)]
p1 <- DimPlot(testSub,group.by = "label1")
p3 <- DimPlot(testSub,group.by = "label3")
p1|p3
DimPlot(test,group.by = "scDblFinder_class")


testSub <- FindClusters(testSub, resolution = 0.5,random.seed = 12341)
labelTest <- testSub$seurat_clusters
testSub <- FindClusters(testSub, resolution = 0.5,random.seed = 12344541)
labelTest2 <- testSub$seurat_clusters


ari(labelTest,labelTest2)


sampleCell <- sample(colnames(testSub),size = 3000)
diffCell <- setdiff(colnames(testSub),sampleCell)
dblCell <- colnames(test)[test$scDblFinder_class == "doublet"]
diffSubCell <-  sample(diffCell,size = length(dblCell))
test1Cell <- c(sampleCell,diffSubCell)
test2Cell <- c(sampleCell,dblCell)
test3Cell <- sampleCell

subset1 <- test[,test1Cell]
subset1 <- runSeurat(subset1)
labelSub1 <- subset1$seurat_clusters
subset2 <- test[,test2Cell]
subset2 <- runSeurat(subset2)
labelSub2 <- subset2$seurat_clusters
subset3 <- test[,test3Cell]
subset3 <- runSeurat(subset3)
labelSub3 <- subset3$seurat_clusters

ari(labelSub1[test3Cell],labelSub3[test3Cell])
ari(labelSub2[test3Cell],labelSub3[test3Cell])

nmi(labelSub1[test3Cell],labelSub3[test3Cell])
nmi(labelSub2[test3Cell],labelSub3[test3Cell])




runDbl <- function(seurat,name){
  seurat <- dbFinderFun(seurat)
  seuratSub <- seurat[,seurat$scDblFinder_class=="singlet"]
  dblLength <- sum(seurat$scDblFinder_class!="singlet")
  seuratLength <- dim(seuratSub)[2]
  sampleSize = seuratLength-dblLength
  
  sampleCell <- sample(colnames(seuratSub),size = sampleSize)
  diffCell <- setdiff(colnames(seuratSub),sampleCell)
  dblCell <- colnames(seurat)[seurat$scDblFinder_class == "doublet"]
  dblCombineCell <- c(sampleCell,dblCell)
  
  subset1 <- seurat[,dblCombineCell]
  subset1 <- runSeurat(subset1)
  labelSub1 <- subset1$seurat_clusters
  subset2 <- seuratSub
  subset2 <- runSeurat(subset2)
  labelSub2 <- subset2$seurat_clusters
  subset3 <- seurat[,sampleCell]
  subset3 <- runSeurat(subset3)
  labelSub3 <- subset3$seurat_clusters
  
  ari1 <- ari(labelSub1[sampleCell],labelSub3[sampleCell])
  ari2 <- ari(labelSub2[sampleCell],labelSub3[sampleCell])
  
  nmi1 <- nmi(labelSub1[sampleCell],labelSub3[sampleCell])
  nmi2 <- nmi(labelSub2[sampleCell],labelSub3[sampleCell])
  metric <- data.frame(c(ari1,ari2),c(nmi1,nmi2))
  DimPlot(seurat,group.by = "scDblFinder_class")
  ggsave(paste0("revision/results/doublet/",name,".pdf"))
  return(metric)
}

metric <- runDbl(test,"test")



dataList <- list.files("data/finalData/",full.names = T)

nameList <- dataList %>% gsub("data/finalData//","",.) %>% gsub(".Rds","",.) 

dblPipeline <- function(path){
  name <- path %>% gsub("data/finalData//","",.) %>% gsub(".Rds","",.) 
  seurat <- readRDS(path)
  metrics <- runDbl(seurat,name = name)
}
metricList <- list()
for (i in dataList){
  metricList[[i]] <- dblPipeline(i)
}
metricList

seurat <- readRDS(i)
seurat <- dbFinderFun(seurat)

testSub <- seurat[,seurat$scDblFinder_class=="singlet"]

metricList_bk <- metricList
metricList <- list()
for (i in dataList) {
  tryCatch(
    {
      metricList[[i]] <- dblPipeline(i)
    },
    error = function(e) {
      # Handle the error (e) here
      message("Error occurred for item: ", i)
      message("Error message: ", e$message)
      # Optionally, you could assign NULL or some default value to metricList[[i]] in case of error
      metricList[[i]] <- NULL
    }
  )
}
metricList[[1]][,1]

metricAri <- lapply(metricList,function(x) x[,1])
metricAriDf <- do.call(rbind,metricAri) %>% as.data.frame()
metricAriDf$ratio <-  metricAriDf$V1/(metricAriDf$V2)
  
hist(metricAriDf$ratio,breaks = 10)

result <- t.test(metricAriDf$V1, metricAriDf$V2, paired = TRUE)
  
metricAriDf <- metricAriDf %>% mutate(ID = (1:dim(metricAriDf)[1]))
ggplot(metricAriDf, aes(x = ID)) +
  geom_line(aes(y = V1, color = "Before"), size = 1) +
  geom_line(aes(y = V2, color = "After"), size = 1) +
  geom_point(aes(y = V1, color = "Before"), size = 3) +
  geom_point(aes(y = V2, color = "After"), size = 3) +
  labs(x = "Pair ID", y = "Value", title = "Line Plot of Paired Data") +
  scale_color_manual(values = c("Before" = "blue", "After" = "red")) +
  theme_minimal()

ggpaired(metricAriDf, 
         cond1      = "V1", 
         cond2      = "V2",
         color      = "condition", 
         line.color = "gray", 
         line.size  = 0.4,
         palette    = "jco")+ 
  annotate("text", x = 1.5, y = 1.03, 
           label = "p-value = 0.76", 
           size = 5, color = "black", fontface = "italic")+
  ggtitle("ARI metric")
ggsave("revision/results/doublet/ari_metric_paired.pdf")
metricNmi <- lapply(metricList,function(x) x[,2])
metricNmiDf <- do.call(rbind,metricNmi) %>% as.data.frame()
t.test(metricNmiDf$V1, metricNmiDf$V2, paired = TRUE)
ggpaired(metricNmiDf, 
         cond1      = "V1", 
         cond2      = "V2",
         color      = "condition", 
         line.color = "gray", 
         line.size  = 0.4,
         palette    = "jco")+ 
  annotate("text", x = 1.5, y = 1.03, 
           label = "p-value = 0.83", 
           size = 5, color = "black", fontface = "italic")+
  ggtitle("NMI metric")

ggsave("revision/results/doublet/NMI_metric_paired.pdf")

write.csv(metricAriDf,"revision/process/doublet/ariDF.csv")
write.csv(metricNmiDf,"revision/process/doublet/nmiDF.csv")
