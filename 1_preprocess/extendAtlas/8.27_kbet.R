# devtools::install_github('theislab/kBET')
# devtools::install_local("../../../soft/kBET-master.zip")
library(kBET)
library(Seurat)
library(Matrix)
library(ggplot2)
# environment-----------------

outdir <- "result/8.27_kbet/"
dir.create(outdir)
limbMerge <- readRDS("../important_processed_data/8.27_merge_78w.Rds")
for (i in unique(limbMerge$Project)){
  study=i
  target=limbMerge[,limbMerge$Project==study]
  data = GetAssayData(object = target, slot = "counts",assay = "RNA")
  data <- t(data)
  batch <- target$Sample
  subset_size <- 1000/length(batch) #subsample to 10% of the data
  
  if (length(batch)>2000){
    subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
    batch.estimate <- kBET(data[subset_id,], batch[subset_id], plot=TRUE, do.pca = TRUE, dim.pca = 2)
  }
  else{
    batch.estimate <- kBET(data, batch, plot=TRUE, do.pca = TRUE, dim.pca = 15)
  }
  saveRDS(batch.estimate,paste0(outdir,study,"_estimateResult.Rds"))
  plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                    each=length(batch.estimate$stats$kBET.observed)), 
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
    labs(x='Test', y='Rejection rate',title='kBET test results') +
    theme_bw() +  
    scale_y_continuous(limits=c(0,1))
  ggsave(paste0(outdir,study,"_reject_box.pdf"))
}


library(ggplot2)
kbetFile <- list.files(path = "result/8.27_kbet/",pattern = "*.Rds",full.names = T)
kbetName <- gsub("result/8.27_kbet//","",kbetFile)
kbetName <- gsub("_estimateResult.Rds","",kbetName)
kbetList <- lapply(kbetFile,readRDS)
names(kbetList) <- kbetName
kbet <- lapply(kbetList,function(x) x$summary$kBET.observed)
kbetDf <- do.call(cbind,kbet)
boxplot(kbetDf)
ggplot2::geom_boxplot(kbetDf)
library(reshape2)

# Specify id.vars: the variables to keep but not split apart on
kbetlong=melt(kbetDf)[,-1]
ggplot(kbetlong,aes(x=Var2,y=value))+geom_boxplot()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = "bold"))+
  xlab("Project")+ylab("rejection rate")
ggsave("result/8.27_qc/overall.pdf",width = 12,height = 6)  





# create batch metadata---------------------------------------------------
source("../function/seurat_utils.R")
library(dplyr)
coreData <- readxl::read_xlsx("data/8.18_sampleMeta.xlsx")
limbMerge <- createSampleMeta(limbMerge,coreData,"batch")
limbMerge <- createSampleMeta(limbMerge,coreData,"Project")
limbMerge <- createSampleMeta(limbMerge,coreData,"Origin")
# limbMerge$short_id <- as.character(limbMerge$short_id)
# limbMerge$short_id[limbMerge$Project=="CalvariaP4_Ayturk"] <- "CAy"
# limbMerge$short_id[limbMerge$Project=="CranioSoxc_Angelozzi"] <- "CAn"
# bm <- RenameCells(limbMerge[,limbMerge$Project=="BMSC-Specification_Kishor"],add.cell.id="BMSCSpecification_Kishor_")
# View(colnames(limbMerge)%>%as.data.frame())
# colName <- colnames(limbMerge[,limbMerge$Project!="BMSC-Specification_Kishor"])
# colName <- c(colName,colnames(bm))
# test <- RenameCells(limbMerge,new.names=colName)
# limbMerge <- test
# rm(test)
write.csv(limbMerge@meta.data,"../unimportant_processed_data/8.25_integrateAll/8.27_metadata.csv")
SaveH5Seurat(limbMerge,"../data/annodata/1.18_limbMerge.h5Seurat")
Convert("../data/annodata/1.18_limbMerge.h5Seurat",dest="h5ad")

saveRDS(limbMerge,"../data/annodata/1.18_osteo_merge.Rds")

