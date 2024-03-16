setwd("../1.13_bmsc_aging/")
source("../0_utils/seurat_utils")
file <- list.dirs("data",recursive = F)
count10x <- lapply(file,Read10X)
fileName <- file%>%
  gsub("data/","",.)
names(count10x) <- fileName

limbSeurat <- map2(count10x,fileName,function(x,y) CreateSeuratObject(x,min.cells = 3, min.features = 500,project = y))
lapply(limbSeurat, dim)


limbSeurat <- lapply(limbSeurat,qcFun)
limbMerge <- merge(limbSeurat[[1]],limbSeurat[2:14])

limbMerge <- runharmony(limbMerge)
FeaturePlot(limbMerge,c("Sp7","Runx2","Acan","Sox9","Aspn","Postn","Alpl","Smpd3","Cxcl12"))
FeaturePlot(limbMerge,"Prrx1",label=T)
limbSub <- subset(limbMerge,idents=c("0","15","9"))
new.id <- paste0("BmscAging_Young_",fileName)
names(new.id) <- fileName
limbSub <- renameFun(limbSub,new.id)
limbSub <- DietSeurat(limbSub)
saveRDS(limbSub,"../data/cleandata/BmscAging_Young.Rds")



