setwd("~/Desktop/disk1/limb/")
library(Seurat)
peristeal <- readRDS("data/1.6_peristeal_injury.Rds")
FeaturePlot(peristeal,c("Sp7","Bglap","Postn"))
DimPlot(peristeal,group.by = "orig.ident")
peristealSce <- as.SingleCellExperiment(peristeal)
peristealSce@assays@data$counts
zellkonverter::writeH5AD(peristealSce,"data/24.7.28_peristeal.h5ad")
# 

peristeal2 <- readRDS("data/finalData/Periosteal2018_Shawon.Rds")
FeaturePlot(peristeal2,c("Acan","Col2a1","Sp7","Alpl","Postn","S100a4","Ctsk","Ly6a"),ncol = 4)
ggsave("revision/results/peristeal/7.29_featureplot.png",width = 10,height = 4)
