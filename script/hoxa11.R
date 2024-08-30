setwd("~/Desktop/disk1/limb/")
wt_integrate <- readRDS("important_processed_data/5.4_wtintegrate_full_seurat.Rds")
FeaturePlot(wt_integrate,"Hoxa11",split.by = "Organ",raster = T)
dir.create("revision/results/hox")
ggsave("revision/results/hox/8.13_hox11_exp.pdf",width = 10,height = 3)

FeaturePlot(wt_integrate,"Foxp2",split.by = "Organ",raster = T)
ggsave("revision/results/hox/8.13_foxp2_exp.pdf",width = 10,height = 3)

FeaturePlot(wt_integrate,"Tbx15",split.by = "Organ",raster = T)
FeaturePlot(wt_integrate,"Msx2",split.by = "Organ",raster = T)
ggsave("revision/results/hox/8.13_msx2_exp.pdf",width = 10,height = 3)
FeaturePlot(wt_integrate,"Fos",raster = T)
