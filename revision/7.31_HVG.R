library(scran)
test <- readRDS("data/finalData/Ablation_Matsushita_con1.Rds")
test <- FindVariableFeatures(test)
varGen1 <- VariableFeatures(test)

testSce <- as.SingleCellExperiment(test)
dec.default <- modelGeneVar(testSce)
dec.noweight <- modelGeneVar(testSce, density.weights=FALSE)

fit.default <- metadata(dec.default)
plot(fit.default$mean, fit.default$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression") 
curve(fit.default$trend(x), col="dodgerblue", add=TRUE, lwd=2)
fit.noweight <- metadata(dec.noweight)
curve(fit.noweight$trend(x), col="red", add=TRUE, lwd=2)
legend("topleft", col=c("dodgerblue", "red"), legend=c("Default", "No weight"), lwd=2)

dec.cv2.test <- modelGeneCV2(testSce)
fit.cv2.pbmc <- metadata(dec.cv2.test)
plot(fit.cv2.pbmc$mean, fit.cv2.pbmc$cv2, log="xy")
curve(fit.cv2.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)
decVar <- dec.cv2.test[order(dec.cv2.test$ratio, decreasing=TRUE),]
FDRgene <- rownames(decVar)[decVar$FDR<0.05]
intersect(FDRgene,varGen1) %>% length()


library(clusterProfiler)
gobp <- read.gmt("revision/data/m5.all.v2023.2.Mm.symbols.gmt")
geneSetdf <-  gobp[gobp$term=="GOBP_OSTEOBLAST_DIFFERENTIATION",]
geneSet1 <- gobp[gobp$term=="GOBP_OSTEOBLAST_DIFFERENTIATION",]$gene
geneSet2 <- gobp[gobp$term=="GOBP_OSTEOBLAST_DEVELOPMENT",]$gene
geneSet3 <- gobp[gobp$term=="GOBP_REGULATION_OF_OSTEOBLAST_DIFFERENTIATION",]$gene
geneSet4 <- gobp[gobp$term=="GOBP_POSITIVE_REGULATION_OF_OSTEOBLAST_PROLIFERATION",]$gene

intersect(varGen1,geneSet1) %>% length
intersect(varGen1,geneSet2) %>% length
intersect(varGen1,geneSet3) %>% length
intersect(varGen1,geneSet4) %>% length


metaFeature <- test@assays$RNA@meta.features
metaFeature <- metaFeature %>% 
  arrange(desc(vst.variance.standardized))
overlap <- as.integer(metaFeature %>% rownames %in% geneSet1)
cumulative_sum <- cumsum(overlap)

# Plot
plot(cumulative_sum, type = 'o', col = 'blue', 
     xlab = 'Position in Ordered Array', ylab = 'Cumulative Overlap',
     main = 'Cumulative Overlap Plot')
plot_data <- data.frame(
  Position = seq_along(metaFeature %>% rownames),
  CumulativeOverlap = cumulative_sum
)

# Plot using ggplot2
ggplot(plot_data, aes(x = Position, y = CumulativeOverlap)) +
  geom_line(color = 'blue') +         # Line plot
  geom_point(color = 'blue') +        # Points at each position
  labs(
    x = 'Position in Ordered Array', 
    y = 'Cumulative Overlap', 
    title = 'Cumulative Overlap Plot'
  ) +
  theme_minimal()  +
  geom_vline(xintercept = 2000,color = "red")

plot_data <- plot_data %>% 
  mutate(slope = CumulativeOverlap/Position)

ggplot(plot_data, aes(x = Position, y = slope)) +
  geom_line(color = 'blue') +         # Line plot
  geom_point(color = 'blue') +        # Points at each position
  labs(
    x = 'Position in Ordered Array', 
    y = 'Cumulative Overlap', 
    title = 'Cumulative Overlap Plot'
  ) +
  theme_minimal()  +
  geom_vline(xintercept = 2000,color = "red")

dir.create("revision/results/hvg")
hvgFun <- function(seurat,geneSet,name, outDir, assay = "RNA"){
  seurat <- FindVariableFeatures(seurat)
  metaFeature <- seurat@assays[[assay]]@meta.features
  metaFeature <- metaFeature %>% 
    arrange(desc(vst.variance.standardized))
  overlap <- as.integer(metaFeature %>% rownames %in% geneSet)
  cumulative_sum <- cumsum(overlap)
  plot_data <- data.frame(
    Position = seq_along(metaFeature %>% rownames),
    CumulativeOverlap = cumulative_sum
  )
  
  plot_data <- plot_data %>% 
    mutate(slope = CumulativeOverlap/Position)
  
  # Plot using ggplot2
  ggplot(plot_data, aes(x = Position, y = CumulativeOverlap)) +
    geom_line(color = 'blue') +         # Line plot
    geom_point(color = 'blue') +        # Points at each position
    labs(
      x = 'Position in Ordered Array', 
      y = 'Cumulative Overlap', 
      title = name
    ) +
    geom_vline(xintercept = 2000,color = "red")+
    xlim(c(0,20000))+
    theme_minimal() 
  
  ggsave(paste0(outDir,"/cum_",name,".pdf"))
  ggsave(paste0(outDir,"/cum_",name,".png"))
  
  ggplot(plot_data, aes(x = Position, y = slope)) +
    geom_line(color = 'blue') +         # Line plot
    geom_point(color = 'blue') +        # Points at each position
    labs(
      x = 'Position in Ordered Array', 
      y = 'Cumulative Overlap', 
      title = name
    ) +
    geom_vline(xintercept = 2000,color = "red")+
    xlim(c(0,20000))+
    ylim(0,0.3)+
    theme_minimal() 
  ggsave(paste0(outDir,"/slope_",name,".pdf"))
  ggsave(paste0(outDir,"/slope_",name,".png"))
  
  write.csv(metaFeature,paste0(outDir,"/feature_",name,".csv"))
  write.csv(plot_data,paste0(outDir,"/cumtable_",name,".csv"))
  
}

wt_integrate <- readRDS("important_processed_data/5.4_wtintegrate_full_seurat.Rds")

projectName <- unique(wt_integrate$Project)
for (i in projectName){
  print(paste0("working on ",i))
  seuratSubset <- wt_integrate[,wt_integrate$Project==i]
  hvgFun(seuratSubset,geneSet = geneSet1,name = i,outDir = "revision/results/hvg/",assay = "originalexp")
}



#hvgFun(test,geneSet = geneSet1,name = "Ablation_Matsushita_con1",outDir = "revision/results/hvg/")



fileName <- list.files(path = "revision/results/hvg/",pattern = "cumtable*",full.names = T)
files <- lapply(fileName,read.csv,row.names = 1)
fileNameShort <- fileName %>% gsub("revision/results/hvg//cumtable_","",.) %>% gsub(".csv","",.)
df <- lapply(files,function(x) x[2000,2])
newdf <- df %>% unlist
newdf <- data.frame(name = fileNameShort,slope = newdf)
newdf$slope <- newdf$slope/length(geneSet1)
newdf$name <- reorder(newdf$name, newdf$slope)
ggplot(newdf,aes(x = name, y = slope))+
  geom_point()+geom_bar(stat = 'identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1,size = 10,face = "bold"))
ggsave("revision/results/hvg/overlap.pdf",width = 6,height = 4)
        

backGround <- rownames(seuratSubset)
backGroundDf <- data.frame(term = "Background",gene = backGround)
backGroundDf <- rbind(backGroundDf,geneSetdf)
enrichRes <- enricher(varGen1,TERM2GENE = backGroundDf,universe =backGround) %>% as.data.frame()

fileNameHVG <- list.files(path = "revision/results/hvg/",pattern = "feature_*",full.names = T)
geneMeta1 <- read.csv(fileNameHVG[[1]],row.names = 1)

enrichDfAll <- lapply(fileNameHVG,function(x) {
  genes <- gsub("revision/results/hvg//feature_","",x) %>% gsub(".csv","",.)
  geneMeta <- read.csv(x, row.names = 1)
  geneHVG <- rownames(geneMeta)[geneMeta$vst.variable]
  enrichRes <- enricher(geneHVG,TERM2GENE = backGroundDf,universe =backGround) %>% as.data.frame()
  rownames(enrichRes) <- genes
  return(enrichRes)
})
enrichDfs <- do.call(rbind,enrichDfAll)

rownames(enrichDfs) <- projectName
enrichDfsSub <- enrichDfs %>% rownames_to_column("project") %>% dplyr::select(c("project","p.adjust")) %>% 
  mutate(logp = -log10(p.adjust))
enrichDfsSub$project <- reorder(enrichDfsSub$project, enrichDfsSub$logp)
ggplot(enrichDfsSub,aes(x = project, y = logp))+
  geom_point()+geom_bar(stat = 'identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1,size = 10,face = "bold"))
ggsave("revision/results/hvg/overlap_enrich.pdf",width = 6,height = 4)



