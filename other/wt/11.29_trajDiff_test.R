peudobulk <- zellkonverter::readH5AD("processed_data/11.27_gene_null/null_2_pseudobulk.h5ad")
bulkdata=peudobulk@assays@data$X
bulkmeta=peudobulk@colData%>%as.data.frame()
library(edgeR)

design <-model.matrix(~0+Group+time, data=bulkmeta)
colnames(bulkdata) <- rownames(design)
fit <-glmQLFit(bulkdata, design)

libSize=rep(1,dim(design)[1])
dge = DGEList(counts=bulkdata,lib.size = libSize)
dge = calcNormFactors(dge, method="none")
dge = estimateDisp(dge, design)
fit = glmQLFit(dge, design, robust=TRUE)
my.contrasts <-makeContrasts(Group=" GroupGroup1-GroupGroup2",levels=design)
qlf <-glmQLFTest(fit, contrast=my.contrasts)
tt=topTags(qlf, sort.by ="none", n=Inf)
tt <- as.data.frame(tt)
