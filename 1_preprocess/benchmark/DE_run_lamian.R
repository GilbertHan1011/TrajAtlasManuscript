rm(list=ls())
null1 <- readH5AD("processed_data/11.27_gene_null/null_model_1.h5ad")

null1_sub <- null1[,null1@colData$random_sample %in%c(1:5,25:29)]

null_col <- null1_sub@colData%>%as.data.frame()
expr <- null1_sub@assays@data$Magic
rownames(expr) <- rownames(null1)
colnames(expr) <- colnames(null1_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

dptRnk <- rank(dptData)
Res <- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics=Res$statistics
saveRDS(Res,"processed_data/11.27_gene_null/11.28_lamian_null1.Rds")


null2 <- readH5AD("processed_data/11.27_gene_null/null_model_trend.h5ad")

null_sub <- null2[,null2@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$magic_data
rownames(expr) <- rownames(null2)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res2 <- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics=Res$statistics
statics2 <- Res2$statistics




null3 <- readH5AD("processed_data/11.27_gene_null/null_model_mean.h5ad")

null_sub <- null3[,null2@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null3)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res3 <- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics3 <- Res3$statistics
statics=Res$statistics
statics2 <- Res2$statistics



null4 <- readH5AD("processed_data/11.27_gene_null/null_model_peak.h5ad")

null_sub <- null4[,null2@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null4)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res4<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics4 <- Res4$statistics
statics=Res$statistics
statics2 <- Res2$statistics






null5 <- readH5AD("processed_data/11.27_gene_null/null_random_trend.h5ad")

null_sub <- null5[,null5@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null5)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res5<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics5 <- Res5$statistics



null6 <- readH5AD("processed_data/11.27_gene_null/null_random_trend_percent6.h5ad")

null_sub <- null6[,null6@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null6)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res6<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics6 <- Res6$statistics
debug(plotXDEHm)



null7 <- readH5AD("processed_data/11.27_gene_null/12.1_trend_finely.h5ad")

null_sub <- null7[,null7@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null7)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res7<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics6 <- Res6$statistics
debug(plotXDEHm)





## get differential dynamic genes statistics
stat <- Res2$statistics
stat <- stat[order(stat[, 1],-stat[, 3]),]
## identify XDE genes with FDR.overall < 0.05 cutoff
diffgene <-
  rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05, ])

## population fit
Res2$populationFit <-
  getPopulationFit(Res2, gene = diffgene, type = 'variable')
## clustering
Res2$covariateGroupDiff <-
  getCovariateGroupDiff(testobj = Res2, gene = diffgene)
Res2$cluster <-
  clusterGene(Res2, gene = diffgene, type = 'variable', k = 5)

lamianHm_res2 <- plotXDEHm(
  Res2,
  cellWidthTotal = 180,
  cellHeightTotal = 350,
  subsampleCell = FALSE,
  sep = ':.*'
)


#== null  trend without mean-------------------


null8 <- readH5AD("processed_data/11.27_gene_null/null_model_trend_without_mean_v8.h5ad")

null_sub <- null8[,null8@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null8)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res8<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics8 <- Res8$statistics


null9 <- readH5AD("processed_data/11.27_gene_null/null_model_mean_finely_v9.h5ad")

null_sub <- null9[,null9@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null9)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res9<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics9 <- Res9$statistics



null10 <- readH5AD("processed_data/11.27_gene_null/null_shuffle_model_mean_finely_v10.h5ad")

null_sub <- null10[,null10@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null10)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res10<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics10 <- Res10$statistics



null11 <- readH5AD("processed_data/11.27_gene_null/null_shuffle_model_trend_withoutmean_finely_v11.h5ad")

null_sub <- null11[,null11@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null11)
colnames(expr) <- colnames(null_sub)

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res11<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics11 <- Res11$statistics




null12 <- readH5AD("processed_data/11.27_gene_null/null_shuffle_model_peak_finely_v12.h5ad")

null_sub <- null12[,null12@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null12)
colnames(expr) <- colnames(null_sub)

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()
dptData=dptData*100

dptRnk <- rank(dptData)

rownames(expr)
Res12<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptRnk,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics12 <- Res12$statistics




## get differential dynamic genes statistics
stat <- Res12$statistics
stat <- stat[order(stat[, 1],-stat[, 3]),]
## identify XDE genes with FDR.overall < 0.05 cutoff
diffgene <-
  rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05, ])

## population fit
Res12$populationFit <-
  getPopulationFit(Res12, gene = diffgene, type = 'variable')
## clustering
Res12$covariateGroupDiff <-
  getCovariateGroupDiff(testobj = Res12, gene = diffgene)

debug(clusterGene)
Res12$cluster <-
  clusterGene(Res12, gene = diffgene, type = 'variable', k = 5)

lamianHm_res12 <- plotXDEHm(
  Res12,
  cellWidthTotal = 180,
  cellHeightTotal = 350,
  subsampleCell = FALSE,
  sep = ':.*'
)









#==9 again-----------
null9_2 <- readH5AD("processed_data/11.27_gene_null/null_model_mean_finely_v9.h5ad")

null_sub <- null9_2[,null9_2@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null9_2)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res9_2<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
statics9_2 <- Res9_2$statistics





null10 <- readH5AD("processed_data/11.27_gene_null/3.7_null_shuffle_model_mean_finely_v10.h5ad")

null_sub <- null10[,null10@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null10)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res10_2<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
end_time_lamian <- Sys.time()
res10Time=end_time_lamian-start_time_lamian


#== time-test---------------




null10 <- readH5AD("processed_data/11.27_gene_null/3.7_time_500gene_v11.h5ad")

null_sub <- null10[,null10@colData$random_sample %in%c(1:10,20:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null10)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res11<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
end_time_lamian <- Sys.time()
res11Time=end_time_lamian-start_time_lamian

#== time 2 200 gene


null10 <- readH5AD("processed_data/11.27_gene_null/3.7_time_2000gene_v12.h5ad")

null_sub <- null10[,null10@colData$random_sample %in%c(1:10,20:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null10)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res12<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)
end_time_lamian <- Sys.time()
res12Time=end_time_lamian-start_time_lamian




null <- readH5AD("processed_data/11.27_gene_null/3.8_gene_test.h5ad")

null_sub <- null[,null@colData$random_sample %in%c(1:10,20:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res13<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)

statics13<- Res13$statistics

end_time_lamian <- Sys.time()
res13Time=end_time_lamian-start_time_lamian


#== 14-----------------------------------
null <- readH5AD("processed_data/11.27_gene_null/3.8_gene_test_v14.h5ad")

null_sub <- null[,null@colData$random_sample %in%c(1:10,20:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res14<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)

statics14<- Res14$statistics

end_time_lamian <- Sys.time()
res14Time=end_time_lamian-start_time_lamian


#==15--------------------------------------------
null <- readH5AD("processed_data/11.27_gene_null/3.8_gene_test_v15.h5ad")

null_sub <- null[,null@colData$random_sample %in%c(1:5,25:29)]

null_col <- null_sub@colData%>%as.data.frame()
expr <- null_sub@assays@data$data
rownames(expr) <- rownames(null)
colnames(expr) <- colnames(null_sub)
# expr=expr[rowSums(expr)!=0,]
# leprMeta <- leprMeta[,colnames(expr)]
# leprMeta <- leprMeta[(1:nrow(leprMeta)) %% 4 == 0, ]
# expr <- expr[,rownames(leprMeta)]

cellanno <- data.frame("Cell"=rownames(null_col),
                       "Sample"=as.character(null_col$random_sample))
rownames(cellanno) <- cellanno$Cell

designPre <- table(null_col$random_sample,null_col$Group)%>%as.data.frame%>%filter(Freq!=0)

design=model.matrix(~designPre$Var2)
rownames(design) <- designPre$Var1
dptData <- null_col$time
names(dptData) <- rownames(null_col)

cellanno$Sample <- as.character(cellanno$Sample)
start_time_lamian <- Sys.time()

Res15<- lamian_test(
  expr = expr,
  cellanno = cellanno,
  pseudotime = dptData,
  design = design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 100,
  ## This is for permutation test only. 
  ## We suggest that users use default permuiter = 100.
  ## Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores.fit = 10,
  ncores = 10,
  verbose.output = TRUE
)

statics15<- Res15$statistics

end_time_lamian <- Sys.time()
res15Time=end_time_lamian-start_time_lamian


time1 <- as.numeric(res10Time)*60
time2 <- as.numeric(res11Time)*60
time3 <- as.numeric(res12Time)*60*60
timeDf <- data.frame(c(time1,time2,time3))
trajDiffTime <- read.csv("processed_data/24.3.6_benchmark/time_trajDiff.csv")
timeSum <- data.frame(timeDf,trajDiffTime$X)
colnames(timeSum) <- c("Lamian","trajDiff")
rownames(timeSum) <- c("1e6","1e7","1e8")

timeSum <- rownames_to_column(timeSum,"Size")
timeSum_long <- pivot_longer(timeSum,-Size,names_to = "methods", values_to = "time" )

ggplot(timeSum_long,aes(x=Size,y=time,fill=methods))+
  geom_bar(position="dodge", stat="identity")+theme_bw(base_size = 20)
ggsave("result/3.8_benchmark/benchmark_time.pdf",width = 7,height = 5)


test1_trajdiff <- read.csv("processed_data/11.27_gene_null/3.9_static/mean_shuffle_v10_res.csv")
test2_trajdiff <- read.csv("processed_data/11.27_gene_null/3.9_static/test_trajdiff_v8.csv")
test3_trajdiff <- read.csv("processed_data/11.27_gene_null/3.9_static/genes_v14_res.csv")

trajdiff1 <- sum(test1_trajdiff$overall_gene_p<0.01)/length(test1_trajdiff$overall_gene_p)
trajdiff2 <- sum(test2_trajdiff$overall_gene_p<0.01)/length(test2_trajdiff$overall_gene_p)
trajdiff3 <- sum(test3_trajdiff$overall_gene_p<0.01)/length(test3_trajdiff$overall_gene_p)

statics10 <- as.data.frame(statics10)
lamian1 <- sum(statics10$fdr.overall<0.01)/length(statics10$fdr.overall)
statics8 <- as.data.frame(statics8)
lamian2 <- sum(statics8$fdr.overall<0.01)/length(statics8$fdr.overall)
statics14 <- as.data.frame(statics14)
lamian3 <- sum(statics14$fdr.overall<0.01)/length(statics14$fdr.overall)


trajdiff1 <- sum(test1_trajdiff$overall_gene_p<0.01)/length(test1_trajdiff$overall_gene_p)
trajdiff2 <- sum(test2_trajdiff$overall_gene_p<0.01)/length(test2_trajdiff$overall_gene_p)
trajdiff3 <- sum(test3_trajdiff$overall_gene_p<0.01)/length(test3_trajdiff$overall_gene_p)

test <- data.frame(c(lamian1,lamian2,lamian3),c(trajdiff1,trajdiff2,trajdiff3))
rownames(test) <- c("Mean","Trend","Generality")
colnames(test) <- c("Lamian","TrajDiff")


geneSum <- rownames_to_column(test,"Attr")
geneSum_long <- pivot_longer(geneSum,-Attr,names_to = "methods", values_to = "time" )

geneSum_long$Attr <- factor(geneSum_long$Attr,c("Mean","Trend","Generality"))
ggplot(geneSum_long,aes(x=Attr,y=time,fill=methods))+
  geom_bar(position="dodge", stat="identity")+theme_bw(base_size = 20)


ggsave("result/3.8_benchmark/benchmark_attr.pdf",width = 7,height = 5)
staticsList <- lapply(ls(pattern = "statics*"),get)
names(staticsList) <- ls(pattern = "statics*")[1:15]
saveRDS(staticsList,"processed_data/11.27_gene_null/3.9_static/staticsList.Rds")
write.csv(timeSum,"processed_data/11.27_gene_null/3.9_static/time.csv")
write.csv(geneSum,"processed_data/11.27_gene_null/3.9_static/attr.csv")
