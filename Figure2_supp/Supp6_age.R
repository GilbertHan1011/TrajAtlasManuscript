library('variancePartition')
library('edgeR')
library('BiocParallel')
library('plyr')
library(UpSetR)
bulk <- read.csv("../3.9_wt_integrate/result/9.30_pesudobulk/matrix_C7_new.csv",header = T,row.names = 1)
countMatrix <- t(as.matrix(bulk))
meta <- read.csv("../3.9_wt_integrate/result/9.30_pesudobulk/meta_c7.csv",row.names = 1)
meta$Organ <- as.factor(meta$Organ) # ensure it is a factor
meta$Organ<- relevel(meta$Organ, ref="Head")
meta$Age <-  factor(meta$Age,levels = c("Organogenesis stage", "Fetal stage", "Postnatal", 
                                                 "Young Adult","Adult",  "Old"))
meta$Age_num <- as.numeric(meta$Age)
meta$Age[meta$Age=="Young Adult"] <- "Yong_adult"
metaSub=meta[meta$C7_named=="Lepr+ BMSC"&meta$Organ=="Limb_adult",]
matrixSub=countMatrix[,meta$C7_named=="Lepr+ BMSC"&meta$Organ=="Limb_adult"]

Group <-factor(paste(metaSub$Age,metaSub$Age,sep="."))
geneExpr = DGEList(matrixSub,group = Group)
keep <-filterByExpr(geneExpr)
geneExpr<-geneExpr[keep,,keep.lib.sizes=FALSE]
geneExpr = calcNormFactors(geneExpr, method="TMM")

design <-model.matrix(~0+Group)
my.contrasts <-makeContrasts(
  yong_post=GroupYong_adult.Yong_adult-GroupPostnatal.Postnatal,
  adult_yong=GroupAdult.Adult-GroupYong_adult.Yong_adult,
  old_adult=GroupOld.Old-GroupYong_adult.Yong_adult,
  levels = design
)

geneExpr = estimateDisp(geneExpr, design)
fit = glmQLFit(geneExpr, design,robust = T)
deg_1 <- topTags(glmQLFTest(fit, contrast=my.contrasts[,1]), sort.by = "none", n=Inf)%>%as.data.frame()
deg_2 <- topTags(glmQLFTest(fit, contrast=my.contrasts[,2]), sort.by = "none", n=Inf)%>%as.data.frame()
deg_3 <- topTags(glmQLFTest(fit, contrast=my.contrasts[,3]), sort.by = "none", n=Inf)%>%as.data.frame()




# 
matrixSub=countMatrix[,meta$C7_named=="Lepr+ BMSC"&meta$Organ=="Limb_adult"]
metaSub=meta[meta$C7_named=="Lepr+ BMSC"&meta$Organ=="Limb_adult",]
rm(Group)
Group <- metaSub$Age_num
geneExpr = DGEList(matrixSub,group = Group)
keep <-filterByExpr(geneExpr)
geneExpr<-geneExpr[keep,,keep.lib.sizes=FALSE]
vobjDream = voomWithDreamWeights(geneExpr, form, meta)
vobjDreamList[[i]]=vobjDream
# geneExpr = calcNormFactors(geneExpr, method="TMM")
# 
# design <-model.matrix(~Group)

# 
# geneExpr = estimateDisp(geneExpr, design)
# fit = glmQLFit(geneExpr, design,robust = T)

sigAge <- topTags(glmQLFTest(fit), sort.by = "none", n=Inf)%>%as.data.frame()
upGene1 <- deg_1%>%
  dplyr::filter(FDR<0.05&logFC>0.5)%>%
  rownames()
upGene2 <- deg_2%>%
  dplyr::filter(FDR<0.05&logFC>0.5)%>%
  rownames()
upGene3 <- deg_3%>%
  dplyr::filter(FDR<0.05&logFC>0.5)%>%
  rownames()

upGene1 <- deg_1%>%
  dplyr::filter(FDR<0.05&logFC>0)%>%
  rownames()
upGene2 <- deg_2%>%
  dplyr::filter(FDR<0.05&logFC>0)%>%
  rownames()
upGene3 <- deg_3%>%
  dplyr::filter(FDR<0.05&logFC>0)%>%
  rownames()

UpSetR::upset(fromList(list(list1=upGene1,list2=upGene2,list3=upGene3)))


downGene1 <- deg_1%>%
  dplyr::filter(FDR<0.05&logFC< -0.5)%>%
  rownames()
downGene2 <- deg_2%>%
  dplyr::filter(FDR<0.05&logFC< -0.5)%>%
  rownames()
downGene3 <- deg_3%>%
  dplyr::filter(FDR<0.05&logFC< -0.2)%>%
  rownames()
UpSetR::upset(fromList(list(list1=downGene1,list2=downGene2,list3=downGene3)))
geneTab <- c(upGene1,upGene2,upGene3)%>%table%>%as.data.frame()
up_intersect <- geneTab%>%filter(Freq>1)%>%dplyr::select(".")%>%unlist()%>%as.character()

write.csv(geneTab,"result/12.16_age_lepr_stem/up_gene_list.csv")
geneTab <- c(downGene1,downGene2,downGene3)%>%table%>%as.data.frame()
down_intersect <- geneTab%>%filter(Freq>1)%>%dplyr::select(".")%>%unlist()%>%as.character()
write.csv(geneTab,"result/12.16_age_lepr_stem/down_gene_list.csv")


#== heatmap---------------------------


diffGene <- c(upGene1,upGene2,upGene3,downGene1,downGene2,downGene3)
matrixSubCpm <- cpm(matrixSub)
exprHM <- matrixSubCpm[diffGene,]

ageSplit <- metaSub$Age
exprHM <- t(scale(t(exprHM)))
Heatmap(exprHM,column_split = ageSplit,show_row_names = F,show_column_names = F,km=10)

intersectGene <-c(up_intersect,down_intersect)


hm <- Heatmap(exprHM[intersectGene,],column_split = ageSplit,show_row_names = F,show_column_names = F,km = 4)

hm <- Heatmap(exprHM,column_split = ageSplit,show_row_names = F,show_column_names = F,km=10)
hm <- draw(hm)
clusterlist = row_order(hm)
countData <- hm@ht_list$matrix_18@matrix
clu_df <- lapply(names(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(countData[clusterlist[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)
rownames(clu_df) <- clu_df$GeneID

write.csv(clu_df,"result/12.16_age_lepr_stem/cluster.csv")


UpSetR::upset(fromList(list(list1=downGene1,list2=downGene2,list3=downGene3)))

form <- "~ Organ + Age + Stage +Origin+Bone.Forming.Methods+(1|Project)" 
vobjDreamList=list()
for (i in unique(meta$C7_named)){
  metaSub=meta[meta$C7_named==i,]
  matrixSub=countMatrix[meta$C7_named==i,]
  geneExpr = DGEList( matrixSub)
  param = SnowParam(progressbar=TRUE) 
  register(param)
  geneExpr = calcNormFactors( geneExpr ) 
  vobjDream = voomWithDreamWeights(geneExpr, form, meta)
  vobjDreamList[[i]]=vobjDream
}
fitmm <- lapply(vobjDreamList,dream,formula=form,data=meta)
geneExpr$samples
output_dir="../3.9_wt_integrate/result/10.2_dream_deg/"
write.fit(fitmm$Chondro, file=paste(output_dir, "/chondro_output.tsv", sep=""))
write.fit(fitmm$Fibroblast, file=paste(output_dir, "/Fibroblast_output.tsv", sep=""))
write.fit(fitmm$`Lepr+ BMSC`, file=paste(output_dir, "/lepr_output.tsv", sep=""))
write.fit(fitmm$`Ly6a+ MSC`, file=paste(output_dir, "/ly6a_output.tsv", sep=""))
write.fit(fitmm$MSC, file=paste(output_dir, "/msc_output.tsv", sep=""))
write.fit(fitmm$Ob, file=paste(output_dir, "/ob_output.tsv", sep=""))
write.fit(fitmm$Pericyte, file=paste(output_dir, "/pericyte_output.tsv", sep=""))
for (i in unique(meta$C7_named)){
  metaSub=meta[meta$C7_named==i,]
  matrixSub=countMatrix[meta$C7_named==i,]
  geneExpr = DGEList( matrixSub)
  geneExpr = calcNormFactors( geneExpr ) 
  write.csv(geneExpr,paste0(output_dir,i,"_geneExpr.csv"))
}


geneExpr = DGEList( countMatrix ) #[isexpr,] ) # creates object in certain format, i.e. a counts object 
# ($counts), that has the count matrix, and a samples object ($samples) that has samples as rows,
# and for each sample: group assignments, lib.sizes, and norm factors
param = SnowParam(progressbar=TRUE) 
register(param)
geneExpr = calcNormFactors( geneExpr ) 
vobjDream = voomWithDreamWeights(geneExpr, form, meta)

#== decouper-------------------
library(decoupleR)
gmt_to_decoupler <- function(pth) {
  pathways <- list()
  
  con <- file(pth, "r")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    fields <- strsplit(line, "\t")[[1]]
    name <- fields[[1]]
    genes <- fields[-1]
    pathways[[name]] <- genes
  }
  close(con)
  
  data <- do.call(rbind.data.frame, lapply(names(pathways), function(name) {
    data.frame(geneset = rep(name, length(pathways[[name]])), genesymbol = pathways[[name]], stringsAsFactors = FALSE)
  }))
  
  return(data)
}
test <- gmt_to_decoupler("data/m2.all.v2023.2.Mm.symbols (1).gmt")
testorder <- grep("www",test$genesymbol)
M2gmt <- test[-testorder,]

test <- run_aucell(
  matrixSub,
  M2gmt,
  .source="geneset",
  .target="genesymbol"
)
sample_acts_mat <- test %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()
# Scale per feature

sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))
sample_acts_mat_sub <- t(sample_acts_mat)
index_osteo <- grep("OSTEO",rownames(sample_acts_mat_sub))
index_fatty <- grep("ADIPO",rownames(sample_acts_mat_sub))

# Plot
pheatmap(sample_acts_mat[,1:10], border_color = NA, color=my_color, breaks = my_breaks) 
Heatmap(t(sample_acts_mat[,c(index_osteo,index_fatty)]),column_split = ageSplit,show_row_names = T,show_column_names = F)




form <- "~  Age_num + Tissue.Specific." 
vobjDreamList=list()


geneExpr = DGEList( matrixSub)

geneExpr = calcNormFactors( geneExpr ) 
vobjDream = voomWithDreamWeights(geneExpr, form, metaSub)
fitmm <- dream(vobjDream,formula=form,data=metaSub)
test <- eBayes(fitmm)
write.fit(fitmm,"result/12.16_age_lepr_stem/age_fit.csv")
fitmmTab <- topTable(test, number=Inf)
vobjDreamList <- list()
vobjDreamList[[1]]=vobjDream
test2 <- decideTests(fitmm)
test2 <- as.data.frame(test2)
fitmmTab <- topTable(test, number=Inf,coef = 1)
fitmm_list <- lapply(vobjDreamList,dream,formula=form,data=metaSub)
write.fit(test2 ,"result/12.16_age_lepr_stem/age_fit.csv")


#== dreamlet analysis------------------
library(Seurat)
library(dreamlet)
wtMeta <- readRDS("../important_processed_data/9.10_wt_meta.Rds")
wt_intergrate$Age <- wtMeta$Age

sa


wt_intergrate$Age <- factor(wt_intergrate$Age,levels = c("Organogenesis stage", "Fetal stage", "Postnatal", 
                                                         "Young Adult","Adult",  "Old"))
wt_intergrate$Age_num <- as.numeric(wt_intergrate$Age)
stemFilter <- wt_intergrate[,wt_intergrate$C7_named%in%c("MSC", "Ly6a+ MSC", "Lepr+ BMSC", "Fibroblast", "Chondro")]
filterCell1 <- colnames(stemFilter)[stemFilter$Organ!="Limb_adult"&stemFilter$C7_named=="Lepr+ BMSC"]
filterCell2 <- colnames(stemFilter)[stemFilter$C19_named=="Pre-ob"]
filterCell3 <- colnames(stemFilter)[stemFilter$Organ=="Limb_adult"&stemFilter$C7_named=="MSC"]
filterCell <- unique(c(filterCell1,filterCell2,filterCell3))
stemFilter <- stemFilter[,!colnames(stemFilter) %in% filterCell]
sce <- as.SingleCellExperiment(stemFilter)
#sce$Age_num <- wt_intergrate$Age_num 
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",
                            cluster_id = "C7_named",
                            sample_id = "Sample",
                            verbose = FALSE
)
assayNames(pb)

res.proc <- processAssays(pb, ~Age_num, min.count = 5)
plotVoom(res.proc)
vp.lst <- fitVarPart(res.proc, ~Age_num)
genes <- vp.lst$gene[2:4]
plotPercentBars(vp.lst[vp.lst$gene %in% genes, ])
plotVarPart(vp.lst, label.angle = 60)
res.dl <- dreamlet(res.proc, ~Age_num)
coefNames(res.dl)
res.dl
plotVolcano(res.dl, coef = "Age_num")
genes <- c("Apoe", "Bglap","Sp7","Fbn2","Cd200","Aspn","Postn","Cxcl12","Mfap4","Ovol2")
plotGeneHeatmap(res.dl, coef = "Age_num", genes = genes)
topGene <- topTable(res.dl, coef = "Age_num",number = Inf)%>%as.data.frame()

diffGene <- topGene%>%filter(adj.P.Val < 0.00000000001)%>%
  select("ID")%>%unlist()%>%
  unique()
plotGeneHeatmap(res.dl, coef = "Age_num", genes = diffGene[1:100])
topGeneSplit <- topGene %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(up_down = ifelse(logFC > 0, "up", "down"))%>%
  mutate(assay_updown = paste0(assay,"_",up_down))%>%
  filter(!assay%in%c("Ob","Pericyte"))%>%
  group_split(assay_updown, keep = FALSE)

topTableSplit <- split(topGene,topGene$assay)
View(topTableSplit$`Lepr+ BMSC`)

plotGeneHeatmap(res.dl, coef = "Age_num", genes = genes)
topGeneName <- lapply(topGeneSplit,select, ID)
topGeneName <- lapply(topGeneName,unlist)
names(topGeneName) <- sapply(topGeneSplit,function(x)  unique(paste0(x$assay,"_",x$up_down)))
UpSetR::upset(fromList(topGeneName))






tab <- topGene
assays = assayNames(res.dl)
tab <- tab[, c("assay", "ID", "z.std", 
                                "adj.P.Val"), drop = FALSE]
tab$ID <- factor(tab$ID)
tab <- as.data.frame(tab[tab$assay %in% assays, , drop = FALSE])
tab <- droplevels(tab)
tab$assay <- factor(tab$assay, assays)
if (nrow(tab) == 0) 
  stop("No assays retained")
assay <- ID <- z.std <- adj.P.Val <- NULL
tab <- as.data.frame(complete(tab, assay, ID))
zmax <- max(abs(tab$z.std), na.rm = TRUE)

ncol <- length(unique(tab$assay))
nrow <- length(unique(tab$ID))
aspect.ratio <- nrow/ncol
tab_bk <- tab
tab <- tab_bk
tab$z.std[tab$adj.P.Val>0.05]=NA
tab <- tab[,-4]
df_wide <- pivot_wider(tab, names_from = assay, values_from = z.std)
df_wide <- df_wide[rowSums(is.na(df_wide))<4,]
df_wide[is.na(df_wide)] <- 0
same_rows <- apply(df_wide[2:6], 1, function(row) all(row >= 0,na.rm = T)| all(row <= 0,na.rm = T))
gene_up <- df_wide$ID[apply(df_wide[2:6], 1, function(row) all(row >= 0,na.rm = T))]
gene_down <- df_wide$ID[apply(df_wide[2:6], 1, function(row) all(row <= 0,na.rm = T))]
geneOther <- setdiff(df_wide$ID, c(gene_up,gene_down))

geneLengthDf <- data.frame("concordant up-reg"= length(gene_up),
                           "concordant down-reg"= length(gene_down),
                           "divergent reg"= length(geneOther))
geneLengthDf <- t(geneLengthDf)
geneLengthDf <- geneLengthDf%>%as.data.frame%>%rownames_to_column("text")
# Plotting with ggplot2
colnames(geneLengthDf) <- c("Category", "Count")
ggplot(geneLengthDf, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +  # You can choose a different theme based on your preference
  labs(title = "Distribution of Gene Categories", y = "Count", x = "Category") +scale_fill_npg()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust text size here
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.title = element_text(size = 14)  # Adjust legend title size
  )
ggsave("result/12.19_figure2_OPC_plot/12.23_gene_reg_cat.pdf",width = 8,height = 8)

write.csv(gene_up,"result/12.16_age_lepr_stem/common_up.csv")
write.csv(gene_down,"result/12.16_age_lepr_stem/common_down.csv")





library(zenith)
C2 <- get_MSigDB("C2",org = "mmu",to = "SYMBOL")
res_zenith <- zenith_gsa(res.dl, coef = "Age_num", C2)
plotZenithResults(res_zenith, 5, 1)
gs <- unique(res_zenith$Geneset[res_zenith$FDR < 0.01])
df <- res_zenith[res_zenith$Geneset %in% gs, ]
plotZenithResults(df, Inf, Inf)


df_wide_pre <- res_zenith[c("assay","Geneset","delta")]
df_function_wide <- pivot_wider(df_wide_pre, names_from = assay, values_from = delta)
df_function_wide$mean <- rowMeans(df_function_wide[2:6],na.rm = T)
write.csv(df_function_wide,"result/12.16_age_lepr_stem/C2_function.csv")
go.gs <- get_GeneOntology("BP", to = "SYMBOL",org = "mmu")

# Run zenith gene set analysis on result of dreamlet
res_zenith_go <- zenith_gsa(res.dl, coef = "Age_num", go.gs)
plotZenithResults(df, Inf, Inf)

plotZenithResults(res_zenith, 5, 1)
select_name <- c("M26899_REACTOME_TRANSLESION_SYNTHESIS_BY_POLH","M27673_REACTOME_MITOTIC_SPINDLE_CHECKPOINT",
                 "M5583_REACTOME_MRNA_SPLICING_MINOR_PATHWAY","M766_KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM","M2059_FERRANDO_HOX11_NEIGHBORS",
                 "M23_PID_WNT_NONCANONICAL_PATHWAY","M6866_MISSIAGLIA_REGULATED_BY_METHYLATION_DN",
                 "M2144_DEMAGALHAES_AGING_UP","M1597_BURTON_ADIPOGENESIS_PEAK_AT_2HR","M54_PID_IL12_2PATHWAY","M27195_REACTOME_SIGNALING_BY_LEPTIN",
                 "M1309_GROSS_HYPOXIA_VIA_HIF1A_DN","M39645_WP_OXIDATIVE_STRESS_RESPONSE","M2602_BIOCARTA_RANKL_PATHWAY")

returnZenithDf <- function(df,gs){
  df <- res_zenith
  delta = se = PValue = tstat = assay = FDR = Geneset = NULL

  df$tstat = with(df, delta/se)

  grd = expand.grid(assay = unique(df$assay), coef = unique(df$coef))
  # gs = apply(grd, 1, function(x) {
  #   df_sub = df[(df$assay == x[1]) & (df$coef == x[2]), 
  #   ]
  #   tstat_sort = sort(df_sub$tstat)
  #   cutoff1 = ifelse(nbottom > 0, tstat_sort[nbottom], -Inf)
  #   cutoff2 = ifelse(ntop > 0, tail(tstat_sort, ntop)[1], 
  #                    Inf)
  #   idx = (df_sub$tstat <= cutoff1) | (df_sub$tstat >= cutoff2)
  #   df_sub$Geneset[idx]
  # })
  # gs = unique(unlist(gs))
  M = dcast(df[df$Geneset %in% gs, ], assay + coef ~ Geneset, 
            value.var = "tstat")
  annot = M[, seq(1, 2)]
  M = as.matrix(M[, -seq(1, 2)])
  rownames(M) = annot$assay
  success = tryCatch({
    hcl2 <- hclust(dist(t(M)))
    TRUE
  }, error = function(e) FALSE)
  if (!success) {
    M_zero = M
    i = which(is.na(M_zero))
    if (length(i) > 0) 
      M_zero[i] = 0
    hcl2 <- hclust(dist(t(M_zero)))
  }
  
  data = df[df$Geneset %in% gs, ]
  data$Geneset = factor(data$Geneset, hcl2$labels[hcl2$order])
  data = as.data.frame(complete(data, Geneset, assay))

  return(data)
  
}


selectZenith <- returnZenithDf(res_zenit,select_name)
ggplot(selectZenith, aes( assay,Geneset, fill = tstat, 
                 label = ifelse(FDR < 0.05, "*", ""))) + geom_tile() + 
  theme_classic() + scale_fill_gradient2("t-statistic", 
                                         low = "blue", mid = "white", high = "red") + 
  geom_text(vjust = 1, hjust = 0.5) + xlab("Gene sets") + 
  xlab("")

ggsave("result/12.19_figure2_OPC_plot/12.23_C2_select_pathway.pdf",width = 8,height = 8)
#== decouper-------------------------

library(decoupleR)
gmt_to_decoupler <- function(pth) {
  pathways <- list()
  
  con <- file(pth, "r")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    fields <- strsplit(line, "\t")[[1]]
    name <- fields[[1]]
    genes <- fields[-1]
    pathways[[name]] <- genes
  }
  close(con)
  
  data <- do.call(rbind.data.frame, lapply(names(pathways), function(name) {
    data.frame(geneset = rep(name, length(pathways[[name]])), genesymbol = pathways[[name]], stringsAsFactors = FALSE)
  }))
  
  return(data)
}
test <- gmt_to_decoupler("data/m2.all.v2023.2.Mm.symbols (1).gmt")
testorder <- grep("www",test$genesymbol)
M2gmt <- test[-testorder,]






leprPbMatrix <- pb@assays@data@listData$`Lepr+ BMSC`%>%as.data.frame()

test <- run_aucell(
  leprPbMatrix,
  M2gmt,
  .source="geneset",
  .target="genesymbol"
)
sample_acts_mat <- test %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()
sample_acts_mat <- t(sample_acts_mat)
sample_acts_mat <- t(scale(t(sample_acts_mat)))
selectRows <- grep("STEM_CELL",rownames(sample_acts_mat))
Heatmap(sample_acts_mat[selectRows,])

saveRDS(pb,"processed_data/12.23_fig2_supp/pseudobulk.Rds")
saveRDS(res_zenith,"processed_data/12.23_fig2_supp/MsigDb_zenit_C2.Rds")
saveRDS(tab_bk,"processed_data/12.23_genetab.Rds")
