sce=zellkonverter::readH5AD("processed_data/11.23_process_null/11.23_process_null.h5ad")
sce=sce[,meta$Stage!="Injury"]
meta=as.data.frame(sce@colData@listData)
anno_null <- data.frame("Cell"=colnames(sce),"Sample"=meta$Sample)
pseudo_null <- meta$dpt_pred
names(pseudo_null) <- colnames(sce)
designPre <- table(meta$Sample,meta$Stage_shuffle)%>%as.data.frame%>%filter(Freq!=0)
designVar_null <- designPre$Var2
design_null=model.matrix(~designVar_null)
rownames(design_null) <- designPre$Var1
Res_null <-
  cellPropTest(
    cellanno = anno_null,
    pseudotime =pseudo_null,
    design =design_null,
    ncores = 15,
    test.type = 'Variable',
    testvar = 2
  )

#== half-----------------------
anno_null <- data.frame("Cell"=colnames(sce),"Sample"=meta$Sample)
pseudo_null <- meta$dpt_pred
pseudo_null <- sample(pseudo_null)
pseudo_null_bk=pseudo_null
names(pseudo_null[1:floor(length(pseudo_null)/2)]) <- sample(names(pseudo_null[1:floor(length(pseudo_null)/2)]))
designPre <- table(meta$Sample,meta$Stage)%>%as.data.frame%>%filter(Freq!=0)
designPre$Var2 <- as.character(designPre$Var2)
designVar_null <- designPre$Var2
design_null=model.matrix(~designVar_null)
rownames(design_null) <- designPre$Var1
Res_null <-
  cellPropTest(
    cellanno = anno_null,
    pseudotime =pseudo_null,
    design =design_null,
    ncores = 15,
    test.type = 'Variable',
    testvar = 2
  )
Res_list=list()
pseudo_list=list()
for (i in 0:10){
  pseudo_null <- pseudo_null_bk
  i_chr <- as.character(i)
  pseudoName <- names(pseudo_null_bk)
  names(pseudo_null)[1:floor(length(pseudo_null)*i/10)] <- sample(pseudoName[1:floor(length(pseudo_null_bk)*i/10)])
  pseudo_list[[i_chr]] <- pseudo_null
  Res_null <-
    cellPropTest(
      cellanno = anno_null,
      pseudotime =pseudo_null,
      design =design_null,
      ncores = 15,
      test.type = 'Variable',
      testvar = 2
    )

  Res_list[[i_chr]]=Res_null$statistics
}

p_list=sapply(Res_list,FUN = function(x) x[[1]])
plot(-log(p_list))
z_list=sapply(Res_list, FUN = function(x) x[[2]])
plot(z_list)
lamian_composition_p=data.frame("p_val"=p_list,'z_val'=z_list)
write.csv(lamian_composition_p,"processed_data/11.23_process_null/11.23_lamian_composition_p.csv")

pseudoName=sapply(pseudo_list, names)%>%as.data.frame()
rownames(pseudoName) <- names(pseudo_null_bk)
write.csv(pseudoName,"processed_data/11.23_process_null/11.23_lamian_composition_name.csv")





#== condiment------------------
head(pseudo_null)
pseudo_null[rownames(meta2)]
meta=as.data.frame(sce@colData@listData)
rownames(meta) <- colnames(sce)

library(condiments)
meta$Stage
prog_res_test <- progressionTest(pseudotime =as.matrix(pseudo_null),cellWeights=as.matrix(rep(1,length(pseudo_null))),conditions = meta$Stage, global = TRUE, lineages = TRUE)
prog_res_list=list()
for (i in 0:10){
  pseudo_null <- pseudo_null_bk
  i_chr <- as.character(i)
  pseudo_null <- pseudo_list[[i_chr]]
  pseudo_null <- pseudo_null[rownames(meta)]%>%as.matrix()
  prog_res_loop <- progressionTest(pseudotime =as.matrix(pseudo_null),cellWeights=as.matrix(rep(1,length(pseudo_null))),conditions = meta$Stage, global = TRUE, lineages = TRUE)
  
  
  prog_res_list[[i_chr]]=prog_res_loop
}
condiment_static_list=sapply(prog_res_list,function(x) x[[2]])
condiment_p_list=sapply(prog_res_list,function(x) x[[3]])
condiment_metrics=data.frame("p_val"=condiment_p_list,'statics'=condiment_static_list)
write.csv(condiment_metrics,"processed_data/11.23_process_null/11.23_condiment_metrics.csv")

#== shuffle sample but not pseudotime-----------------
designPre <- table(meta$Sample,meta$Stage)%>%as.data.frame%>%filter(Freq!=0)
sampleTable <- designPre[c(1,2)]
colnames(sampleTable) <-c("Sample","Stage")
set.seed(1)
sampleTable$Stage_shuffle <- sample(sampleTable$Stage)
stage_shuffle <- sampleTable$Stage_shuffle
names(stage_shuffle) <- sampleTable$Sample
meta$Sample <- as.character(meta$Sample)
condition_shuffle <- stage_shuffle[meta$Sample]
shuffle_condiment_list=list()
condition_list=list()
for (i in 1:10){
  set.seed(i)
  sampleTable$Stage_shuffle <- sample(sampleTable$Stage)
  stage_shuffle <- sampleTable$Stage_shuffle
  names(stage_shuffle) <- sampleTable$Sample
  meta$Sample <- as.character(meta$Sample)
  condition_shuffle <- stage_shuffle[meta$Sample]
  condition_shuffle <- factor(condition_shuffle,levels = unique(condition_shuffle))
  condition_list[[i]] <- condition_shuffle
  prog_res_loop <- progressionTest(pseudotime =as.matrix(meta$dpt_pred),cellWeights=as.matrix(rep(1,length(pseudo_null))),conditions =condition_shuffle, global = TRUE, lineages = TRUE)
  shuffle_condiment_list[[i]] <- prog_res_loop
  
}
shuffle_condiment_static_list=sapply(shuffle_condiment_list,function(x) x[[2]])
shuffle_condiment_p_list=sapply(shuffle_condiment_list,function(x) x[[3]])
condiment_metrics_shuffle=data.frame("p_val"=shuffle_condiment_p_list,'statics'=shuffle_condiment_static_list)
write.csv(condiment_metrics_shuffle,"processed_data/11.23_process_null/11.23_condiment_metrics_shuffle.csv")



condition_df=sapply(condition_list, function(x) x)%>%as.data.frame()
rownames(condition_df) <- colnames(sce)

write.csv(condition_df,"processed_data/11.23_process_null/11.23_shuffle_condition.csv")



#== lamian_shuffle--------------------------
pseudo_null <- meta$dpt_pred
names(pseudo_null) <- colnames(sce)
Res_list_shuffle <- list()
for (i in 1:10){
  meta$Stage_shuffle <- condition_list[[i]]
  designPre <- table(meta$Sample,meta$Stage_shuffle)%>%as.data.frame%>%filter(Freq!=0)
  designPre$Var2 <- as.character(designPre$Var2)
  designVar_null <- designPre$Var2
  design_null=model.matrix(~designVar_null)
  rownames(design_null) <- designPre$Var1
  i_chr <- as.character(i)
  Res_null <-
    cellPropTest(
      cellanno = anno_null,
      pseudotime =pseudo_null,
      design =design_null,
      ncores = 15,
      test.type = 'Variable',
      testvar = 2
    )
  
  Res_list_shuffle[[i_chr]]=Res_null$statistics
}
Res_list_shuffle

meta$Stage_shuffle <- condition_list[[5]]
designPre <- table(meta$Sample,meta$Stage_shuffle,meta$Stage)%>%as.data.frame%>%filter(Freq!=0)
designPre$Var2 <- as.character(designPre$Var2)
designVar_null <- designPre$Var2

meta$random_sample <- "None"
meta$random_sample[meta$Stage=="Development"] <- sample(1:20, sum(meta$Stage=="Development"), replace = T)
meta$random_sample[meta$Stage!="Development"] <- sample(21:40, sum(meta$Stage!="Development"), replace = T)


#== condiments 2-----------------------
shuffle_condiment_list2 <- list()
condition_list2 <- list()
designPre <- table(meta$random_sample,meta$Stage)%>%as.data.frame%>%filter(Freq!=0)
sampleTable <- designPre[c(1,2)]
colnames(sampleTable) <-c("Sample","Stage")
for (i in 1:10){
  set.seed(i)
  sampleTable$Stage_shuffle <- sample(sampleTable$Stage)
  stage_shuffle <- sampleTable$Stage_shuffle
  names(stage_shuffle) <- sampleTable$Sample
  condition_shuffle <- stage_shuffle[meta$random_sample]
  condition_shuffle <- factor(condition_shuffle,levels = unique(condition_shuffle))
  condition_list2[[i]] <- condition_shuffle
  prog_res_loop <- progressionTest(pseudotime =as.matrix(meta$dpt_pred),cellWeights=as.matrix(rep(1,length(pseudo_null))),conditions =condition_shuffle, global = TRUE, lineages = TRUE)
  shuffle_condiment_list2[[i]] <- prog_res_loop
  
}


shuffle_condiment_static_list2=sapply(shuffle_condiment_list2,function(x) x[[2]])
shuffle_condiment_p_list2=sapply(shuffle_condiment_list2,function(x) x[[3]])
condiment_metrics_shuffle2=data.frame("p_val"=shuffle_condiment_p_list2,'statics'=shuffle_condiment_static_list2)
write.csv(condiment_metrics_shuffle2,"processed_data/11.23_process_null/11.23_condiment_metrics_shuffle2_newassign.csv")

condition_df2=sapply(condition_list, function(x) x)%>%as.data.frame()
rownames(condition_df2) <- colnames(sce)
condition_df2$random_sample <- meta$random_sample
write.csv(condition_df2,"processed_data/11.23_process_null/11.23_shuffle_condition_newassign.csv")

#== lamian------------------
pseudo_null <- meta$dpt_pred
names(pseudo_null) <- colnames(sce)
anno_new <- data.frame("Cell"=colnames(sce),"Sample"=meta$random_sample)
Res_list_shuffle_newassign <- list()
for (i in 1:10){
  meta$Stage_shuffle <- condition_list[[i]]
  designPre <- table(meta$random_sample,meta$Stage_shuffle)%>%as.data.frame%>%filter(Freq!=0)
  designPre$Var2 <- as.character(designPre$Var2)
  designVar_null <- designPre$Var2
  design_null=model.matrix(~designVar_null)
  rownames(design_null) <- designPre$Var1
  i_chr <- as.character(i)
  Res_null <-
    cellPropTest(
      cellanno = anno_new,
      pseudotime =pseudo_null,
      design =design_null,
      ncores = 15,
      test.type = 'Variable',
      testvar = 2
    )
  
  Res_list_shuffle_newassign[[i_chr]]=Res_null$statistics
}

p_list_null=sapply(Res_list_shuffle_newassign,FUN = function(x) x[[1]])
plot(-log(p_list_null))
z_list_null=sapply(Res_list_shuffle_newassign, FUN = function(x) x[[2]])
plot(z_list_null)
lamian_composition_null=data.frame("p_val"=p_list_null,'z_val'=z_list_null)
write.csv(lamian_composition_null,"processed_data/11.23_process_null/11.23_lamian_composition_null.csv")



meta$Stage_shuffle <- condition_list[[5]]
designPre <- table(meta$Sample,meta$Stage_shuffle,meta$Stage)%>%as.data.frame%>%filter(Freq!=0)
designPre$Var2 <- as.character(designPre$Var2)
designVar_null <- designPre$Var2

meta$random_sample <- "None"
meta$random_sample[meta$Stage=="Development"] <- sample(1:20, sum(meta$Stage=="Development"), replace = T)
meta$random_sample[meta$Stage!="Development"] <- sample(21:40, sum(meta$Stage!="Development"), replace = T)

