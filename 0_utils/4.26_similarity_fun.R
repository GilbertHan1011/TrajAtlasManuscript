downsample_and_predict <- function(reflatent_seurat, scANVI, subsample_size,suggestion=TRUE,threshold=0.2) {
  require(randomForest)
  cell_ranks <- levels(reflatent_seurat)
  downsampled_cells <- c()
  for (i in 1:length(cell_ranks)) {
    current_cells <- WhichCells(reflatent_seurat, idents = cell_ranks[i])
    if (min(length(current_cells), subsample_size) < subsample_size) {
      warning(paste0("The current ident ", cell_ranks[i], " is smaller than sample size, adjust to the ident size"))
    }
    downsampled_cells <- c(downsampled_cells, sample(current_cells, size = min(length(current_cells), subsample_size), replace = FALSE))
  }
  downsample_seurat <- reflatent_seurat[, downsampled_cells]
  
  rf_model <- randomForest(downsample_seurat@reductions[[scANVI]]@cell.embeddings, y = Idents(downsample_seurat))
  votes <- rf_model$votes
  identmatrix <- as.data.frame(Idents(downsample_seurat))
  names(identmatrix) <- "ident"
  votes <- cbind(votes, identmatrix)
  summary <- votes %>%
    group_by(ident) %>%
    summarize(across( everything(),.fns = mean))%>%
    column_to_rownames("ident")%>%
    as.matrix()
  summary <- apply(summary, 2, function(x) (x)/max(x))
  diag(summary) <- 0
  summary <- t(summary)
  merge_submit(summary,threshold)
  return(summary)
  
}


#== multilevel_funciton---------------------------------------
multilevel_predict <- function(reflatent_seurat, scANVI, level,pred_level,subsample_size,suggestion=TRUE,threshold=0.2) {
  require(randomForest)
  downsampledf <- data.frame()
  pred <- c()
  for (i in level){
    Idents(reflatent_seurat) <- reflatent_seurat@meta.data[[i]]
    cell_ranks <- levels(reflatent_seurat)
    downsampled_cells <- c()
    for (i in 1:length(cell_ranks)) {
      current_cells <- WhichCells(reflatent_seurat, idents = cell_ranks[i])
      if (min(length(current_cells), subsample_size) < subsample_size) {
        warning(paste0("The current ident ", cell_ranks[i], " is smaller than sample size, adjust to the ident size"))
      }
      downsampled_cells <- c(downsampled_cells, sample(current_cells, size = min(length(current_cells), subsample_size), replace = FALSE))
    }
    downsample_seurat <- reflatent_seurat[, downsampled_cells]
    downsampledf <- rbind(downsampledf,downsample_seurat@reductions[[scANVI]]@cell.embeddings)
    pred <- c(pred,as.character(Idents(downsample_seurat)))
  }
  pred <- as.factor(pred)
  rf_model <- randomForest(downsampledf, y = pred)
  prediction <- predict(rf_model,reflatent_seurat@reductions[[scANVI]]@cell.embeddings,type = "prob")
  identmatrix <- as.data.frame(reflatent_seurat@meta.data[[pred_level]])
  names(identmatrix) <- "ident"
  prediction <- cbind(prediction, identmatrix)
  
  return(prediction)
}


#== predict ident function----------------------

ident_predict <- function(prediction,pred_level,threshold){
  identmatrix <- as.data.frame(pred_level)
  names(identmatrix) <- "ident"
  prediction <- cbind(prediction, identmatrix)
  summary <- prediction %>%
    group_by(ident) %>%
    summarize(across( everything(),.fns = mean))%>%
    column_to_rownames("ident")%>%
    as.matrix()
  summary <- t(summary)
  
  # identify indices of elements that meet the condition
  idx <- which(summary > threshold, arr.ind = TRUE)
  
  # extract row and column names
  rows <- rownames(summary)[idx[, 1]]
  cols <- colnames(summary)[idx[, 2]]
  
  # print row and column names
  for (i in 1:length(rows)) {
    cat("Value", summary[rows[i], cols[i]], "in", cols[i], "column and", rows[i], "row exceeds the threshold of", threshold, "\n")
  }
  
  return(summary)
}

#== this function is designed to reminder user what cluster can be merged-------------
merge_submit <- function(mat,threshold){
  for (i in rownames(mat)){
    for (j in colnames(mat)){
      if (mat[i,j]>threshold){
        message(paste0(i," and ",j,' should be merged'))
      }
    }
  }
}

#== this function is to show what ident have been changed after update annotation------------------
show_ident_change <- function(seurat, previous_anno, latter_anno) {
  df <- seurat@meta.data[c(previous_anno, latter_anno)]
  df_ident <- df %>%
    group_by({{previous_anno}}) %>%
    unique() %>%
    column_to_rownames(previous_anno)
  for (i in rownames(df_ident)) {
    if (df_ident[i, ] != i) {
      print(paste0("Change ", i, " to ", df_ident[i, ]))
    } else {
      df_ident <- subset(df_ident, row.names(df_ident) != i)
    }
  }
  return(df_ident)
}
