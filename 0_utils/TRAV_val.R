validate_local <- function (dataset, RAVmodel, method = "pearson", maxFrom = "PC", 
          level = "max", scale = FALSE) 
{
  if (!is.list(dataset)) {
    if (ncol(dataset) < 8) {
      stop("Provide a study with at least 8 samples.")
    }
  }
  else {
    if (any(lapply(dataset, ncol) < 8)) {
      stop("Provide a study with at least 8 samples.")
    }
    if (level == "all") {
      stop("'level = \"all\"' is not available for a list of datasets.")
    }
  }
  sw <- silhouetteWidth(RAVmodel)
  cl_size <- S4Vectors::metadata(RAVmodel)$size
  avgLoading <- SummarizedExperiment::assay(RAVmodel)
  if (maxFrom == "PC") {
    gene_common <- intersect(rownames(avgLoading), rownames(dataset))
    model_test <- RcppML::nmf(t(dataset[gene_common,]), 15, verbose = F, seed = 1234)
    h <- model_test$h
    colnames(h) <- rownames(dataset)
    colFilter <- rowSums(h)!=0
    h <- h[colFilter,]
    loading_cor <- abs(stats::cor(avgLoading[gene_common,], 
                                  t(h[,gene_common]), use = "pairwise.complete.obs", 
                                  method = method))
    x <- loading_cor

    z <- apply(x, 1, max) %>% as.data.frame
    z$PC <- apply(x, 1, which.max)
    colnames(z)[1] <- "score"
    z$sw <- sw
    z$cl_size <- cl_size
    z$cl_num <- readr::parse_number(rownames(z))
    res <- z
    return(res)

    }
}
