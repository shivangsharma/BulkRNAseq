# Run DESeq2 algorithm of the data ---------------------------------------------
run_deseq2 <- function(mat_data,
                       meta_data,
                       comparison_variable,
                       factors_comparison_variable,
                       log2FC_threshold = 2,
                       padj_threshold = 0.05,
                       recalc_sf = FALSE) {

  deseq2Data <- DESeqDataSetFromMatrix(
    countData = round(mat_data),
    colData = meta_data,
    design = ~ Treatment)
  deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 3, ]
  register(MulticoreParam(80))
  if (recalc_sf == TRUE) {

    sf <- calc_size_factors(mat_data)
    sizeFactors(deseq2Data) <- sf

  }
  deseq2Data <- DESeq(deseq2Data)

  deseq2Results <- results(deseq2Data, contrast = c(comparison_variable, factors_comparison_variable))
  DESeq2ResDF <- as.data.frame(deseq2Results)
  #summary(deseq2Results)

  vol_plot <- volcano_plot(DESeq2ResDF = DESeq2ResDF, log2FC_threshold = log2FC_threshold, padj_threshold = padj_threshold, factors = factors_comparison_variable)
  return(list(deseq2_obj = deseq2Data,
              deseq2_res = DESeq2ResDF,
              volcano_plot = vol_plot))

}

calc_size_factors <- function(mat_data, pctile = 75) {


  median_mat_data <- apply(mat_data, MARGIN = 1, FUN = median)
  pctile_val <- quantile(x = median_mat_data, probs = c(0.75))
  pctile_genes <- median_mat_data[median_mat_data >= pctile_val] %>% names()
  red_mat_data <- mat_data[pctile_genes, ]
  row_med_red_mat_data <- apply(red_mat_data, MARGIN = 1, FUN = function(x) {

    log() %>% mean() %>% exp() %>% return()

  }) %>% as.vector(mode = "matrix")
  red_mat_data <- sweep(X = red_mat_data,
                        MARGIN = 2,
                        STATS = row_med_red_mat_data,
                        FUN = "/")
  sf <- apply(X = red_mat_data, MARGIN = 2, FUN = median)
  return(sf)

}
