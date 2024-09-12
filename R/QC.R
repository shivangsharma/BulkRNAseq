counts_distribution <- function(counts_data,
                                meta_data,
                                category,
                                category_order) {

  counts_melted_data <- melt.data.table(data = counts_data,
                                        id.vars = "gene",
                                        variable.name = "sample",
                                        value.name = "count")
  meta_data <- as.data.frame(meta_data)
  rownames(meta_data) <- meta_data[, 1]
  meta_data[, category] <- factor(x = meta_data[, category],
                                  levels = category_order)
  counts_melted_data[, category] <- meta_data[counts_melted_data[, "sample"] %>%
                                                as.matrix() %>%
                                                as.vector(), category]
  log_counts_melted_data <- counts_melted_data
  log_counts_melted_data[, "count"] <- log1p(log_counts_melted_data[, "count"])
  plt <- density_plot(melted_data = log_counts_melted_data,
                      category = category)
  return(plt)

}
pca_plot <- function(counts_data) {

  mat_data <- counts_data[, -c(1)] %>% as.matrix() %>% t()
  colnames(mat_data) <- counts_data[["gene"]]

  pca_mat_data <- PCA(mat_data, scale.unit = TRUE, ncp = 100, graph = TRUE)
  plt_pca <- PCA_plot(pca_mat_data = pca_mat_data,
                      meta_data = meta_data,
                      PCs = list(c(1, 2), c(3, 4), c(5, 6)),
                      color_by = "Treatment",
                      shape_by = ,
                      fill_by = ,
                      label = ,
                      pca_plot_only = FALSE)
  plt[["ScreePlot"]] <- fviz_screeplot(pca_mat_data, addlabels = TRUE, ylim = c(0, 50))
  plt[["PC1_vs_PC2"]] <- PCA_plot(pca_mat_data, meta_data)

}
