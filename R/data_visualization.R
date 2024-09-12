density_plot <- function(melted_data, category) {

  cols <- vector(mode = "character", length = 0)
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  for(i in c(1:nrow(qual_col_pals))) {

    cols <- c(cols, brewer.pal(n = qual_col_pals[i, ]$maxcolors,
                               name = rownames(qual_col_pals)[i]))

  }
  n_cols <- melted_data[[(category)]] %>% levels() %>% length()
  col_vector <- cols[n_cols]
  ggplot(data = melted_data,
         mapping = aes(x = count, color = !!sym(category))) +
    geom_density()

}

PCA_plot <- function(pca_mat_data,
                     meta_data,
                     PCs = list(c(1, 2), c(3, 4), c(5, 6)),
                     color_by = NA,
                     shape_by = NA,
                     fill_by = NA,
                     label_by = NA,
                     pca_plot_only = FALSE) {

  pca_df <- pca_mat_data$ind$coord %>% as.data.frame()
  pca_df <- cbind(pca_df, meta_data[, -c(1)])

  plt_PCA <- vector(mode = "list", length = length(PCs))
  for(axes in c(1:length(PCs))) {

    plt_PCA[[paste0("PC", PCs[[axes]][1], "_vs_", PCs[[axes]][2])]] <-
      ggplot(data = pca_df) +
      geom_point(mapping = aes(x = !!sym(paste0("Dim.", PCs[[axes]][1])),
                               y = !!sym(paste0("Dim.", PCs[[axes]][2])),
                               color = !!sym(color_by),
                               fill = !!sym(fill_by),
                               shape = !!sym(shape_by),
                               label = !!sym(label_by)
                               )
                 ) +
      geom_text(hjust = 0, nudge_x = 0.05, size = 3) +
      geom_text_repel() +
      labs(x = PCs[[axes]][2],
           y = PCs[[axes]][1],
           title = paste0("PC", PCs[[axes]][1], "_vs_", PCs[[axes]][2])
           )
      theme_classic()
  }

  if(pca_plot_only) {

    return(plt_PCA)

  } else {

    plt_PCA[["Screeplot"]] <- fviz_screeplot(pca_mat_data, addlabels = TRUE)


  }

  design <- "
  12
  33
  45
  66
  78
"
  plt <- FP[[1]] + FP[[1]] +  plot_spacer() + FP[[4]] + FP[[5]] + plot_spacer() + FP[[6]] + FP[[6]] + plot_layout(design = design, guides = 'collect', heights = c(10, -1.5, 10, -1.5, 10)) & theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"))
  plt

}

RLE_plot <- function(data, show_legend = FALSE, show_outliers = FALSE) {
  # Provie a matrix of genesx samples, with gene names as rownames and
  # sample names as column names
  library(matrixStats)
  if (is.matrix(data)) {

    RLE_mat <- data

  } else {

    RLE_mat <- as.matrix(data)

  }
  RLE_mat <- sweep(x = RLE_mat, MARGIN = 1, STATS = matrixStats::rowMedians(RLE_mat), FUN = "-")
  RLE_df <- as.data.frame(RLE_mat)
  RLE_df <- melt(RLE_df, variable.name = "Patient_ID", value.name = "RLE")

  whiskers <- RLE_df %>%
    group_by(Patient_ID) %>%
    summarize(lower = quantile(RLE, 0.25) - 1.5*IQR(RLE),
              upper = quantile(RLE, 0.75) + 1.5*IQR(RLE),
              median = median(RLE),
              RLE = mean(RLE))
  whiskers$discrete_x <- c(1:dim(whiskers)[1])
  ylim_abs <- max(abs(whiskers$lower), whiskers$upper)

  plot_output <- ggplot(RLE_df, aes(x = Patient_ID, y = RLE, fill = Patient_ID)) +
    scale_x_discrete() +
    geom_segment(data = whiskers, aes(x = discrete_x - 0.2, xend = discrete_x + 0.2, y = upper, yend = upper)) +
    geom_segment(data = whiskers, aes(x = discrete_x - 0.2, xend = discrete_x + 0.2, y = lower, yend = lower)) +
    geom_linerange(data = whiskers, aes(x = Patient_ID, ymin = lower, ymax = upper),
                   size = 0.2, linetype = "dashed") +
    geom_boxplot(width = 0.4, color = "black", show.legend = show_legend,
                 outlier.shape = if (show_outliers) "circle" else NA, coef = 0) +
    scale_y_continuous(limits = c(-ylim_abs, ylim_abs)) +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = "black", linewidth = 2),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
    ) +
    labs(x = "Samples", y = "Relative expression")
  print(plot_output)

}

# Supply DESeq2 output
volcano_plot <- function(DESeq2ResDF, log2FC_threshold, padj_threshold, factors, genes_of_interest = c()) {

  # Provide a data frame with gene names as row names and c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj") as column names
  max.limit = ceiling(max(DESeq2ResDF$log2FoldChange))
  min.limit = floor(min(DESeq2ResDF$log2FoldChange))
  limit = min(c(abs(max.limit),abs(min.limit)))
  DESeq2ResDF$gene_symbol <- rownames(DESeq2ResDF)
  upregulated_label <- paste0("Adj. p < ", as.character(padj_threshold), " and FC > ", as.character(2^log2FC_threshold))
  downregulated_label <- paste0("Adj. p < ", as.character(padj_threshold), " and FC < ", as.character(2^(-log2FC_threshold)))
  no_change_label <- "No change"
  DESeq2ResDF$diffExp <- no_change_label
  DESeq2ResDF$diffExp[DESeq2ResDF$log2FoldChange > log2FC_threshold &
                        DESeq2ResDF$padj < padj_threshold] <- upregulated_label
  DESeq2ResDF$diffExp[DESeq2ResDF$log2FoldChange < -log2FC_threshold &
                        DESeq2ResDF$padj < padj_threshold] <- downregulated_label
  DESeq2ResDF$diffExp <- factor(x = DESeq2ResDF$diffExp, levels = c(downregulated_label, no_change_label, upregulated_label))

  DESeq2ResDF$Gene <- ""
  if (length(genes_of_interest) == 0) {

    up_down_genes_loc <- which(DESeq2ResDF$diffExp %in% c(downregulated_label, upregulated_label))
    DESeq2ResDF$Gene[up_down_genes_loc] <- DESeq2ResDF$gene_symbol[up_down_genes_loc]

  } else {

    genes_of_interest_rows <- DESeq2ResDF$gene_symbol %in% genes_of_interest
    DESeq2ResDF$Gene[genes_of_interest_rows] <- DESeq2ResDF$gene_symbol[genes_of_interest_rows]

  }

  xlabeltext <- paste0('Up in ', factors[2], ' <--- log2(FC) ---> Up in ', factors[1])
  cols_assignment <- c("#242daf", "gray", "#B2182B")
  names(cols_assignment) <- c(downregulated_label, no_change_label, upregulated_label)
  p <- ggplot(data=DESeq2ResDF, aes(x = log2FoldChange,
                                    y = -log10(padj),
                                    color = diffExp,
                                    label = Gene)) +
    geom_point() +
    lims(x = c(- limit, limit)) +
    theme_minimal() +
    geom_text_repel(max.overlaps = 10000, segment.color = ifelse(length(genes_of_interest) == 0, NA, "black"), color = "black") +
    scale_color_manual(values = cols_assignment) +
    geom_vline(xintercept=c(-log2FC_threshold, log2FC_threshold),
               col="black",
               linetype = "dashed") +
    geom_hline(yintercept=-log10(padj_threshold),
               col="black",
               linetype = "dashed") +
    labs(x = xlabeltext, y = '-log10(adjusted p)', color = 'Adj. p and FC:') +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=16),
          axis.title=element_text(size=16),
          legend.position="bottom")
  # p2 <- ggplot_build(p)
  # breaks_x <- p2$layout$panel_params[[1]]$x$breaks
  # breaks_x <- breaks_x[!is.na(breaks_x)]
  # breaks_y <- p2$layout$panel_params[[1]]$y$breaks
  # breaks_y <- breaks_y[!is.na(breaks_y)]
  # p <- p +
  #   scale_x_continuous(breaks = sort(c(-log2FC_threshold, log2FC_threshold, breaks_x))) +
  return(p)

}


# Supply dataframe with each column containing the variables to be plotted
# (x-axis) and one additional column (1st column) for groups in box
# plot (if applicable)
# Currently, grouped plots only are supported
box_plot <- function(df_expr,
                     x_axis_name,
                     y_axis_name,
                     group_col_name = "",
                     col_fill,
                     show_outliers = FALSE,
                     show_legend = FALSE,
                     group_plot = TRUE,
                     whiskers_rep = c("IQR", "sd", "95CI", "sm"),
                     p_val = FALSE) {

  if (group_plot == TRUE) {

    df_expr <- reshape2::melt(data = df_expr,
                              id = x_axis_name,
                              variable.name = group_col_name,
                              value.name = y_axis_name)

    whiskers <- df_expr %>%
      group_by(!!sym(x_axis_name), !!sym(group_col_name)) %>%
      summarize(lower = quantile(!!sym(y_axis_name), 0.25) - 1.5*IQR(!!sym(y_axis_name)),
                upper = quantile(!!sym(y_axis_name), 0.75) + 1.5*IQR(!!sym(y_axis_name)),
                median = median(!!sym(y_axis_name)),
                Mean = mean(!!sym(y_axis_name))) %>%
      as.data.frame()
    box_width_per_box <- length(unique(whiskers[, group_col_name]))

  } else {

    whiskers <- df_expr %>%
      group_by(!!sym(x_axis_name)) %>%
      summarize(lower = quantile(!!sym(y_axis_name), 0.25) - 1.5*IQR(!!sym(y_axis_name)),
                upper = quantile(!!sym(y_axis_name), 0.75) + 1.5*IQR(!!sym(y_axis_name)),
                median = median(!!sym(y_axis_name)),
                Mean = mean(!!sym(y_axis_name))) %>%
      as.data.frame()
    box_width_per_box <- 1

  }

  box_width <-0.8
  line_width <- box_width
  discrete_x <- vector(mode = "numeric")
  unique_x_labels <- whiskers[, x_axis_name] %>% as.character() %>% unique()
  x_center_loc <- c(1:length(unique_x_labels))
  for (n_label in x_center_loc) {

    x_label <- unique_x_labels[n_label]
    whiskers_subset <- whiskers[whiskers[, x_axis_name] %in% x_label, ]
    n_boxes <- dim(whiskers_subset)[1]

    if (n_boxes%%2 == 1) {

      lower_seq <- seq(from = n_label, length.out = (n_boxes - 1)/2, by = - box_width/box_width_per_box)
      upper_seq <- seq(from = n_label, length.out = (n_boxes - 1)/2, by = box_width/box_width_per_box)
      discrete_x <- c(discrete_x, c(lower_seq, n_label, upper_seq))

    } else {

      lower_seq <- seq(from = n_label - box_width/2/box_width_per_box, length.out = n_boxes/2, by = - box_width/box_width_per_box)
      upper_seq <- seq(from = n_label + box_width/2/box_width_per_box, length.out = n_boxes/2, by = box_width/box_width_per_box)
      discrete_x <- c(discrete_x, c(lower_seq, upper_seq))

    }

  }
  whiskers$discrete_x <- discrete_x
  ylim_max <- max(whiskers$upper)
  ylim_min <- min(whiskers$lower)
  ylim_max <- ylim_max + (ylim_max - ylim_min)*0.1
  ylim_min <- ylim_min - (ylim_max - ylim_min)*0.1

  fill_var <- ifelse(group_plot == TRUE, group_col_name, x_axis_name)
  plot_output <- ggplot(df_expr,
                        aes(x = !!sym(x_axis_name),
                            y = !!sym(y_axis_name),
                            fill = !!sym(fill_var)
                            )
                        ) +
    scale_x_discrete() +
    geom_segment(data = whiskers,
                 aes(x = discrete_x - box_width/4/box_width_per_box,
                     xend = discrete_x + box_width/4/box_width_per_box,
                     y = upper,
                     yend = upper),
                 linewidth = line_width,
                 inherit.aes = FALSE) +
    geom_segment(data = whiskers,
                 aes(x = discrete_x - box_width/4/box_width_per_box,
                     xend = discrete_x + box_width/4/box_width_per_box,
                     y = lower,
                     yend = lower),
                 linewidth = line_width,
                 inherit.aes = FALSE) +
    geom_linerange(data = whiskers,
                   aes(x = discrete_x,
                       ymin = lower,
                       ymax = upper),
                   linewidth = line_width,
                   linetype = "dashed",
                   inherit.aes = FALSE) +
    geom_boxplot(width = box_width,
                 lwd = line_width,
                 color = "black",
                 show.legend = show_legend,
                 outlier.shape = ifelse(show_outliers, "circle", NA),
                 coef = 0,
                 fatten = box_width/2) +
    scale_fill_manual(values = col_fill) +
    scale_y_continuous(limits = c(ylim_min, ylim_max)) +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white",
                                          color = "black",
                                          linewidth = 2)
    )
  print(plot_output)

}

umap_plot <- function(object, dims = c(1, 2), cells = NULL, cols = NULL,
                      pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL,
                      shape.by = NULL, order = NULL, shuffle = FALSE, seed = 1,
                      label = FALSE, label.size = 4, label.color = "black", label.box = FALSE,
                      repel = FALSE, alpha = 1, cells.highlight = NULL, cols.highlight = "#DE2D26",
                      sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE,
                      raster = NULL, raster.dpi = c(512, 512)) {



}
