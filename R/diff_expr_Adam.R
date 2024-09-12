# Run DESeq2 algorithm of the data ---------------------------------------------
run_deseq2 <- function(mat_data,
                       meta_data,
                       comparison_variable,
                       factors_comparison_variable,
                       log2FC_threshold = 2,
                       padj_threshold = 0.05) {

  deseq2Data <- DESeqDataSetFromMatrix(
    countData = round(mat_data),
    colData = meta_data,
    design = ~ Treatment)
  deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 3, ]
  register(MulticoreParam(80))
  deseq2Data <- DESeq(deseq2Data)

  deseq2Results <- results(deseq2Data, contrast = c(comparison_variable, factors_comparison_variable))
  DESeq2ResDF <- as.data.frame(deseq2Results)
  #summary(deseq2Results)

  vol_plot <- volcano_plot(DESeq2ResDF = DESeq2ResDF, log2FC_threshold = log2FC_threshold, padj_threshold = padj_threshold, factors = factors_comparison_variable)
  return(list(deseq2_obj = deseq2Data,
              deseq2_res = DESeq2ResDF,
              volcano_plot = vol_plot))

}

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
