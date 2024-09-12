# install.packages("msigdbr")
# BiocManager::install("org.Hs.eg.db")
#
# Warning messages:
#   1: In install.packages(...) :
#   installation of package ‘igraph’ had non-zero exit status
# 2: In install.packages(...) :
#   installation of package ‘tidygraph’ had non-zero exit status
# 3: In install.packages(...) :
#   installation of package ‘ggraph’ had non-zero exit status
# 4: In install.packages(...) :
#   installation of package ‘enrichplot’ had non-zero exit status
# 5: In install.packages(...) :
#   installation of package ‘clusterProfiler’ had non-zero exit status
# 6: In install.packages(update[instlib == l, "Package"], l, contriburl = contriburl,  :
#                          installation of package ‘hdf5r’ had non-zero exit status
#                        7: In install.packages(update[instlib == l, "Package"], l, contriburl = contriburl,  :
#                                                 installation of package ‘igraph’ had non-zero exit status
create_gsea_input_file <- function(counts_data, file_loc) {

  sample_names <- colnames(counts_data)[-c(1)]
  description <- grch38$description[match(counts_data$gene, grch38$symbol)]
  counts_data_gsea <- cbind(counts_data, Description = description)
  setcolorder(x = counts_data_gsea, neworder = c("gene", "Description", sample_names))
  colnames(counts_data_gsea)[c(1)] <- c("Name")
  write.table(x = counts_data_gsea,
              file = file_loc,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)

}

create_cls_file <- function(meta_data, group_col_name, file_loc) {

  line1 <- c(nrow(meta_data),
             levels(meta_data[[group_col_name]]) %>% length(),
             1) %>% paste(., collapse = " ")
  line2 <- c("#", levels(meta_data[[group_col_name]])) %>%
    str_replace_all(string = ., pattern = " ", replacement = "_") %>%
    paste(., collapse = " ")
  line3 <- meta_data[[group_col_name]] %>%
    str_replace_all(string = ., pattern = " ", replacement = "_") %>%
    paste(., collapse = " ")
  writeLines(text = c(line1, line2, line3), con = file_loc, sep = "\n")

}
