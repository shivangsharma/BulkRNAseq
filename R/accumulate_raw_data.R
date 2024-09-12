bulkRNAseq <- setClass(Class = "bulkRNAseq", slots = c(counts = "data.table", meta.data = "data.frame"))

# Read data into a data table
rsem_output_directory <- "/media/singhlab/B338-CE9D/Cell_line_RNAseq/data/RSEM_quantification"
data_list <- read_rsem_samples(working_dir = rsem_output_directory,
                               sample_names = NA,
                               file_names = NA,
                               isoforms = FALSE)
counts_data <- rsem_samples_list_to_data_table(data_list = data_list,
                                               isoforms = FALSE,
                                               counts_type = "expected_count")
counts_data[, "id"] <-
  ensg_to_gene(id_vec = str_remove(string = counts_data[, "id"] %>%
                                     as.matrix() %>%
                                     as.vector(),
                                   pattern = "\\..*")
  )
counts_data <- combine_duplicate_genes(counts_data)

# Read meta data into a data frame
meta_data <- read_xlsx(path = "/media/singhlab/B338-CE9D/Cell_line_RNAseq/data/23257-01-Sample Key.xlsx")
order_column <- list("Cell line" = c("LnCAP", "DU145", "U87"),
                     "Treatment" = c("Wild Type", "Scramble", "CD276 KO"),
                     "Replicate" = c("1", "2", "3"))
meta_data <- convert_groups_to_factors(meta_data, order_column)
obj <- new(Class = "bulkRNAseq", counts = counts_data, meta.data = meta_data)
saveRDS(object = obj, file = file.path(cd..(rsem_output_directory), "raw_counts_meta_data.rds"))
