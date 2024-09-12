#' Convert Tempus raw count data to standard matrix for input in DESeq2
convert_tempus_to_matrix <- function(filepath, col_name) {

  # Read Tempus file containing normalized raw counts --------------------------
  obj <- read.csv(file = filepath, header = TRUE) %>%
    as.data.frame()
  patient_ids <- unique(obj$partner_patient_id)


  # Convert into data frame with each patient as a separate column -------------
  obj <- obj[, c("ensembl_gene", "partner_patient_id", col_name)]
  obj <- dcast(obj, formula = ensembl_gene ~ partner_patient_id,
               value.var = col_name)


  # Convert gene names from ENSG format to HGNC (common names) format ----------
  obj$Gene <- convert_ensg_to_hgnc(obj$ensembl_gene)


  # Rearrange data frame -------------------------------------------------------
  obj <- obj[, c("ensembl_gene", "Gene", patient_ids)]


  # Use gene names as rownames and delete Gene column --------------------------
  rownames(obj) <- obj$Gene


  # Create data matrix ---------------------------------------------------------
  rawCounts <- as.matrix(obj[, patient_ids])

}

#' Convert gene names from ENSG format to HGNC format

convert_ensg_to_hgnc <- function(genes_ensg) {

  # Download Ensembl mapping of gene names from ENSG to HGNC format ------------
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  ensg_hgnc_df <- getBM(filters = c("ensembl_gene_id"),
                  attributes= c("ensembl_gene_id", "hgnc_symbol", "description"),
                  values = genes_ensg, mart= mart)


  # Use Ensembl reference to convert ENSG to HGNC in data frame ----------------
  genes_hgnc <- ensg_hgnc_df$hgnc_symbol
  names(genes_hgnc) <- ensg_hgnc_df$ensembl_gene_id
  gene_names <- genes_hgnc[genes_ensg]


  # Replace NAs and "" by ENSG format gene names -------------------------------
  gene_names[gene_names %in% c(NA, "")] <- genes_ensg[gene_names %in% c(NA, "")]


  # Replace recurring Gene names with their ENSG IDs ---------------------------
  # This step should be done only after removing NAs and "" gene names
  gene_names[duplicated(gene_names) %>% which()] <-
    genes_ensg[duplicated(gene_names) %>% which()]


  gene_names

}

read_tcga_dataset <- function(working_dir,
                              unzip,
                              biospecimen_data = FALSE,
                              clinical_data = FALSE,
                              sample_sheet = FALSE,
                              datatype) {

  samples_path <- working_dir

  if(unzip == TRUE) {

    data_zipped <- list.files(pattern = "gdc_download_.*",
               path = working_dir,
               full.names = TRUE)
    samples_path <- file.path(working_dir, basename(data_zipped) %>% str_remove(pattern = ".tar.gz"))
    untar(tarfile = data_zipped, exdir = file.path(working_dir, basename(data_zipped) %>% str_remove(pattern = ".tar.gz")))

  }

  samples <- list.dirs(path = samples_path, full.names = TRUE) %>%
    list.files(path = ., pattern = ".*mirnas.quantification.txt", full.names = TRUE)
  samples_expr <- lapply(X = samples, datatype = datatype, FUN = function(x, datatype) {

    if (datatype == "mRNA") {

      temp <- read.table(file = x, sep = "\t", skip = 6)
      colnames(temp) <- c("gene_id", "gene_name", "gene_type", "unstranded", "stranded_first", "stranded_second", "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded")
      temp <- temp[, c("gene_id", "unstranded")]

      } else if (datatype == "miRNA") {

        temp <- read.table(file = x, sep = "\t", header = TRUE)
        temp <- temp[, c("miRNA_ID", "read_count")]

    }

  })

  names(samples_expr) <- basename(samples)
  samples_expr_names <- names(samples_expr) %>% str_extract(string = ., pattern = "[^.]+")
  if (datatype == "mRNA") {

    by_list <- rep(x = "gene_id", times = length(samples_expr_names) - 1)
    dt <- reduce2(samples_expr, by_list, full_join)
    colnames(dt) <- c("gene_id", samples_expr_names)
    dt$gene_name <- convert_ensg_to_hgnc(dt[, 1] %>% str_extract(pattern = "[^.]+"))
    dt <- dt[, c("gene_id", "gene_name", samples_expr_names)]

  } else if (datatype == "miRNA") {

    by_list <- rep(x = "miRNA_ID", times = length(samples_expr_names) - 1)
    dt <- reduce2(samples_expr, by_list, full_join)
    colnames(dt) <- c("miRNA_ID", samples_expr_names)

  }

  return(dt)

}

read_gtex_dataset <- function(working_dir, unzip) {

  samples <- list.files(path = working_dir,
                        pattern = ".*gct.gz",
                        full.names = TRUE)
  dt <- read.table(file = samples, skip = 2)
  colnames(dt) <- dt[1, ]
  colnames(dt)[c(2, 3)] <- c("gene_id", "gene_name")
  dt <- dt[-1, ]
  dt$id <- NULL
  rownames(dt) <- NULL
  return(dt)

}

cd.. <- function(dir, base = FALSE) {

  if (base == FALSE) {

    new_dir <- str_remove(string = dir, pattern = paste0("/", basename(dir)))

  } else {

    # Need to find right regex
    new_dir <- str_remove(string = dir, pattern = "")

  }

  return(new_dir)

}

read_rsem_samples <- function(working_dir,
                              sample_names = NA,
                              file_names = NA,
                              isoforms = FALSE) {

  if (isoforms == TRUE) {

    extension <- "-trimmed.isoforms.results"

  } else {

    extension <- "-trimmed.genes.results"

  }

  if (!is.na(sample_names)) {

    file_locs <- paste0(sample_names, extension) %>%
      file.path(working_dir, .)

  } else if (!is.na(file_names)) {

    file_locs <- file.path(sample_names, file_names)
    sample_names <- str_remove(string = file_names, pattern = extension)

  } else {

    file_locs <- list.files(path = working_dir,
                            pattern = paste0("*", extension),
                            full.names = TRUE)
    sample_names <- basename(file_locs) %>%
      str_remove(string = ., pattern = extension)
    sprintf("Found %i samples", length(sample_names)) %>% message()

  }

  data_list <- vector(mode = "list", length = length(file_locs))
  names(data_list) <- sample_names
  for (file_loc in file_locs) {

    sample_name <- basename(file_loc) %>%
      str_remove(string = ., pattern = extension)
    sprintf("Reading sample: %s", sample_name) %>% message()
    data_list[[sample_name]] <- read.table(file = file_loc,
                                           header = TRUE,
                                           sep = "\t")

  }
  return(data_list)

}

rsem_samples_list_to_data_table <- function(data_list,
                                            isoforms = FALSE,
                                            counts_type = c("expected_count",
                                                            "TPM",
                                                            "FPKM")) {

  if (isoforms == TRUE) {

    pivot_column <- "transcript_id"

  } else {

    pivot_column <- "gene_id"

  }

  all_species <- Reduce(f = union,
                        x = lapply(X = data_list, FUN = "[", , pivot_column))
  data_list_expanded <- lapply(X = data_list, FUN = function(x) {

    sample_species <- x[, pivot_column]
    missing_sample_species <- all_species[!(all_species %in% sample_species)]
    x_append <- matrix(0,
                       nrow = length(missing_sample_species),
                       ncol = ncol(x)) %>% data.table()
    x_append[, 1] <- missing_sample_species
    colnames(x_append) <- colnames(x)
    x <- rbind(x, x_append)
    x <- x[match(x[, pivot_column], all_species), ]

  })
  counts_table <- matrix(nrow = length(all_species),
                         ncol = length(data_list_expanded) + 1) %>%
    data.table()
  colnames(counts_table) <- c("id", names(data_list_expanded))
  counts_table[, "id"] <- all_species
  for (sample in names(data_list_expanded)) {

    counts_table[, sample] <- data_list_expanded[[sample]][, counts_type] %>%
      as.data.table()

  }
  return(counts_table)

}

ensg_to_gene <- function(id_vec) {

  gene_dict <- grch38$symbol
  names(gene_dict) <- grch38$ensgene
  gene_vec <- gene_dict[id_vec]
  return(gene_vec)

}

combine_duplicate_genes <- function(counts_data) {
  # Provide a data table with first column containing gene names and the name of
  # the first column should be "id". The data table is supplied as an argument
  # to counts_data variable. The output will be a data table with first column
  # containing gene names and named as "gene". The rows with same gene names
  # will be summed together.
  simplified_counts_data <- counts_data %>%
    group_by(c(id)) %>%
    summarize(across(.cols = -id, .fns = sum))
  colnames(simplified_counts_data)[1] <- "gene"
  simplified_counts_data <-
    simplified_counts_data[-(simplified_counts_data$gene %in% c("")), ]
  simplified_counts_data <-
    simplified_counts_data[!is.na(simplified_counts_data$gene), ]
  return(simplified_counts_data %>% as.data.table())

}
