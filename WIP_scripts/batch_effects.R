library(limma)
library(sva)
library(ruv)

# Batch correct with ComBat-seq
batches <- sapply(X = sample_metadata$Dataset, FUN = switch, "GSE176031" = 3, "GSE185344" = 5, "GSE193337" = 6, "GSE172357" = 4, "GSE120716" = 1, "GSE141445" = 2, USE.NAMES = F)
groups <- sapply(X = sample_metadata$Sample_type, FUN = switch, "Normal" = 1, "Benign" = 2, "BPH" = 3, "Tumor" = 4, USE.NAMES = F)
corrected_data <- ComBat_seq(counts = data_PCgenes[, c(-1)], batch = batches)

# Batch correct with limma



# Batch correct with RUV


# Batch correct with conos

data_pc_norm <- data_pc
data_pc_norm[, c(-1)] <- t(t(data_pc[, c(-1)]) / colSums(data_pc[, c(-1)]))*1e6
data_pc_norm_scale <- data_pc_norm
data_pc_norm_scale[, c(-1)] <- scale(x = t(data_pc_norm[, c(-1)]), 
                            center = TRUE, 
                            scale = TRUE) %>% t()

df_patients <- read_xlsx(path = file.path(working_dir, "Patient_categories.xlsx"), col_names = TRUE) %>% as.data.frame()
df_patients$Sample_type <- as.factor(df_patients$Sample_type)
