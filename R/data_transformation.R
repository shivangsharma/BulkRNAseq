unfactor_data <- function(data, column_vec = NA) {

  if (is.data.frame(data)) {

    if (is.na(column_vec)) {

      data[, 1] <- data[, 1] %>% as.matrix() %>% as.vector()
      return(data)

    } else {

      for(column in column_vec) {

        data[, column] <- data[, column] %>% as.matrix() %>% as.vector()

      }
      return(data)

    }

  } else if (is.data.table(data)) {

    if (is.na(column_vec)) {

      data[, 1] <- data[, 1] %>% as.matrix() %>% as.vector()
      return(data)

    } else {

      for(column in column_vec) {

        data[, column] <- data[, column] %>% as.matrix() %>% as.vector()

      }
      return(data)

    }

  } else if (is.matrix(data)) {

    rows_n <- nrow(data)
    data <- matrix(data = as.character(data), nrow = rows_n)

  } else if (is.factor(data)) {

    data <- data %>% as.matrix() %>% as.vector()
    return(data)

  } else if (!is.data.frame(data) & !is.data.table(data) & is.list(data)) {

    for (column in names(column_vec)) {

      data[[column]] <- data[[column]] %>%  as.matrix() %>% as.vector()

    }

  } else {

    stop("Unidentified data type. Please supply data frame, data table, matrix, list or vector type variable as data")

  }

}

df_transpose <- function(df) {

  col_names <- rownames(df)
  row_names <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  df <- as.data.frame(df, row.names = row_names)
  colnames(df) <- col_names
  return(df)

}

convert_groups_to_factors <- function(meta_data, order_column) {

  for (column in names(order_column)) {

    meta_data[, column] <- factor(meta_data[[column]],
                                  levels = order_column[[column]])

  }

  return(meta_data)

}

subset_data <- function(counts_data, meta_data, subset_column) {

  for (column in names(subset_column)) {

    ids <- meta_data[[column]] %>%
      unfactor_data() %in%
      subset_column[[column]]
    meta_data_sub <- meta_data[ids, ]
    counts_data_sub <- subset(counts_data, , c(TRUE, ids))

  }
  return(list("Counts data" = counts_data_sub, "Meta data" = meta_data_sub))

}
