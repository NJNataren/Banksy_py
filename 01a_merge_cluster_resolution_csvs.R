#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
})

option_list <- list(
  make_option(
    "--input_dir",
    type = "character",
    help = "Directory containing split script 00 cluster assignment CSV files."
  ),
  make_option(
    "--dataset_name",
    type = "character",
    help = "Dataset/sample name used to select files and name the merged output."
  ),
  make_option(
    "--cluster_prefix",
    type = "character",
    help = "Shared prefix for resolution columns, e.g. labels_scaled_gaussian_pc30_nc0.20_r"
  ),
  make_option(
    "--output_csv",
    type = "character",
    default = NULL,
    help = "Merged output CSV path. If omitted, a resolution-labelled filename is written in input_dir."
  ),
  make_option(
    "--file_pattern",
    type = "character",
    default = NULL,
    help = "Optional regex for candidate CSV basenames. Defaults to dataset clustree CSV files."
  ),
  make_option(
    "--cell_id_col",
    type = "character",
    default = "index",
    help = "Cell ID column used to join split files [default %default]."
  )
)

stop_if_missing <- function(value, option_name) {
  if (is.null(value) || is.na(value) || identical(value, "")) {
    stop("Missing required option: ", option_name, call. = FALSE)
  }
}

drop_csv_index_columns <- function(data) {
  index_like <- names(data) %in% c("", "...1", "X", "X1")
  if (any(index_like)) {
    data <- data[, !index_like, drop = FALSE]
  }
  data
}

read_cluster_file <- function(path, cell_id_col, cluster_prefix) {
  data <- read_csv(path, show_col_types = FALSE, name_repair = "unique_quiet")
  data <- drop_csv_index_columns(data)

  cluster_cols <- names(data)[startsWith(names(data), cluster_prefix)]
  if (length(cluster_cols) == 0) {
    return(NULL)
  }
  if (length(cluster_cols) > 1) {
    message("Skipping already-merged or multi-resolution CSV: ", path)
    return(NULL)
  }

  if (!cell_id_col %in% names(data)) {
    stop("Missing cell ID column '", cell_id_col, "' in ", path, call. = FALSE)
  }

  if (anyDuplicated(data[[cell_id_col]]) > 0) {
    stop("Duplicate cell IDs in ", path, call. = FALSE)
  }

  data %>%
    select(all_of(c(cell_id_col, cluster_cols))) %>%
    mutate(across(all_of(cluster_cols), as.character))
}

sort_cluster_columns <- function(cluster_cols, cluster_prefix) {
  resolution_labels <- str_remove(cluster_cols, fixed(cluster_prefix))
  resolution_values <- suppressWarnings(as.numeric(resolution_labels))

  if (all(!is.na(resolution_values))) {
    cluster_cols[order(resolution_values)]
  } else {
    warning(
      "Could not parse all resolution suffixes as numeric; using lexical column order.",
      call. = FALSE
    )
    sort(cluster_cols)
  }
}

merge_cluster_tables <- function(tables, paths, cell_id_col, cluster_prefix) {
  merged <- tables[[1]]
  reference_cells <- merged[[cell_id_col]]

  if (length(tables) > 1) {
    for (i in seq(2, length(tables))) {
      current <- tables[[i]]
      current_cells <- current[[cell_id_col]]

      if (!setequal(reference_cells, current_cells)) {
        missing_from_current <- setdiff(reference_cells, current_cells)
        extra_in_current <- setdiff(current_cells, reference_cells)
        stop(
          "Cell ID sets differ in ", paths[[i]], ". Missing from current: ",
          length(missing_from_current), "; extra in current: ", length(extra_in_current),
          call. = FALSE
        )
      }

      current <- current[match(reference_cells, current[[cell_id_col]]), , drop = FALSE]
      overlap_cols <- intersect(
        names(merged)[startsWith(names(merged), cluster_prefix)],
        names(current)[startsWith(names(current), cluster_prefix)]
      )

      for (col_name in overlap_cols) {
        if (!identical(merged[[col_name]], current[[col_name]])) {
          stop("Conflicting assignments for duplicate cluster column: ", col_name, call. = FALSE)
        }
      }

      current_new_cols <- setdiff(names(current), c(cell_id_col, overlap_cols))
      merged <- bind_cols(merged, current[, current_new_cols, drop = FALSE])
    }
  }

  cluster_cols <- sort_cluster_columns(names(merged)[startsWith(names(merged), cluster_prefix)], cluster_prefix)
  merged %>% select(all_of(c(cell_id_col, cluster_cols)))
}

make_default_output_csv <- function(input_dir, dataset_name, cluster_cols, cluster_prefix) {
  resolution_labels <- str_remove(cluster_cols, fixed(cluster_prefix))
  resolution_suffix <- paste(resolution_labels, collapse = "_")
  file.path(
    input_dir,
    paste0(dataset_name, "_cell_cluster_id_across_clustering_res_", resolution_suffix, ".csv")
  )
}

opts <- parse_args(OptionParser(option_list = option_list))

stop_if_missing(opts$input_dir, "--input_dir")
stop_if_missing(opts$dataset_name, "--dataset_name")
stop_if_missing(opts$cluster_prefix, "--cluster_prefix")
stop_if_missing(opts$cell_id_col, "--cell_id_col")

if (!dir.exists(opts$input_dir)) {
  stop("Input directory does not exist: ", opts$input_dir, call. = FALSE)
}

file_pattern <- opts$file_pattern
if (is.null(file_pattern) || is.na(file_pattern) || identical(file_pattern, "")) {
  file_pattern <- paste0(
    "^",
    opts$dataset_name,
    "_cell_cluster_id_across_clustering_res_.*[.]csv$"
  )
}

candidate_paths <- list.files(
  opts$input_dir,
  pattern = file_pattern,
  full.names = TRUE
)
candidate_paths <- candidate_paths[!str_detect(basename(candidate_paths), "_merge_summary[.]csv$")]

if (!is.null(opts$output_csv) && !is.na(opts$output_csv) && !identical(opts$output_csv, "")) {
  candidate_paths <- setdiff(normalizePath(candidate_paths, mustWork = FALSE), normalizePath(opts$output_csv, mustWork = FALSE))
}

if (length(candidate_paths) == 0) {
  stop("No candidate CSV files matched pattern '", file_pattern, "' in ", opts$input_dir, call. = FALSE)
}

message("Candidate CSV files:")
message(paste("  -", candidate_paths, collapse = "\n"))

tables <- list()
used_paths <- character()
for (path in candidate_paths) {
  table <- read_cluster_file(path, opts$cell_id_col, opts$cluster_prefix)
  if (!is.null(table)) {
    tables[[length(tables) + 1]] <- table
    used_paths <- c(used_paths, path)
  }
}

if (length(tables) == 0) {
  stop("No cluster columns matched prefix '", opts$cluster_prefix, "' in candidate CSV files.", call. = FALSE)
}

merged <- merge_cluster_tables(tables, used_paths, opts$cell_id_col, opts$cluster_prefix)
cluster_cols <- names(merged)[startsWith(names(merged), opts$cluster_prefix)]

output_csv <- opts$output_csv
if (is.null(output_csv) || is.na(output_csv) || identical(output_csv, "")) {
  output_csv <- make_default_output_csv(opts$input_dir, opts$dataset_name, cluster_cols, opts$cluster_prefix)
}

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write_csv(merged, output_csv)

summary_data <- tibble(
  column = cluster_cols,
  resolution = str_remove(cluster_cols, fixed(opts$cluster_prefix)),
  n_cells = nrow(merged),
  n_clusters = vapply(cluster_cols, function(col_name) n_distinct(merged[[col_name]], na.rm = TRUE), integer(1))
)
summary_csv <- str_replace(output_csv, "[.]csv$", "_merge_summary.csv")
write_csv(summary_data, summary_csv)

message("Merged ", length(used_paths), " CSV files.")
message("Detected cluster columns:")
message(paste("  -", cluster_cols, collapse = "\n"))
message("Wrote merged CSV: ", output_csv)
message("Wrote merge summary: ", summary_csv)
