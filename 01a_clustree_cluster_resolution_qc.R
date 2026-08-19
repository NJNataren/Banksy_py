#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(clustree)
  library(stringr)
  library(jsonlite)
  library(ggrepel)
})

opt <- function(name, type = "character", default = NULL, help) {
  make_option(
    paste0("--", name),
    type = type,
    default = default,
    help = help
  )
}

option_list <- list(
  opt(
    "cluster_csv",
    help = "CSV from script 00 or script 06 with one row per cell and one BANKSY cluster-label column per resolution."
  ),
  opt("dataset_name", help = "Dataset/sample name used in output filenames."),
  opt(
    "cluster_prefix",
    help = "Shared prefix for resolution columns, e.g. labels_scaled_gaussian_pc30_nc0.20_r"
  ),
  opt(
    "cluster_suffix",
    default = "",
    help = paste(
      "Optional suffix after the numeric resolution in cluster columns.",
      "Use for script 06 recluster outputs, e.g.",
      "_recluster_filtered_qc_v1_qc_pass_only_smoke [default %default]."
    )
  ),
  opt("output_dir", help = "Directory for clustree QC outputs."),
  opt(
    "qc_config",
    default = NULL,
    help = "Optional script 01 QC JSON config containing cluster_col and new_labels for annotated clustree output."
  ),
  opt("width", type = "double", default = 12, help = "Plot width in inches [default %default]."),
  opt("height", type = "double", default = 8, help = "Plot height in inches [default %default]."),
  opt("dpi", type = "integer", default = 300, help = "PNG resolution in dots per inch [default %default].")
)

get_options <- function() {
  if (interactive()) {
    message("Interactive mode detected; using CK_skin_res local fallback values.")
    return(list(
      cluster_csv = "data/xenium/processed/vbct/CK_skin_res/CK_skin_res_cell_cluster_id_across_clustering_res_0.70_0.80_0.90_1.00.csv",
      dataset_name = "CK_skin_res",
      cluster_prefix = "labels_scaled_gaussian_pc30_nc0.20_r",
      cluster_suffix = "",
      output_dir = "data/xenium/output/vbct/CK_skin_res/clustree_qc",
      qc_config = "config/01_QC/vbct/CK_skin_res.json",
      width = 12,
      height = 8,
      dpi = 300
    ))
  }

  parse_args(OptionParser(option_list = option_list))
}

stop_if_missing <- function(value, option_name) {
  if (is.null(value) || is.na(value) || identical(value, "")) {
    stop("Missing required option: ", option_name, call. = FALSE)
  }
}

# Clustree expects columns named as <prefix><numeric_resolution>. Script 06
# recluster tables append a run-specific suffix after the resolution, so strip
# only the prefix/suffix and keep the middle as the resolution label.
extract_resolution_label <- function(cluster_col, cluster_prefix, cluster_suffix = "") {
  resolution_label <- str_remove(cluster_col, fixed(cluster_prefix))
  if (!is.null(cluster_suffix) && !is.na(cluster_suffix) && !identical(cluster_suffix, "")) {
    resolution_label <- substr(
      resolution_label,
      start = 1,
      stop = nchar(resolution_label) - nchar(cluster_suffix)
    )
  }
  resolution_label
}

make_cluster_column_map <- function(cluster_data, cluster_prefix, cluster_suffix = "") {
  # The optional suffix lets the same script handle original script 00 cluster
  # CSVs and script 06 recluster CSVs without changing the clustree prefix.
  candidate_cols <- names(cluster_data)[startsWith(names(cluster_data), cluster_prefix)]
  if (!is.null(cluster_suffix) && !is.na(cluster_suffix) && !identical(cluster_suffix, "")) {
    candidate_cols <- candidate_cols[endsWith(candidate_cols, cluster_suffix)]
  }

  if (length(candidate_cols) == 0) {
    suffix_message <- ifelse(
      is.null(cluster_suffix) || is.na(cluster_suffix) || identical(cluster_suffix, ""),
      "",
      paste0("' and suffix '", cluster_suffix)
    )
    stop(
      "No cluster columns matched prefix '", cluster_prefix, suffix_message,
      "'. Available columns: ",
      paste(names(cluster_data), collapse = ", "),
      call. = FALSE
    )
  }

  column_map <- tibble(
    original_column = candidate_cols,
    resolution = vapply(
      candidate_cols,
      extract_resolution_label,
      character(1),
      cluster_prefix = cluster_prefix,
      cluster_suffix = cluster_suffix
    ),
    resolution_value = suppressWarnings(as.numeric(resolution)),
    clustree_column = paste0(cluster_prefix, resolution)
  )

  if (any(is.na(column_map$resolution_value))) {
    warning(
      "Could not parse all resolution labels as numeric; using lexical column order.",
      call. = FALSE
    )
    column_map <- column_map %>% arrange(resolution)
  } else {
    column_map <- column_map %>% arrange(resolution_value)
  }

  if (any(duplicated(column_map$clustree_column))) {
    duplicated_cols <- column_map$original_column[duplicated(column_map$clustree_column)]
    stop(
      "Multiple input columns map to the same clustree column: ",
      paste(duplicated_cols, collapse = ", "),
      call. = FALSE
    )
  }

  column_map
}

prepare_clustree_data <- function(cluster_data, column_map, cluster_prefix) {
  # Copy selected cluster labels into temporary normalized columns. These are
  # the only BANKSY-prefixed columns that clustree should see.
  plot_data <- cluster_data
  for (idx in seq_len(nrow(column_map))) {
    original_col <- column_map$original_column[[idx]]
    clustree_col <- column_map$clustree_column[[idx]]
    plot_data[[clustree_col]] <- as.factor(cluster_data[[original_col]])
  }

  # Remove original suffixed columns after copying them. If they remain,
  # clustree tries to parse both original and temporary columns as resolutions.
  prefix_cols <- names(plot_data)[startsWith(names(plot_data), cluster_prefix)]
  extra_cluster_cols <- setdiff(prefix_cols, column_map$clustree_column)
  if (length(extra_cluster_cols) > 0) {
    plot_data <- plot_data %>% select(-all_of(extra_cluster_cols))
  }

  plot_data
}

sort_cluster_columns <- function(cluster_cols, cluster_prefix, cluster_suffix = "") {
  resolution_labels <- vapply(
    cluster_cols,
    extract_resolution_label,
    character(1),
    cluster_prefix = cluster_prefix,
    cluster_suffix = cluster_suffix
  )
  resolution_values <- suppressWarnings(as.numeric(resolution_labels))

  if (all(!is.na(resolution_values))) {
    cluster_cols[order(resolution_values)]
  } else {
    warning(
      "Could not parse all resolution labels as numeric; using lexical column order.",
      call. = FALSE
    )
    sort(cluster_cols)
  }
}

summarise_cluster_columns <- function(cluster_data, column_map) {
  bind_rows(lapply(seq_len(nrow(column_map)), function(idx) {
    original_col <- column_map$original_column[[idx]]
    clustree_col <- column_map$clustree_column[[idx]]
    values <- cluster_data[[clustree_col]]

    tibble(
      original_column = original_col,
      clustree_column = clustree_col,
      resolution = column_map$resolution[[idx]],
      resolution_value = column_map$resolution_value[[idx]],
      n_cells = sum(!is.na(values)),
      n_clusters = n_distinct(values, na.rm = TRUE)
    )
  }))
}

read_annotation_config <- function(qc_config, cluster_data, cluster_prefix, column_map) {
  # Annotation configs may refer to either the original CSV column or the
  # normalized clustree column; map both forms back to the same resolution.
  if (is.null(qc_config) || is.na(qc_config) || identical(qc_config, "")) {
    return(NULL)
  }

  if (!file.exists(qc_config)) {
    stop("QC config does not exist: ", qc_config, call. = FALSE)
  }

  config <- fromJSON(qc_config)
  if (is.null(config$cluster_col) || is.null(config$new_labels)) {
    stop("QC config must contain cluster_col and new_labels: ", qc_config, call. = FALSE)
  }

  annotation_row <- column_map %>%
    filter(
      original_column == config$cluster_col |
        clustree_column == config$cluster_col
    )

  if (nrow(annotation_row) != 1) {
    stop(
      "Configured annotation cluster_col is not present in selected cluster columns: ",
      config$cluster_col,
      call. = FALSE
    )
  }

  annotation_col <- annotation_row$clustree_column[[1]]
  resolution_label <- annotation_row$resolution[[1]]
  resolution_value <- annotation_row$resolution_value[[1]]
  if (is.na(resolution_value)) {
    stop("Could not parse annotation resolution from cluster_col: ", config$cluster_col, call. = FALSE)
  }

  label_values <- unlist(config$new_labels, use.names = TRUE)
  cluster_counts <- cluster_data %>%
    count(cluster_id = as.character(.data[[annotation_col]]), name = "n_cells")

  annotation_data <- tibble(
    cluster_col = config$cluster_col,
    clustree_col = annotation_col,
    resolution = resolution_label,
    resolution_value = resolution_value,
    cluster_id = names(label_values),
    annotation_label = as.character(label_values)
  ) %>%
    left_join(cluster_counts, by = "cluster_id") %>%
    arrange(suppressWarnings(as.numeric(cluster_id)))

  list(
    qc_config = qc_config,
    cluster_col = config$cluster_col,
    clustree_col = annotation_col,
    resolution_label = resolution_label,
    resolution_value = resolution_value,
    annotation_data = annotation_data
  )
}

make_annotated_plot <- function(base_plot, plot_data, cluster_prefix, annotation_info) {
  tree_layout <- as.data.frame(clustree(plot_data, prefix = cluster_prefix, return = "layout"))

  annotation_layout <- tree_layout %>%
    mutate(
      resolution_value = suppressWarnings(as.numeric(as.character(.data[[cluster_prefix]]))),
      cluster_id = as.character(cluster)
    ) %>%
    filter(abs(resolution_value - annotation_info$resolution_value) < 1e-8) %>%
    left_join(annotation_info$annotation_data, by = c("resolution_value", "cluster_id")) %>%
    filter(!is.na(annotation_label)) %>%
    mutate(annotation_plot_label = str_replace_all(annotation_label, "_", " "))

  if (nrow(annotation_layout) == 0) {
    warning(
      "No clustree nodes matched annotation labels for ",
      annotation_info$cluster_col,
      call. = FALSE
    )
    return(NULL)
  }

  base_plot +
    geom_label_repel(
      data = annotation_layout,
      aes(x = x, y = y, label = annotation_plot_label),
      inherit.aes = FALSE,
      size = 2.3,
      linewidth = 0.2,
      label.padding = grid::unit(0.12, "lines"),
      box.padding = grid::unit(0.45, "lines"),
      point.padding = grid::unit(0.45, "lines"),
      min.segment.length = 0,
      segment.alpha = 0.65,
      segment.size = 0.25,
      force = 2,
      force_pull = 0.6,
      max.overlaps = Inf,
      direction = "both",
      seed = 123,
      alpha = 0.94
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 45, 10, 45)) +
    labs(subtitle = paste(
      "Cluster prefix:", cluster_prefix,
      "| annotations:", annotation_info$cluster_col
    ))
}

opts <- get_options()

stop_if_missing(opts$cluster_csv, "--cluster_csv")
stop_if_missing(opts$dataset_name, "--dataset_name")
stop_if_missing(opts$cluster_prefix, "--cluster_prefix")
stop_if_missing(opts$output_dir, "--output_dir")
if (is.null(opts$cluster_suffix) || is.na(opts$cluster_suffix)) {
  opts$cluster_suffix <- ""
}

if (!file.exists(opts$cluster_csv)) {
  stop("Cluster CSV does not exist: ", opts$cluster_csv, call. = FALSE)
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading cluster assignments: ", opts$cluster_csv)
cluster_data <- read_csv(opts$cluster_csv, show_col_types = FALSE)

# Build a small lookup table from input columns to the normalized columns used
# for plotting and summaries. With an empty suffix this is effectively a no-op.
column_map <- make_cluster_column_map(cluster_data, opts$cluster_prefix, opts$cluster_suffix)
cluster_cols <- column_map$original_column

message("Detected cluster columns:")
message(paste("  -", cluster_cols, collapse = "\n"))
if (!identical(opts$cluster_suffix, "")) {
  message("Using cluster suffix: ", opts$cluster_suffix)
  message("Temporary clustree columns:")
  message(paste("  -", column_map$clustree_column, collapse = "\n"))
}

plot_data <- prepare_clustree_data(cluster_data, column_map, opts$cluster_prefix)

# Write the column map alongside basic counts so later review can trace each
# plotted resolution back to its original script 00/script 06 CSV column.
summary_data <- summarise_cluster_columns(plot_data, column_map)
summary_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_input_columns.csv"))
write_csv(summary_data, summary_path)

annotation_info <- read_annotation_config(opts$qc_config, plot_data, opts$cluster_prefix, column_map)
annotation_path <- NULL
if (!is.null(annotation_info)) {
  annotation_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_annotation_labels.csv"))
  write_csv(annotation_info$annotation_data, annotation_path)
}

tree_plot <- clustree(plot_data, prefix = opts$cluster_prefix) +
  labs(
    title = paste(opts$dataset_name, "BANKSY Clustree Resolution QC"),
    subtitle = paste("Cluster prefix:", opts$cluster_prefix, "| cluster suffix:", opts$cluster_suffix)
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

sc3_stability_plot <- clustree(
  plot_data,
  prefix = opts$cluster_prefix,
  node_colour = "sc3_stability"
) +
  labs(
    title = paste(opts$dataset_name, "BANKSY Clustree SC3 Stability QC"),
    subtitle = paste("Cluster prefix:", opts$cluster_prefix, "| cluster suffix:", opts$cluster_suffix),
    colour = "SC3 stability"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

annotated_plot <- NULL
if (!is.null(annotation_info)) {
  annotated_plot <- make_annotated_plot(tree_plot, plot_data, opts$cluster_prefix, annotation_info)
}

png_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_resolution_qc.png"))
pdf_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_resolution_qc.pdf"))
sc3_png_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_resolution_qc_sc3_stability.png"))
sc3_pdf_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_resolution_qc_sc3_stability.pdf"))
annotated_png_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_resolution_qc_annotated.png"))
annotated_pdf_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_resolution_qc_annotated.pdf"))

ggsave(png_path, tree_plot, width = opts$width, height = opts$height, dpi = opts$dpi)
ggsave(pdf_path, tree_plot, width = opts$width, height = opts$height)
ggsave(sc3_png_path, sc3_stability_plot, width = opts$width, height = opts$height, dpi = opts$dpi)
ggsave(sc3_pdf_path, sc3_stability_plot, width = opts$width, height = opts$height)
if (!is.null(annotated_plot)) {
  ggsave(annotated_png_path, annotated_plot, width = opts$width, height = opts$height, dpi = opts$dpi)
  ggsave(annotated_pdf_path, annotated_plot, width = opts$width, height = opts$height)
}

message("Wrote clustree PNG: ", png_path)
message("Wrote clustree PDF: ", pdf_path)
message("Wrote SC3 stability clustree PNG: ", sc3_png_path)
message("Wrote SC3 stability clustree PDF: ", sc3_pdf_path)
if (!is.null(annotated_plot)) {
  message("Wrote annotated clustree PNG: ", annotated_png_path)
  message("Wrote annotated clustree PDF: ", annotated_pdf_path)
}
message("Wrote input summary: ", summary_path)
if (!is.null(annotation_path)) {
  message("Wrote annotation labels: ", annotation_path)
}
