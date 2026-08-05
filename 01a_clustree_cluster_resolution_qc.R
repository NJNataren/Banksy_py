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

option_list <- list(
  make_option(
    "--cluster_csv",
    type = "character",
    help = "CSV from script 00 with one row per cell and one BANKSY cluster-label column per resolution."
  ),
  make_option(
    "--dataset_name",
    type = "character",
    help = "Dataset/sample name used in output filenames."
  ),
  make_option(
    "--cluster_prefix",
    type = "character",
    help = "Shared prefix for resolution columns, e.g. labels_scaled_gaussian_pc30_nc0.20_r"
  ),
  make_option(
    "--output_dir",
    type = "character",
    help = "Directory for clustree QC outputs."
  ),
  make_option(
    "--qc_config",
    type = "character",
    default = NULL,
    help = "Optional script 01 QC JSON config containing cluster_col and new_labels for annotated clustree output."
  ),
  make_option(
    "--width",
    type = "double",
    default = 12,
    help = "Plot width in inches [default %default]."
  ),
  make_option(
    "--height",
    type = "double",
    default = 8,
    help = "Plot height in inches [default %default]."
  ),
  make_option(
    "--dpi",
    type = "integer",
    default = 300,
    help = "PNG resolution in dots per inch [default %default]."
  )
)

get_options <- function() {
  if (interactive()) {
    message("Interactive mode detected; using CK_skin_res local fallback values.")
    return(list(
      cluster_csv = "data/xenium/processed/vbct/CK_skin_res/CK_skin_res_cell_cluster_id_across_clustering_res_0.70_0.80_0.90_1.00.csv",
      dataset_name = "CK_skin_res",
      cluster_prefix = "labels_scaled_gaussian_pc30_nc0.20_r",
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

summarise_cluster_columns <- function(cluster_data, cluster_cols, cluster_prefix) {
  bind_rows(lapply(cluster_cols, function(col_name) {
    resolution <- str_remove(col_name, fixed(cluster_prefix))
    values <- cluster_data[[col_name]]

    tibble(
      column = col_name,
      resolution = resolution,
      n_cells = sum(!is.na(values)),
      n_clusters = n_distinct(values, na.rm = TRUE)
    )
  }))
}

read_annotation_config <- function(qc_config, cluster_data, cluster_prefix) {
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

  if (!config$cluster_col %in% names(cluster_data)) {
    stop(
      "Configured annotation cluster_col is not present in cluster CSV: ",
      config$cluster_col,
      call. = FALSE
    )
  }

  resolution_label <- str_remove(config$cluster_col, fixed(cluster_prefix))
  resolution_value <- suppressWarnings(as.numeric(resolution_label))
  if (is.na(resolution_value)) {
    stop("Could not parse annotation resolution from cluster_col: ", config$cluster_col, call. = FALSE)
  }

  label_values <- unlist(config$new_labels, use.names = TRUE)
  cluster_counts <- cluster_data %>%
    count(cluster_id = as.character(.data[[config$cluster_col]]), name = "n_cells")

  annotation_data <- tibble(
    cluster_col = config$cluster_col,
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

if (!file.exists(opts$cluster_csv)) {
  stop("Cluster CSV does not exist: ", opts$cluster_csv, call. = FALSE)
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading cluster assignments: ", opts$cluster_csv)
cluster_data <- read_csv(opts$cluster_csv, show_col_types = FALSE)

cluster_cols <- names(cluster_data)[startsWith(names(cluster_data), opts$cluster_prefix)]
cluster_cols <- sort_cluster_columns(cluster_cols, opts$cluster_prefix)

if (length(cluster_cols) == 0) {
  stop(
    "No cluster columns matched prefix '", opts$cluster_prefix, "'. Available columns: ",
    paste(names(cluster_data), collapse = ", "),
    call. = FALSE
  )
}

message("Detected cluster columns:")
message(paste("  -", cluster_cols, collapse = "\n"))

plot_data <- cluster_data %>%
  mutate(across(all_of(cluster_cols), as.factor))

summary_data <- summarise_cluster_columns(plot_data, cluster_cols, opts$cluster_prefix)
summary_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_input_columns.csv"))
write_csv(summary_data, summary_path)

annotation_info <- read_annotation_config(opts$qc_config, plot_data, opts$cluster_prefix)
annotation_path <- NULL
if (!is.null(annotation_info)) {
  annotation_path <- file.path(opts$output_dir, paste0(opts$dataset_name, "_clustree_annotation_labels.csv"))
  write_csv(annotation_info$annotation_data, annotation_path)
}

tree_plot <- clustree(plot_data, prefix = opts$cluster_prefix) +
  labs(
    title = paste(opts$dataset_name, "BANKSY Clustree Resolution QC"),
    subtitle = paste("Cluster prefix:", opts$cluster_prefix)
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
    subtitle = paste("Cluster prefix:", opts$cluster_prefix),
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
