# Setup a local R environment for Xenium BANKSY clustree QC.

required_packages <- c(
  "renv",
  "clustree",
  "optparse",
  "readr",
  "dplyr",
  "ggplot2",
  "stringr",
  "jsonlite",
  "ggrepel"
)

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

if (file.exists("renv.lock")) {
  message("Found renv.lock; restoring recorded R package versions.")
  renv::restore(prompt = FALSE)
} else {
  message("No renv.lock found; initializing renv and installing clustree QC packages.")
  renv::init(bare = TRUE)
  renv::install(setdiff(required_packages, "renv"))
  renv::snapshot(prompt = FALSE)
}

message("Clustree R environment is ready.")
