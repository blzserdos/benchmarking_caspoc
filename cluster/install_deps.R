#!/usr/bin/env Rscript
# =============================================================================
# install_deps.R
# Run ONCE on the Saga login node.
#
# Usage:
#   module load R/4.5.2-gfbf-2025b       # match the version used in submit.sh
#   Rscript cluster/install_deps.R
# =============================================================================

# --- Personal library ---
user_lib <- Sys.getenv("R_LIBS_USER", unset = "~/R/library")
user_lib <- path.expand(user_lib)
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
message("Installing into: ", user_lib)

repos <- c(CRAN = "https://cloud.r-project.org")

# --- CRAN packages ---
cran_pkgs <- c(
  "BiocManager", "devtools",
  "caret", "dplyr", "tibble", "MASS",
  "future", "furrr", "parallelly"
)
to_install <- setdiff(cran_pkgs, rownames(installed.packages(lib.loc = user_lib)))
if (length(to_install)) {
  install.packages(to_install, lib = user_lib, repos = repos)
}

# --- Bioconductor: mixOmics ---
if (!"mixOmics" %in% rownames(installed.packages(lib.loc = user_lib))) {
  BiocManager::install("mixOmics", lib = user_lib, update = FALSE, ask = FALSE)
}

# --- CASPOC (pinned version) ---
CASPOC_REF <- "CASPOC-v1.1"
devtools::install_github(
  "jonathanth/caspoc",
  ref     = CASPOC_REF,
  lib     = user_lib,
  upgrade = "never"
)

# --- Verify ---
required <- c("mixOmics", "caret", "caspoc", "dplyr", "tibble", "MASS",
              "future", "furrr", "parallelly")
missing  <- required[!required %in% rownames(installed.packages(lib.loc = user_lib))]
if (length(missing)) {
  stop("Missing packages after install: ", paste(missing, collapse = ", "))
}

caspoc_ver <- as.character(packageVersion("caspoc", lib.loc = user_lib))
message(sprintf("All dependencies installed. caspoc version: %s (ref: %s)",
                caspoc_ver, CASPOC_REF))
