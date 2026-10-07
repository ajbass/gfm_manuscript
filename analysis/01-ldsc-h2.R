##################################################
### File runs S-LDSC per trait and model ###
## Heritability and enrichment per arm; nothing
## downstream reads it. Completed combinations
## are skipped
##################################################
library(tidyverse)

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) or ANDY (fit stores) to rebuild somewhere else.
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")
ANDY <- Sys.getenv("ANDY", "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
# inputs, read only
SUMSTATS <- Sys.getenv("SUMSTATS",
                       "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
POLYFUN <- Sys.getenv("POLYFUN", "/rds/user/ab3105/hpc-work/polyfun")
FM <- file.path(DATA_GWASAI, "foundational_models")

# The arms, the paper's unless ARMS=all; joint arms write under ldsc/<model>
# and subset arms under ldsc/subset/<model>
ARMS_TSV <- Sys.getenv("ARMS_TSV", "annotations/arms.tsv")
ARM_SET <- Sys.getenv("ARMS", "paper")
arms <- read_tsv(ARMS_TSV, show_col_types = FALSE)
if (ARM_SET == "paper") {
  arms <- arms %>% filter(paper == 1)
}
models <- setNames(file.path(FM, arms$model, arms$dir,
                             paste0(arms$model, ".")),
                   if_else(arms$variant == "joint", arms$model,
                           file.path("subset", arms$model)))

filename <- list.files(file.path(SUMSTATS, "sumstats/UKBB-indep/"),
                       pattern = "parquet", full.names = TRUE)

ldsc_path <- file.path(ANDY, "ldsc/")
weights_path <- file.path(DATA_GWASAI,
                          "baselineLD2.3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.")

args <- commandArgs(trailingOnly = TRUE)
id2 <- as.numeric(args[1])
print(paste("Array ID:", id2))
print(paste("Total files before split:", length(filename)))
id <- (seq_along(filename) %% 25) == id2
filename <- filename[id]
print(paste("Files to process in this job:", length(filename)))

# --ref-ld-chr needs LD scores beside .l2.M
check_model <- function(prefix, chrs = 1:22) {
  miss_ld <- c()
  miss_m <- c()
  for (chr in chrs) {
    ld_ok <- any(file.exists(paste0(prefix, chr, c(".l2.ldscore.gz",
                                                   ".l2.ldscore.parquet",
                                                   ".l2.ldscore"))))
    if (!ld_ok) miss_ld <- c(miss_ld, chr)
    if (!file.exists(paste0(prefix, chr, ".l2.M"))) miss_m <- c(miss_m, chr)
  }
  return(list(miss_ld = miss_ld, miss_m = miss_m))
}

cat("\nLD scores per model\n")
bad <- FALSE
for (model_name in names(models)) {
  chk <- check_model(models[[model_name]])
  n_ok <- 22 - length(chk$miss_ld)
  cat(sprintf("%-30s ldscore %2d/22, .l2.M %2d/22\n",
              model_name, n_ok, 22 - length(chk$miss_m)))
  if (length(chk$miss_ld)) {
    cat("    missing LD scores for chr: ", paste(chk$miss_ld, collapse = ","),
        "\n")
    bad <- TRUE
  }
  if (length(chk$miss_m)) {
    cat("    missing .l2.M for chr:     ", paste(chk$miss_m, collapse = ","),
        "\n")
  }
}
if (bad) {
  stop("LD scores missing: run annotations/02-compute-ldscores.sh first")
}

# S-LDSC per trait and model
for (file in filename) {
  split <- str_remove(basename(file), "\\.sumstats\\.parquet$")

  for (model_name in names(models)) {
    model_path <- models[[model_name]]

    model_out_dir <- file.path(ldsc_path, model_name)
    if (!dir.exists(model_out_dir)) dir.create(model_out_dir, recursive = TRUE)

    out_prefix <- file.path(model_out_dir, split)

    # resume: skip a completed trait/model
    if (file.exists(paste0(out_prefix, ".results"))) {
      cat(paste0("skip (already done): ", split, " | ", model_name, "\n"))
      next
    }

    cat(paste0("\n=== Running trait: ", split, " | Model: ", model_name,
               " ===\n"))

    cmd <- paste(
      paste0("python3 ", POLYFUN, "/ldsc.py"),
      "--h2", file,
      "--ref-ld-chr", model_path,
      "--w-ld-chr", weights_path,
      "--out", out_prefix,
      "--overlap-annot",
      "--not-M-5-50", # M over all SNPs, the MAF > 0.01 convention
      "--print-coefficients"
    )
    system(cmd)
  }
}

cat("\ndone\n")
