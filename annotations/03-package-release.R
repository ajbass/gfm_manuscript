##################################################
### File packages the release for distribution ###
## One tar per model of the release files, the
## column guide and SNP list, the LD-pruned SNPs,
## the AlphaGenome API export trimmed to the
## columns 00 reads, and a manifest with md5 sums
##################################################
library(tidyverse)

REL_MD <- normalizePath("annotations/RELEASE.md", mustWork = TRUE)
PRUNED <- normalizePath("qc_panukk/pruned_snps_0.3.tsv", mustWork = FALSE)

# Cluster paths this analysis was run on.
setwd(Sys.getenv("DATA_GWASAI", "/home/ab3105/rds/hpc-work/data_gwasAI"))
REL <- "./foundational_models/release"
AG_SCORES <- "./foundational_models/alphagenome/scores.csv"
OUT <- Sys.getenv("RELEASE_OUT", "./foundational_models/release_bundle")

AG_TYPES <- c("ATAC", "CAGE", "CHIP_HISTONE", "CHIP_TF", "CONTACT_MAPS",
              "DNASE", "PROCAP", "RNA_SEQ")

AG_NEED <- c("variant_id", paste0(AG_TYPES, "_raw"))
AG_KEEP <- c(AG_NEED, "chrom", "POS", "REF", "ALT")
AG_SPEC <- do.call(cols, c(list(.default = col_character()),
                           setNames(rep(list(col_double()), length(AG_TYPES)),
                                    paste0(AG_TYPES, "_raw"))))

guide <- read_tsv(file.path(REL, "columns.tsv"), show_col_types = FALSE)
models <- unique(guide$file)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## PART 1: one tar per model, so a reader takes only the models they need
for (m in models) {
  f <- file.path(OUT, paste0(m, ".tar"))
  if (file.exists(f)) {
    message(m, ".tar exists")
    next
  }
  files <- sprintf("%s.%d.parquet", m, 1:22)
  miss <- files[!file.exists(file.path(REL, files))]
  if (length(miss)) {
    stop(m, ": release files missing: ", paste(miss, collapse = ", "))
  }
  st <- system2("tar", c("-cf", shQuote(f), "-C", shQuote(REL), files))
  if (st != 0) {
    stop("tar failed for ", m)
  }
  message(sprintf("  %-12s %.2f GB", paste0(m, ".tar"), file.size(f) / 1e9))
}
invisible(file.copy(file.path(REL, c("columns.tsv", "snps.tsv")), OUT,
                    overwrite = TRUE))
invisible(file.copy(REL_MD, OUT, overwrite = TRUE))
if (file.exists(PRUNED)) {
  invisible(file.copy(PRUNED, OUT, overwrite = TRUE))
} else {
  warning("no qc_panukk/pruned_snps_0.3.tsv in the run tree; the bundle ",
          "lacks the LD-pruned SNP list (data/04-pruned-snps.R writes it)")
}

## PART 2: the AlphaGenome export, the one input of 00 with no public source
f_ag <- file.path(OUT, "alphagenome_scores.csv.gz")
if (file.exists(f_ag)) {
  message("alphagenome_scores.csv.gz exists")
} else {
  if (!file.exists(AG_SCORES)) {
    AG_SCORES <- paste0(AG_SCORES, ".gz")
  }
  ag <- read_csv(AG_SCORES, col_types = AG_SPEC, name_repair = "unique_quiet")
  if (!all(AG_NEED %in% names(ag))) {
    stop("scores.csv lacks ",
         paste(setdiff(AG_NEED, names(ag)), collapse = ", "))
  }
  ag <- ag %>% select(any_of(AG_KEEP))
  write_csv(ag, f_ag)
  message(sprintf("  alphagenome_scores.csv.gz: %d variants, %d columns",
                  nrow(ag), ncol(ag)))
  rm(ag)
}

## PART 3: the manifest
fs <- list.files(OUT, full.names = TRUE)
fs <- fs[basename(fs) != "MANIFEST.tsv"]
manifest <- tibble(file = basename(fs), bytes = file.size(fs),
                   md5 = unname(tools::md5sum(fs)))
write_tsv(manifest, file.path(OUT, "MANIFEST.tsv"))
print(manifest, n = Inf)
cat("\nbundle in", OUT, "\n")
