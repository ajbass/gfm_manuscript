##################################################
### File extracts lead-SNP statistics ###
##################################################
library(tidyverse)

# Each model x trait reads a ~500 MB file, so the work is sharded across
# an array:
#   Rscript 06-lead-stats.R <task_id> <n_tasks>
#   Rscript 06-lead-stats.R merge
# VARIANTS= and MODELS= restrict the work.
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1 && args[1] == "merge") "merge" else "extract"
task_id <- if (length(args) >= 1 && MODE == "extract") as.integer(args[1]) else 0L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
KEEP_VARIANTS <- strsplit(Sys.getenv("VARIANTS", "joint,subset"), ",")[[1]]
KEEP_MODELS <- Sys.getenv("MODELS", "")

# Cluster paths this analysis was run on.
BASE <- Sys.getenv("ANDY", "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
RES_ROOT <- file.path(BASE, "sffdr/final_results")

exact_p <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[!is.na(x) & x < .Machine$double.xmin] <- .Machine$double.xmin
  return(x)
}
p_cols <- function(f, which = c("P", "fP", "fpvalues")) {
  hdr <- names(read_tsv(f, n_max = 0, show_col_types = FALSE))
  want <- intersect(which, hdr)
  spec <- replicate(length(want), col_character(), simplify = FALSE)
  names(spec) <- want
  return(do.call(cols, c(spec, list(.default = col_guess()))))
}

OUT_ROOT <- file.path(normalizePath("."), "indep-discoveries-1000G")
OUT_FILE <- file.path(OUT_ROOT, "sigresults_1000G.tsv")
SHARD_DIR <- file.path(OUT_ROOT, "_leadstats_shards")
dir.create(SHARD_DIR, showWarnings = FALSE, recursive = TRUE)
AH_WINDOW <- 100000

if (MODE == "merge") {
  sh <- list.files(SHARD_DIR, pattern = "\\.tsv$", full.names = TRUE)
  if (!length(sh)) {
    stop("no shards in ", SHARD_DIR)
  }
  # the shards carry P and fP exact (2.2250738585072014e-308); read them
  # back the same way
  out <- map_dfr(sh, ~ suppressMessages(read_tsv(.x, col_types = p_cols(.x),
                                                 show_col_types = FALSE)) %>%
                   mutate(across(any_of(c("P", "fP", "fpvalues")), exact_p)))
  n_raw <- nrow(out)

  # Shards must partition the work; unsharded duplicates N-fold.
  # NA-stats rows sort last, so distinct() keeps the row with stats
  key <- c("variant", "model", "trait", "type", "SNP")
  out <- out %>%
    arrange(across(all_of(key)), is.na(Z)) %>%
    distinct(across(all_of(key)), .keep_all = TRUE)
  n_dup <- n_raw - nrow(out)
  if (n_dup > 0) {
    fold <- round(n_raw / nrow(out), 2)
    message("!! ", n_dup, " duplicate rows removed (", fold, "x duplication).")
    message("   Each shard appears to contain the full dataset rather than a slice.")
    message("   Check that the submitted script is the sharded version:")
    message("     grep -n 'task %d/%d' 06-lead-stats.R")
    message("   and that the submit script passes both <task_id> and <n_tasks>.")
  }

  write_tsv(out, OUT_FILE)
  cat(sprintf("merged %d shards -> %s  (%d leads, %d duplicates removed)\n",
              length(sh), OUT_FILE, nrow(out), n_dup))
  print(out %>% count(variant, model, type) %>%
          pivot_wider(names_from = type, values_from = n), n = 40)
  miss <- out %>% filter(is.na(Z))
  if (nrow(miss)) {
    cat(sprintf("\n!! %d leads without statistics:\n", nrow(miss)))
    print(miss %>% count(variant, model, trait), n = 40)
  }
  quit(save = "no")
}

read_clumps <- function(f) {
  d <- suppressWarnings(suppressMessages(read_table(f, show_col_types = FALSE)))
  if (!nrow(d)) {
    return(NULL)
  }
  nm <- names(d)
  return(tibble(CHR = as.integer(d[[intersect(c("#CHROM", "CHR"), nm)[1]]]),
                POS = as.numeric(d[[intersect(c("POS", "BP"), nm)[1]]]),
                SNP = as.character(d[[intersect(c("ID", "SNP"), nm)[1]]])))
}

# Every lead SNP from stage 1
files <- list.files(OUT_ROOT, pattern = "\\.clumps$", recursive = TRUE,
                    full.names = TRUE)
files <- files[!str_detect(files, "_paired|_union")]
# Ignore clumps of '.full' whole-data fits
files <- files[!str_detect(basename(files), "\\.full-")]
if (!length(files)) {
  stop("no stage-1 .clumps under ", OUT_ROOT)
}
message("stage-1 clump files: ", length(files))

leads <- map_dfr(files, function(f) {
  rel <- str_remove(f, fixed(paste0(OUT_ROOT, "/")))
  p <- str_split(rel, "/")[[1]]
  b <- str_split(str_remove(basename(f), fixed(".clumps")), "-")[[1]]
  d <- read_clumps(f)
  if (is.null(d)) {
    return(NULL)
  }
  d %>% mutate(variant = p[1], model = p[2], type = b[length(b)],
               trait = paste(b[seq_len(length(b) - 2)], collapse = "-"))
})
message(sprintf("leads: %d across %d model-directories", nrow(leads),
                n_distinct(leads$model)))

# Summary statistics for those leads
jobs <- leads %>% distinct(variant, model, trait) %>%
  filter(variant %in% KEEP_VARIANTS) %>%
  { if (nzchar(KEEP_MODELS)) filter(.,
                                    model %in% strsplit(KEEP_MODELS,
                                                        ",")[[1]]) else . } %>%
  arrange(variant, model, trait) %>%
  mutate(k = row_number()) %>%
  filter(n_tasks <= 1 | (k %% n_tasks) == task_id)
message(sprintf("task %d/%d: %d model x trait combinations", task_id, n_tasks,
                nrow(jobs)))
if (!nrow(jobs)) {
  message("nothing to do")
  quit(save = "no")
}

t0 <- Sys.time()
stats <- map_dfr(seq_len(nrow(jobs)), function(i) {
  j <- jobs[i, ]
  message(sprintf("  [%3d/%3d %5.1f min] %s | %s | %s", i, nrow(jobs),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  j$variant, j$model, j$trait))
  f <- if (j$variant == "joint") {
    file.path(RES_ROOT, j$model, paste0(j$trait, ".sffdr_final.gz"))
  } else {
    file.path(RES_ROOT, "subset", j$model, paste0(j$trait, ".sffdr_final.gz"))
  }
  if (!file.exists(f)) {
    message("  missing: ", f)
    return(NULL)
  }
  want <- leads %>%
    filter(variant == j$variant, model == j$model, trait == j$trait) %>%
    pull(SNP)
  d <- tryCatch(suppressMessages(read_tsv(f, show_col_types = FALSE,
         col_select = any_of(c("SNP", "CHR", "BP", "A1", "A2", "EAF", "N",
                               "Beta", "se", "Z",
                               "P", "fP", "fpvalues", "echisq")),
         col_types = p_cols(f))) %>%
    mutate(across(any_of(c("P", "fP", "fpvalues")), exact_p)),
       error = function(e) {
         message("  !! unreadable (corrupt?): ", basename(f))
         NULL
       })
  if (is.null(d)) {
    return(NULL)
  }
  fcol <- intersect(c("fP", "fpvalues"), names(d))[1]
  d %>% filter(SNP %in% want) %>%
    mutate(fP = .data[[fcol]], variant = j$variant, model = j$model,
    trait = j$trait) %>%
    select(variant, model, trait, SNP, CHR, BP,
           any_of(c("A1", "A2", "EAF", "N", "Beta", "se")),
           Z, P, fP, any_of("echisq"))
})

# This task's jobs only: without it the merge sees n_tasks copies of
# every lead
out <- leads %>%
  semi_join(jobs %>% select(variant, model, trait),
            by = c("variant", "model", "trait")) %>%
  select(variant, model, trait, type, SNP) %>%
  left_join(stats, by = c("variant", "model", "trait", "SNP"))
message(sprintf("task %d: %d leads across %d model x trait combinations",
                task_id, nrow(out), nrow(jobs)))

# No Classification column: joined on
# lead_SNP by 07-collapse-classify.R

if (n_tasks > 1L) {
  f <- file.path(SHARD_DIR, sprintf("shard_%03d.tsv", task_id))
  write_tsv(out, f)
  cat("\nwrote shard ", f, "  (", nrow(out), " leads)\n", sep = "")
  cat("when every task is done:  Rscript 06-lead-stats.R merge\n")
  quit(save = "no")
}
write_tsv(out, OUT_FILE)
cat("\nwrote ", OUT_FILE, "  (", nrow(out), " leads)\n", sep = "")

cat("\n=== leads per variant/model/type ===\n")
print(out %>% count(variant, model, type) %>%
        pivot_wider(names_from = type, values_from = n), n = 40)

miss <- out %>% filter(is.na(Z))
if (nrow(miss)) {
  cat(sprintf("\n!! %d leads have no statistics attached (clumped SNP absent from the\n",
              nrow(miss)))
  cat("   sfFDR file, or that file is unreadable). Affected model x trait:\n")
  print(miss %>% count(variant, model, trait), n = 40)
}
