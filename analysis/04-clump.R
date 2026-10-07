##################################################
### File clumps sfFDR results into loci ###
##################################################
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 30L

# Paths and clumping parameters
BASE <- Sys.getenv("ANDY", "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")
RES_ROOT <- file.path(BASE, "sffdr/final_results")
# Tree-relative: runs never collide
OUT_ROOT <- file.path(normalizePath("."), "indep-discoveries-1000G")
TMP_ROOT <- file.path(OUT_ROOT, "_sumstats")
CLAIM <- file.path(OUT_ROOT, "_claims")
REF_TMPL <- file.path(DATA_GWASAI,
                      "baselineLD2.3/1000G_EUR_Phase3_plink/1000G.EUR.QC.%d")
PLINK <- "./plink2"

ARMS_TSV <- Sys.getenv("ARMS_TSV", "annotations/arms.tsv")
ARM_SET <- Sys.getenv("ARMS", "paper")

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

arms <- read_tsv(ARMS_TSV, show_col_types = FALSE)
if (ARM_SET == "paper") {
  arms <- arms %>% filter(paper == 1)
}
TYPES <- c("P", "fP")
P1 <- 5e-8
P2 <- 5e-8
R2 <- 0.01
KB <- 5000
MHC <- list(chr = 6, start = 24e6, end = 45e6)
PLINK_MEM <- 8000

walk(c(OUT_ROOT, TMP_ROOT, CLAIM), dir.create, showWarnings = FALSE,
     recursive = TRUE)

# fail now, not 280 jobs later
if (!file.exists(PLINK)) {
  stop("plink binary not found at ", PLINK)
}
ph <- tryCatch(system(paste(PLINK, "--version"), intern = TRUE,
                      ignore.stderr = TRUE),
               error = function(e) character(0),
               warning = function(w) character(0))
if (length(ph)) {
  message("plink: ", ph[1])
} else {
  message("note: could not read plink version; continuing")
}
miss_ref <- (1:22)[!file.exists(sprintf(paste0(REF_TMPL, ".bed"), 1:22))]
if (length(miss_ref)) {
  stop("missing 1000G reference for chr: ", paste(miss_ref, collapse = ","))
}
message("preflight OK")

# Job list
all_dirs <- arms %>%
  transmute(variant, model,
            path = if_else(variant == "joint", file.path(RES_ROOT, model),
                           file.path(RES_ROOT, "subset", model)))
dirs <- all_dirs %>% filter(dir.exists(path))

jobs <- dirs %>%
  mutate(f = map(path,
                 ~ list.files(.x, pattern = "\\.sffdr_final\\.gz$",
                              full.names = TRUE))) %>%
  unnest(f) %>%
  # Exclude '.full' fits: same P, different fP, so they clump as phantom
  # traits
  filter(!str_detect(basename(f), "\\.full\\.")) %>%
  mutate(trait = str_remove(basename(f), fixed(".sffdr_final.gz")),
         jid = paste(variant, model, trait, sep = "__")) %>%
  arrange(variant, model, trait)

cat("\ndirectories found:\n")
print(dirs %>% count(variant, model), n = 40)
missing <- all_dirs %>% filter(!dir.exists(path))
if (nrow(missing)) {
  cat("\nno results directory (skipped):\n")
  print(missing %>% select(variant, model), n = 40)
}

message(sprintf("\ntask %d/%d | %d trait-jobs total across %d directories",
                task_id, n_tasks, nrow(jobs), nrow(dirs)))

# Stagger the start; avoid contending for job 1
start <- if (nrow(jobs)) (task_id * max(1L, nrow(jobs) %/% n_tasks)) %% nrow(jobs) else 0L
order_idx <- if (nrow(jobs)) ((start + seq_len(nrow(jobs)) - 1L) %% nrow(jobs)) + 1L else integer(0)

claim_job <- function(jid) suppressWarnings(dir.create(file.path(CLAIM, jid)))

# One trait, per chromosome and p-value column
run_one <- function(row) {
  out_dir <- file.path(OUT_ROOT, row$variant, row$model)
  tmp_dir <- file.path(TMP_ROOT, row$variant, row$model)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  done <- function(pref) file.exists(paste0(pref, ".clumps")) ||
    file.exists(paste0(pref, ".empty"))

  d <- suppressMessages(read_tsv(row$f, show_col_types = FALSE,
        col_select = any_of(c("SNP", "CHR", "BP", "P", "fP", "fpvalues")),
        col_types = p_cols(row$f))) %>%
    mutate(across(any_of(c("P", "fP", "fpvalues")), exact_p))
  fcol <- intersect(c("fP", "fpvalues"), names(d))[1]
  if (is.na(fcol) || !"P" %in% names(d)) {
    message("  !! missing p-value column: ", row$f)
    return(invisible(NULL))
  }
  d <- d %>%
    transmute(SNP, CHR, BP, P, fP = .data[[fcol]]) %>%
    filter(!(CHR == MHC$chr & BP > MHC$start & BP < MHC$end))

  n_run <- 0L
  for (ty in TYPES) {
    dt <- d %>% filter(is.finite(.data[[ty]]), .data[[ty]] < P1)
    for (ch in 1:22) {
      pref <- file.path(out_dir, sprintf("%s-%d-%s", row$trait, ch, ty))
      if (done(pref)) {
        next
      }
      dc <- dt %>% filter(CHR == ch)
      if (!nrow(dc)) {
        file.create(paste0(pref, ".empty"))
        next
      }

      sf <- file.path(tmp_dir, sprintf("%s-%d-%s.assoc", row$trait, ch, ty))
      write_tsv(dc %>% select(SNP, all_of(ty)), sf)
      cmd <- paste(PLINK,
        "--bfile", sprintf(REF_TMPL, ch), "--clump", shQuote(sf),
        "--clump-snp-field SNP --clump-field", ty,
        "--clump-p1", P1, "--clump-p2", P2, "--clump-r2", R2, "--clump-kb", KB,
        "--rm-dup exclude-all", "--memory", PLINK_MEM, "--out", shQuote(pref),
        ">", shQuote(paste0(pref, ".stdout")), "2>",
        shQuote(paste0(pref, ".stderr")))
      st <- system(cmd)
      n_run <- n_run + 1L
      unlink(paste0(pref, ".failed"))
      if (st != 0) {
        # plink did not finish: record the failure, keep its output
        file.create(paste0(pref, ".failed"))
        message(sprintf("  !! plink exit %d: %s chr%d %s (see .stderr)",
                        st, row$trait, ch, ty))
      } else if (!file.exists(paste0(pref, ".clumps"))) {
        # nothing survived thresholding: record it
        file.create(paste0(pref, ".empty"))
      } else {
        unlink(c(paste0(pref, ".stdout"), paste0(pref, ".stderr")))
      }
      unlink(sf)
    }
  }
  message(sprintf("  done: %s/%s/%s  (%d plink calls)", row$variant, row$model,
                  row$trait, n_run))
}

n_done <- 0L
for (k in order_idx) {
  r <- jobs[k, ]
  if (!claim_job(r$jid)) {
    next
  } # another task has it
  n_done <- n_done + 1L
  message(sprintf("[%s] claim %d: %s | %s | %s",
                  format(Sys.time(), "%H:%M:%S"), n_done, r$variant, r$model,
                  r$trait))
  tryCatch(run_one(r), error = function(e)
    message("  !! ERROR ", r$trait, ": ", conditionMessage(e)))
}
message(sprintf("task %d finished %s | %d jobs processed",
                task_id, format(Sys.time()), n_done))
