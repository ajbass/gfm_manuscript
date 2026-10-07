##################################################
### File classifies loci as shared or distinct ###
##################################################
library(tidyverse)

# Collapsing first makes the counts reconcile with
# locus_counts_1000G.tsv.
#   Rscript 07-collapse-classify.R <id> <n_tasks>
#   Rscript 07-collapse-classify.R merge
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1 && args[1] == "merge") "merge" else "run"
task_id <- if (length(args) >= 1 && MODE == "run" && !is.na(suppressWarnings(as.integer(args[1])))) {
  as.integer(args[1])
} else {
  0L
}
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L

OUT_ROOT <- file.path(normalizePath("."), "indep-discoveries-1000G")
SIG_FILE <- file.path(OUT_ROOT, "sigresults_1000G.tsv")
# Cluster paths this analysis was run on.
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")
REF_TMPL <- file.path(DATA_GWASAI,
                      "baselineLD2.3/1000G_EUR_Phase3_plink/1000G.EUR.QC.%d")
PLINK <- "./plink2"
AH_WINDOW <- 100000
R2_THRESH <- 0.01
WINDOW_KB <- 5000
SHARD_DIR <- file.path(OUT_ROOT, "_collapse_shards")
TMP_ROOT <- file.path(OUT_ROOT, "_collapse_tmp")
dir.create(SHARD_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_ROOT, showWarnings = FALSE, recursive = TRUE)

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

# Collapse leads within 100 kb, most significant as lead; `pv` is the
# column to rank by
collapse_loci <- function(d, pv) {
  if (!nrow(d)) {
    return(d[0, ])
  }
  loci <- d %>% arrange(CHR, BP) %>%
    group_by(CHR) %>%
    mutate(locus = cumsum(c(TRUE, diff(BP) > AH_WINDOW))) %>%
    group_by(CHR, locus) %>%
    slice_min(.data[[pv]], n = 1, with_ties = FALSE) %>% # most significant
    ungroup() %>% select(-locus)
  return(loci)
}

# Merge shards into the classification table
if (MODE == "merge") {
  sh <- list.files(SHARD_DIR, pattern = "\\.tsv$", full.names = TRUE)
  if (!length(sh)) {
    stop("no shards in ", SHARD_DIR)
  }
  loci <- map_dfr(sh, ~ suppressMessages(read_tsv(.x, show_col_types = FALSE)))
  write_tsv(loci, file.path(OUT_ROOT, "collapse_classification.tsv"))

  cmp <- loci %>% count(trait, Classification) %>%
    pivot_wider(names_from = Classification, values_from = n, values_fill = 0)
  for (nm in c("Shared", "Distinct_P", "Distinct_fP")) {
    if (!nm %in% names(cmp)) cmp[[nm]] <- 0L
  }
  cmp <- cmp %>% mutate(n_P = Shared + Distinct_P, n_fP = Shared + Distinct_fP,
                        net_gain = Distinct_fP - Distinct_P,
                        pct_change = round(100 * (n_fP - n_P) / pmax(n_P, 1),
                                           2))
  write_tsv(cmp, file.path(OUT_ROOT, "collapse_comparison.tsv"))

  cat("=== pooled (100 kb-collapsed loci, joint baselineLD) ===\n")
  print(cmp %>% summarise(traits = n(), shared = sum(Shared),
                          only_fP = sum(Distinct_fP), only_P = sum(Distinct_P),
                          net = sum(net_gain), n_P = sum(n_P),
                          n_fP = sum(n_fP)),
        width = Inf)

  sa <- read_tsv(file.path(OUT_ROOT, "locus_counts_1000G.tsv"),
                 show_col_types = FALSE) %>%
    filter(variant == "joint", str_detect(model, "baselineLD2.2$")) %>%
    transmute(ts = str_remove(trait, "^UKBB-[^_]+_"), sa_nP = n_P,
              sa_nfP = n_fP)
  chk <- cmp %>% mutate(ts = str_remove(trait, "^UKBB-[^_]+_")) %>%
    left_join(sa, by = "ts")
  cat(sprintf("\nCONSISTENCY with locus_counts_1000G.tsv (should be exact):\n"))
  cat(sprintf("  n_P  matches: %d/%d    n_fP matches: %d/%d\n",
              sum(chk$n_P == chk$sa_nP, na.rm = TRUE), nrow(chk),
              sum(chk$n_fP == chk$sa_nfP, na.rm = TRUE), nrow(chk)))
  off <- chk %>% filter(n_P != sa_nP | n_fP != sa_nfP) %>%
    select(trait, n_P, sa_nP, n_fP, sa_nfP)
  if (nrow(off)) {
    cat("  mismatches:\n")
    print(off, n = 30)
  }

  cat("\n=== net gain per trait ===\n")
  print(cmp %>% select(trait, Shared, Distinct_P, Distinct_fP, net_gain,
                       pct_change) %>%
          arrange(desc(net_gain)), n = 40)
  quit(save = "no")
}

# Lead sets
sig <- read_tsv(SIG_FILE, col_types = p_cols(SIG_FILE), show_col_types = FALSE) %>%
  mutate(across(any_of(c("P", "fP")), exact_p)) %>%
  filter(variant == "joint", model == "baselineLD2.2", !is.na(CHR), !is.na(BP))
jobs <- sig %>% distinct(trait) %>% arrange(trait) %>%
  mutate(k = row_number()) %>% filter(n_tasks <= 1 | (k %% n_tasks) == task_id)
message(sprintf("task %d/%d: %d traits", task_id, n_tasks, nrow(jobs)))
if (!nrow(jobs)) {
  message("nothing to do")
  quit(save = "no")
}

# r2 between SNP lists; --r2-unphased writes .vcor
ld_r2 <- function(chr, snps, tag) {
  EMPTY <- tibble(A = character(0), B = character(0), r2 = numeric(0))
  if (length(snps) < 2) {
    return(EMPTY)
  }
  tmp <- file.path(TMP_ROOT, sprintf("%s_chr%d_%d", tag, chr, task_id))
  writeLines(unique(snps), paste0(tmp, ".snps"))
  cmd <- paste(PLINK, "--bfile", sprintf(REF_TMPL, chr),
               "--extract", shQuote(paste0(tmp, ".snps")),
               "--r2-unphased", "--ld-window-r2", "0",
               "--ld-window-kb", WINDOW_KB, "--ld-window", 99999,
               "--out", shQuote(tmp), ">", shQuote(paste0(tmp, ".log")), "2>&1")
  system(cmd)
  f <- paste0(tmp, ".vcor")
  if (!file.exists(f)) {
    unlink(list.files(TMP_ROOT, pattern = basename(tmp), full.names = TRUE))
    return(EMPTY)
  }
  d <- suppressMessages(read_table(f, show_col_types = FALSE))
  unlink(list.files(TMP_ROOT, pattern = basename(tmp), full.names = TRUE))
  if (!nrow(d) || !all(c("ID_A", "ID_B", "UNPHASED_R2") %in% names(d))) {
    return(EMPTY)
  }
  return(tibble(A = d$ID_A, B = d$ID_B, r2 = as.numeric(d$UNPHASED_R2)))
}

classify_trait <- function(tr) {
  d <- sig %>% filter(trait == tr)
  # Step 1: collapse each method to 100 kb loci
  P_loci <- collapse_loci(d %>% filter(type == "P") %>% distinct(SNP, CHR, BP,
                                                                 P, fP), "P")
  fP_loci <- collapse_loci(d %>% filter(type == "fP") %>% distinct(SNP, CHR, BP,
                                                                   P, fP), "fP")

  # Step 2: r2-match representatives per chromosome
  out <- map_dfr(sort(unique(c(P_loci$CHR, fP_loci$CHR))), function(ch) {
    Pc <- P_loci %>% filter(CHR == ch)
    fPc <- fP_loci %>% filter(CHR == ch)
    ld <- ld_r2(ch, unique(c(Pc$SNP, fPc$SNP)),
                str_replace_all(tr, "[^A-Za-z0-9]", ""))

    linked <- ld %>% filter(r2 > R2_THRESH) %>%
      filter((A %in% Pc$SNP & B %in% fPc$SNP) | (A %in% fPc$SNP & B %in% Pc$SNP)) %>%
      transmute(P_snp = if_else(A %in% Pc$SNP, A, B),
                fP_snp = if_else(A %in% fPc$SNP, A, B), r2)
    same <- intersect(Pc$SNP, fPc$SNP)
    if (length(same)) linked <- bind_rows(linked,
                                          tibble(P_snp = same, fP_snp = same,
                                                 r2 = 1))

    # One-to-one greedy match by descending r2
    linked <- linked %>% arrange(desc(r2))
    usedP <- character(0)
    usedF <- character(0)
    m <- list()
    for (i in seq_len(nrow(linked))) {
      p <- linked$P_snp[i]
      f <- linked$fP_snp[i]
      if (p %in% usedP || f %in% usedF) {
        next
      }
      usedP <- c(usedP, p)
      usedF <- c(usedF, f)
      m[[length(m) + 1]] <- tibble(P_snp = p, fP_snp = f)
    }
    matched <- if (length(m)) bind_rows(m) else tibble(P_snp = character(0),
                                                       fP_snp = character(0))

    # a shared locus is keyed by its standard lead. fP_lead names the
    # functional lead it was matched to, which can be another SNP
    bind_rows(
      Pc %>% inner_join(matched, by = c("SNP" = "P_snp")) %>%
        transmute(CHR, POS = BP, lead_SNP = SNP, fP_lead = fP_snp,
                  Classification = "Shared"),
      Pc %>% filter(!SNP %in% matched$P_snp) %>%
        transmute(CHR, POS = BP, lead_SNP = SNP, fP_lead = NA_character_,
                  Classification = "Distinct_P"),
      fPc %>% filter(!SNP %in% matched$fP_snp) %>%
        transmute(CHR, POS = BP, lead_SNP = SNP, fP_lead = SNP,
                  Classification = "Distinct_fP"))
  })
  out <- out %>% mutate(trait = tr, model = "baselineLD2.2", variant = "joint")
  return(out)
}

res <- map_dfr(jobs$trait, function(tr) {
  message(sprintf("  [%s] %s", format(Sys.time(), "%H:%M:%S"), tr))
  classify_trait(tr)
})

if (n_tasks > 1L) {
  f <- file.path(SHARD_DIR, sprintf("shard_%03d.tsv", task_id))
  write_tsv(res, f)
  cat("wrote shard ", f, "  (", nrow(res), " rows)\n", sep = "")
  cat("when all done:  Rscript 07-collapse-classify.R merge\n")
} else {
  write_tsv(res, file.path(OUT_ROOT, "collapse_classification.tsv"))
  cat("wrote collapse_classification.tsv  (", nrow(res), " rows)\n", sep = "")
  cat("now run:  Rscript 07-collapse-classify.R merge\n")
}
