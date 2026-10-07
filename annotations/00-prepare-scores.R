##################################################
### File builds the annotation release ###
##################################################
library(tidyverse)
library(arrow)

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) to build somewhere else.
setwd(Sys.getenv("DATA_GWASAI", "/home/ab3105/rds/hpc-work/data_gwasAI"))

# Inputs, all under DATA_GWASAI (downloads and the AlphaGenome API output;
# see the README)
AG_SCORES <- "./foundational_models/alphagenome/scores.csv"
BZ_RAW <- "./foundational_models/borzoi/borzoi_102_annotation_set/annot_raw/borzoi.%d.annot.parquet"
BZ_SNPS <- "./foundational_models/borzoi/borzoi_102_annotation_set/SNP_LIST_1000G.tsv"
GS_SCORES <- "./foundational_models/gpn-star/scores.parquet"
BL_RAW <- "./foundational_models/baselineLD2.2/annot/baselineLD.%d.annot.gz"
# baselineLD v2.3, which carries the two Zoonomia columns
ZO_RAW <- "./foundational_models/zoonomia/annot/baselineLD.%d.annot.gz"
BIM <- "./baselineLD2.3/1000G_EUR_Phase3_plink/1000G.EUR.QC.%d.bim"
# the API export, or a gzipped copy of it
if (!file.exists(AG_SCORES)) {
  AG_SCORES <- paste0(AG_SCORES, ".gz")
}
# The Atlas variant-impact release (unzipped file and .tbi, from
# https://deepmind.google.com/science/alphagenome), its GRCh37 to GRCh38
# map (SNP and BP38 columns) and tabix.
AVI_RAW <- Sys.getenv("AVI_RAW", "./foundational_models/avi/raw")
MAP38 <- Sys.getenv("MAP38", "./baselineLD2.3/panukbb_37_to_38_map.tsv")
TABIX <- Sys.getenv("TABIX", unset = Sys.which("tabix"))

# Outputs, all under the release directory: the SNP list, the per-source
# intermediates, then the release files, one per model and chromosome on
# one SNP set, with their column guide
REL <- "./foundational_models/release"
SNPS <- file.path(REL, "snps.tsv")
SRC_FILE <- file.path(REL, "sources/%s.%d.parquet")
REL_FILE <- file.path(REL, "%s.%d.parquet")
COLUMNS <- file.path(REL, "columns.tsv")

KEYS <- c("SNP", "CHR", "BP", "A1", "A2")
AG_TYPES <- c("ATAC", "CAGE", "CHIP_HISTONE", "CHIP_TF", "CONTACT_MAPS",
              "DNASE", "PROCAP", "RNA_SEQ")
GS_COLS <- c("GPN-STAR-P", "GPN-STAR-M", "GPN-STAR-V")
ZO_COLS <- c("zoonomia_mammals", "zoonomia_primates")
# the API export, typed so no score is guessed from a leading run of
# missing values; the columns not used are read as text
AG_SPEC <- do.call(cols, c(list(.default = col_character(),
                                variant_id = col_character()),
                           setNames(rep(list(col_double()), length(AG_TYPES)),
                                    paste0(AG_TYPES, "_raw"))))
# GPN-Star indicators: the top 0.1% of each score genome-wide, the GPN-Star
# paper's own cut; ties broken by a seeded draw (random_state = 42
# upstream).
GS_QUANTILES <- c(bin01 = 0.001)
GS_TIE_SEED <- 42
GS_BIN <- as.vector(outer(GS_COLS, names(GS_QUANTILES), paste, sep = "_"))
AVI_EXPECT <- c("CHROM", "POS", "REF", "ALT", "raw_score", "PHRED")
# The sources every release SNP must carry; the indicators and AVI follow
SOURCES <- c("baselineLD", "zoonomia", "alphagenome", "borzoi", "gpn_star")

# Log-modulus transform for the Borzoi components
log_modulus <- function(x) sign(x) * log1p(abs(x))

# Standardize with a mean and SD accumulated over the whole genome
standardize <- function(x, mu, sdv) (x - mu) / sdv

src_file <- function(src, chr) sprintf(SRC_FILE, src, chr)
done_src <- function(src) all(file.exists(src_file(src, 1:22)))
read_src <- function(src, chr) read_parquet(src_file(src, chr))

# The downloaded annotation files, CHR and BP as integers so the join
# keeps the .bim's integer keys (readr would guess doubles; PolyFun
# matches on CHR.BP.A1.A2 as text, and 5.0 is not 5)
ANNOT_TYPES <- cols(CHR = col_integer(), BP = col_integer(),
                    SNP = col_character(), .default = col_double())

# The 1000G EUR keys, one row per SNP
read_bim <- function(chr) {
  bim <- read_tsv(sprintf(BIM, chr),
                  col_names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
                  col_types = "icdicc", progress = FALSE) %>%
    select(all_of(KEYS)) %>% distinct(SNP, .keep_all = TRUE)
  return(bim)
}

comp <- function(x) chartr("ACGT", "TGCA", x)
pick_panel_row <- function(d, keys) {
  d <- d %>%
    mutate(file_row = row_number()) %>%
    inner_join(keys %>% select(SNP, BP_p = BP, A1_p = A1, A2_p = A2),
               by = "SNP") %>%
    mutate(is_panel = (A1 == A1_p & A2 == A2_p) | (A1 == A2_p & A2 == A1_p) |
             (comp(A1) == A1_p & comp(A2) == A2_p) |
             (comp(A1) == A2_p & comp(A2) == A1_p)) %>%
    group_by(SNP) %>%
    mutate(first_row = min(file_row)) %>%
    ungroup() %>%
    arrange(SNP, desc(is_panel), file_row) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    mutate(from_other_row = file_row != first_row,
           BP = BP_p, A1 = A1_p, A2 = A2_p) %>%
    select(-BP_p, -A1_p, -A2_p, -file_row, -first_row)
  return(d)
}
allele_tally <- list()

dir.create(file.path(REL, "sources"), showWarnings = FALSE, recursive = TRUE)

## PART 1: the SNPs scored by every model
# Borzoi's list is already the MAF > 0.01 1000G EUR set. GPN-Star and
# AlphaGenome contribute the SNPs they scored without missing values
if (file.exists(SNPS)) {
  message("SNP list exists: ", SNPS)
} else {
  message("=== PART 1: SNP overlap ===")
  bz_snps <- read_tsv(BZ_SNPS, show_col_types = FALSE)$SNP
  gs <- read_parquet(GS_SCORES, col_select = c("rsid", GS_COLS))
  gs_snps <- gs$rsid[complete.cases(gs)]
  rm(gs)
  ag <- read_csv(AG_SCORES, col_select = c("variant_id", "RNA_SEQ_raw"),
                 col_types = AG_SPEC)
  ag_snps <- ag$variant_id[!is.na(ag$RNA_SEQ_raw)]
  rm(ag)
  snps <- Reduce(intersect, list(bz_snps, gs_snps, ag_snps))
  message(sprintf("  borzoi %d, gpn-star %d, alphagenome %d -> %d shared",
                  length(bz_snps), length(gs_snps), length(ag_snps),
                  length(snps)))
  write_tsv(tibble(SNP = snps), SNPS)
}
snps <- read_tsv(SNPS, show_col_types = FALSE)$SNP
message(sprintf("%d SNPs in the overlap (7,674,062 in the paper's run)",
                length(snps)))
bim_all <- map_dfr(1:22, read_bim)

## PART 2a: baselineLD v2.2 on the overlap with the 1000G keys
if (done_src("baselineLD")) {
  message("baselineLD exists")
} else {
  message("=== PART 2a: baselineLD ===")
  for (chr in 1:22) {
    d <- read_tsv(sprintf(BL_RAW, chr), col_types = ANNOT_TYPES,
                  progress = FALSE) %>%
      filter(SNP %in% snps) %>% select(-CM)
    d <- read_bim(chr) %>% inner_join(d, by = c("CHR", "SNP", "BP")) %>%
      filter(!is.na(base)) %>%
      select(all_of(KEYS), base, everything())
    write_parquet(d, src_file("baselineLD", chr))
    message(sprintf("  chr%-2d %d SNPs, %d annotations", chr, nrow(d),
                    ncol(d) - length(KEYS) - 1L))
  }
}

## PART 2b: Zoonomia, the two constraint columns of baselineLD v2.3
if (done_src("zoonomia")) {
  message("zoonomia exists")
} else {
  message("=== PART 2b: zoonomia ===")
  for (chr in 1:22) {
    d <- read_tsv(sprintf(ZO_RAW, chr), col_types = ANNOT_TYPES,
                  progress = FALSE) %>%
      filter(SNP %in% snps) %>%
      select(CHR, SNP, BP, base, all_of(ZO_COLS))
    d <- read_bim(chr) %>% inner_join(d, by = c("CHR", "SNP", "BP")) %>%
      filter(!is.na(base)) %>%
      select(all_of(KEYS), all_of(ZO_COLS))
    write_parquet(d, src_file("zoonomia", chr))
    message(sprintf("  chr%-2d %d SNPs", chr, nrow(d)))
  }
}

## PART 2c: AlphaGenome, log then standardize, genome-wide
# The scores are the maximum absolute effect over the tracks of each
# output type. Seven are log(1 + x); CONTACT_MAPS is log(x), as run
if (done_src("alphagenome")) {
  message("alphagenome exists")
} else {
  message("=== PART 2c: alphagenome ===")
  # scores.csv is a pandas export: an unnamed index column, quantile and
  # raw columns per output type, then variant_id (rsID) and positions
  ag <- read_csv(AG_SCORES, col_types = AG_SPEC, name_repair = "unique_quiet")
  need <- c("variant_id", paste0(AG_TYPES, "_raw"))
  if (!all(need %in% names(ag))) {
    stop("scores.csv lacks ", paste(setdiff(need, names(ag)), collapse = ", "))
  }
  ag <- ag %>%
    select(SNP = variant_id, all_of(paste0(AG_TYPES, "_raw"))) %>%
    filter(!is.na(RNA_SEQ_raw), SNP %in% snps)
  n_na <- colSums(is.na(ag))
  if (any(n_na > 0)) {
    stop("missing scores after the RNA_SEQ filter: ",
         paste(names(n_na)[n_na > 0], collapse = ", "))
  }
  # standardized over the scored SNPs in the overlap, before the keys are
  # attached, as run
  for (ty in AG_TYPES) {
    raw <- ag[[paste0(ty, "_raw")]]
    v <- if (ty == "CONTACT_MAPS") log(raw) else log1p(raw)
    if (any(!is.finite(v))) {
      stop(ty, ": non-finite values after the log transform")
    }
    ag[[ty]] <- standardize(v, mean(v), sd(v))
  }
  ag <- bim_all %>% inner_join(ag, by = "SNP") %>%
    distinct(SNP, .keep_all = TRUE) %>%
    select(all_of(KEYS), all_of(AG_TYPES))
  for (chr in 1:22) {
    d <- ag %>% filter(CHR == chr)
    write_parquet(d, src_file("alphagenome", chr))
    message(sprintf("  chr%-2d %d SNPs", chr, nrow(d)))
  }
  rm(ag)
}

## PART 2d: Borzoi, log-modulus then standardize, genome-wide, two passes
# The moments are accumulated over the rows of the raw files, duplicated
# rsIDs included, as run.
if (done_src("borzoi")) {
  message("borzoi exists")
} else {
  message("=== PART 2d: borzoi ===")
  bz_cols <- setdiff(names(open_dataset(sprintf(BZ_RAW, 1L))$schema),
                     c(KEYS, "base", "CM", "MAF"))
  message("  ", length(bz_cols), " components")
  if (length(bz_cols) != 102) {
    warning("expected 102 Borzoi components, found ", length(bz_cols))
  }
  s1 <- setNames(numeric(length(bz_cols)), bz_cols)
  s2 <- s1
  n <- 0
  for (chr in 1:22) {
    x <- read_parquet(sprintf(BZ_RAW, chr)) %>%
      filter(SNP %in% snps) %>%
      select(all_of(bz_cols)) %>%
      mutate(across(everything(), log_modulus))
    if (anyNA(x)) {
      stop("chr", chr, ": missing Borzoi scores on the overlap")
    }
    s1 <- s1 + colSums(x)
    s2 <- s2 + colSums(x^2)
    n <- n + nrow(x)
    rm(x)
  }
  bz_mu <- s1 / n
  bz_sd <- sqrt((s2 - n * bz_mu^2) / (n - 1))
  n_other <- 0L
  n_none <- 0L
  for (chr in 1:22) {
    # the raw file carries 1000G keys, but for a site with two alternate
    # alleles its first row can be the allele the panel does not carry
    raw <- read_parquet(sprintf(BZ_RAW, chr)) %>% filter(SNP %in% snps)
    d <- pick_panel_row(raw, read_bim(chr))
    other <- sum(d$from_other_row & d$is_panel)
    n_drop <- sum(!d$is_panel)
    d <- d %>% filter(is_panel)
    for (cl in bz_cols) {
      d[[cl]] <- standardize(log_modulus(d[[cl]]), bz_mu[[cl]], bz_sd[[cl]])
    }
    d <- d %>% select(all_of(KEYS), all_of(bz_cols))
    write_parquet(d, src_file("borzoi", chr))
    message(sprintf("  chr%-2d %d SNPs; %d taken from another row of the file, %d dropped (scored for neither panel allele)",
                    chr, nrow(d), other, n_drop))
    n_other <- n_other + other
    n_none <- n_none + n_drop
    rm(raw)
  }
  allele_tally$borzoi <- c(other_row = n_other, none = n_none)
}

## PART 2e: GPN-Star
if (done_src("gpn_star")) {
  message("gpn_star exists")
} else {
  message("=== PART 2e: gpn_star ===")
  gs <- read_parquet(GS_SCORES, col_select = c("rsid", "ref", "alt", GS_COLS)) %>%
    rename(SNP = rsid, A1 = alt, A2 = ref) %>%
    filter(SNP %in% snps)
  if (anyNA(gs[GS_COLS])) {
    stop("missing GPN-Star scores on the overlap")
  }
  # standardized over the SNPs that carry 1000G keys, as run; then, of an
  # rsID's rows, the panel's variant (the file's alleles are GRCh38)
  gs <- bim_all %>% select(SNP, CHR) %>% inner_join(gs, by = "SNP") %>%
    mutate(across(all_of(GS_COLS), ~ standardize(.x, mean(.x), sd(.x)))) %>%
    pick_panel_row(bim_all)
  other <- sum(gs$from_other_row & gs$is_panel)
  n_drop <- sum(!gs$is_panel)
  gs <- gs %>% filter(is_panel)
  allele_tally$gpn_star <- c(other_row = other, none = n_drop)
  message(sprintf("  %d SNPs taken from another row of the file, %d dropped (scored for neither panel allele)",
                  other, n_drop))
  gs <- gs %>% select(all_of(KEYS), all_of(GS_COLS))
  for (chr in 1:22) {
    d <- gs %>% filter(CHR == chr)
    write_parquet(d, src_file("gpn_star", chr))
    message(sprintf("  chr%-2d %d SNPs", chr, nrow(d)))
  }
  rm(gs)
}

## PART 3: GPN-Star indicators
bin_done <- done_src("gpn_star_bin") &&
  identical(setdiff(names(open_dataset(src_file("gpn_star_bin", 1L))$schema),
                    "SNP"), GS_BIN)
if (bin_done) {
  message("gpn_star_bin exists")
} else {
  message("=== PART 3: GPN-Star indicators ===")
  gs_all <- map_dfr(1:22, function(chr) {
    keep <- Reduce(intersect, map(SOURCES, function(src) {
      read_parquet(src_file(src, chr), col_select = "SNP")$SNP
    }))
    read_src("gpn_star", chr) %>% filter(SNP %in% keep) %>%
      select(SNP, CHR, all_of(GS_COLS))
  })
  n_gs <- nrow(gs_all)
  message(sprintf("  %d SNPs genome-wide", n_gs))
  set.seed(GS_TIE_SEED)
  gs_tie <- runif(n_gs)
  for (ql in names(GS_QUANTILES)) {
    n_sel <- as.integer(floor(n_gs * GS_QUANTILES[[ql]]))
    for (sc in GS_COLS) {
      v <- gs_all[[sc]]
      ord <- order(-v, gs_tie)
      s <- logical(n_gs)
      s[ord[seq_len(n_sel)]] <- TRUE
      gs_all[[paste0(sc, "_", ql)]] <- as.integer(s)
      message(sprintf("  %-18s n = %d, threshold %.4f", paste0(sc, "_", ql),
                      n_sel, v[ord[n_sel]]))
    }
  }
  for (chr in 1:22) {
    d <- gs_all %>% filter(CHR == chr) %>% select(-CHR, -all_of(GS_COLS))
    write_parquet(d, src_file("gpn_star_bin", chr))
  }
  rm(gs_all)
}

## PART 4: AVI, the Atlas variant-impact score
avi_fs <- list.files(AVI_RAW, pattern = "\\.(bgz|gz)$", recursive = TRUE,
                     full.names = TRUE)
avi_fs <- avi_fs[!grepl("\\.tbi$", avi_fs)]
if (done_src("avi")) {
  message("avi exists")
} else if (!length(avi_fs) || !nzchar(TABIX) || !file.exists(MAP38)) {
  message("=== PART 4 SKIPPED: AVI needs the Atlas file under ", AVI_RAW,
          ", tabix on PATH (module load htslib) and the 37->38 map at ",
          MAP38, "; the release is built without the AVI column ===")
} else {
  message("=== PART 4: AVI ===")
  f_avi <- avi_fs[which.max(file.size(avi_fs))]
  if (!file.exists(paste0(f_avi, ".tbi"))) {
    stop("no .tbi beside ", f_avi)
  }
  message("  Atlas file: ", f_avi)
  avi_h <- suppressWarnings(system(sprintf("zcat %s 2>/dev/null | head -1",
                                           shQuote(f_avi)),
                                   intern = TRUE))
  avi_cols <- strsplit(sub("^#+\\s*", "", avi_h[1]), "\t")[[1]]
  if (!identical(avi_cols, AVI_EXPECT)) {
    stop("unexpected Atlas header: ", paste(avi_cols, collapse = ", "))
  }

  avi_l1 <- readLines(MAP38, n = 1)
  avi_mp <- if (grepl("\t", avi_l1)) {
    read_tsv(MAP38, show_col_types = FALSE)
  } else {
    read_table(MAP38, show_col_types = FALSE)
  }
  names(avi_mp) <- toupper(names(avi_mp))
  avi_pick <- function(cands) {
    i <- match(cands, names(avi_mp))
    return(names(avi_mp)[i[!is.na(i)][1]])
  }
  avi_csnp <- avi_pick(c("SNP", "RSID", "RS", "ID"))
  avi_cbp <- avi_pick(c("BP38", "POS38"))
  if (any(is.na(c(avi_csnp, avi_cbp)))) {
    stop("map: need SNP and BP38 columns; found: ",
         paste(names(avi_mp), collapse = ", "))
  }
  avi_map <- avi_mp %>%
    transmute(SNP = .data[[avi_csnp]], bp38 = as.integer(.data[[avi_cbp]])) %>%
    distinct(SNP, .keep_all = TRUE)
  message("  map rows: ", nrow(avi_map))
  rm(avi_mp)

  # pass 1: extract raw_score per chromosome
  avi_sc <- vector("list", 22)
  avi_st <- vector("list", 22)
  for (chr in 1:22) {
    sp <- read_parquet(src_file("baselineLD", chr), col_select = KEYS)
    lk <- sp %>% inner_join(avi_map, by = "SNP")
    n_nomap <- nrow(sp) - nrow(lk)
    # Cached per chromosome, so a job resumes; delete extract_chr*.tsv
    # to re-extract
    cache <- file.path(AVI_RAW, sprintf("extract_chr%d.tsv", chr))
    if (!file.exists(cache)) {
      posf <- tempfile()
      tmpo <- paste0(cache, ".tmp")
      writeLines(as.character(sort(unique(lk$bp38))), posf)
      st <- system(sprintf(
        "%s %s chr%d | awk -v pf=%s 'BEGIN{while((getline p < pf)>0) w[p]=1} w[$2]' > %s",
        shQuote(TABIX), shQuote(f_avi), chr, shQuote(posf), shQuote(tmpo)))
      if (st != 0) {
        stop("streaming extraction failed on chr", chr)
      }
      file.rename(tmpo, cache)
      unlink(posf)
    } else {
      message(sprintf("  chr%-2d using cached extraction", chr))
    }
    av <- read_tsv(cache, col_names = AVI_EXPECT, comment = "#",
                   show_col_types = FALSE,
                   col_types = cols(CHROM = "c", POS = "i", REF = "c",
                                    ALT = "c",
                                    raw_score = "d", PHRED = "d")) %>%
      select(bp38 = POS, REF, ALT, raw_score)

    j <- lk %>% inner_join(av, by = "bp38", relationship = "many-to-many")
    n_nopos <- nrow(lk) - n_distinct(j$SNP)
    direct <- j %>%
      filter((REF == A2 & ALT == A1) | (REF == A1 & ALT == A2)) %>%
      distinct(SNP, .keep_all = TRUE)
    # palindromic SNPs are skipped
    flip <- j %>% filter(!SNP %in% direct$SNP, A2 != comp(A1)) %>%
      mutate(f1 = comp(A1), f2 = comp(A2)) %>%
      filter((REF == f2 & ALT == f1) | (REF == f1 & ALT == f2)) %>%
      distinct(SNP, .keep_all = TRUE)
    avi_sc[[chr]] <- bind_rows(direct %>% select(SNP, raw_score),
                               flip %>% select(SNP, raw_score))
    n_unm <- nrow(lk) - n_nopos - nrow(avi_sc[[chr]])
    avi_st[[chr]] <- tibble(chr = chr, spine = nrow(sp), no_map = n_nomap,
                            no_atlas_row = n_nopos, strand_flipped = nrow(flip),
                            allele_unmatched = n_unm,
                            matched = nrow(avi_sc[[chr]]))
    message(sprintf("  chr%-2d spine %6d | no_map %4d | no_row %4d | flipped %4d | unmatched %4d | matched %6d",
                    chr, nrow(sp), n_nomap, n_nopos, nrow(flip), n_unm,
                    nrow(avi_sc[[chr]])))
  }
  avi_st <- bind_rows(avi_st)
  avi_all <- bind_rows(avi_sc)
  cat("\nAVI extraction totals\n")
  print(avi_st %>% summarise(across(-chr, sum)), width = Inf)

  # standardize genome-wide: (raw - mean) / SD
  avi_mu <- mean(avi_all$raw_score)
  avi_sd <- sd(avi_all$raw_score)
  message(sprintf("standardizing: mean %.5f, sd %.5f over %d matched SNPs",
                  avi_mu, avi_sd, nrow(avi_all)))
  avi_all <- avi_all %>% mutate(AVI = (raw_score - avi_mu) / avi_sd) %>%
    select(SNP, AVI)
  avi_imp <- 100 * (sum(avi_st$spine) - nrow(avi_all)) / sum(avi_st$spine)
  message(sprintf("imputing 0 (the mean) for %.3f%% of spine SNPs", avi_imp))
  if (avi_imp > 1) warning("more than 1% imputed -- check the map and allele handling")

  # pass 2: the column, 0 where unmatched
  for (chr in 1:22) {
    d <- read_parquet(src_file("baselineLD", chr), col_select = "SNP") %>%
      left_join(avi_all, by = "SNP") %>%
      mutate(AVI = replace_na(AVI, 0))
    write_parquet(d, src_file("avi", chr))
  }
  cat("\nAVI: matched", nrow(avi_all), "of", sum(avi_st$spine),
      sprintf("SNPs (%.3f%% set to 0),", avi_imp),
      sum(avi_st$strand_flipped), "by strand flip\n")
}

## PART 5: the release files
MODELS <- list(baselineLD = "baselineLD", borzoi = "borzoi",
               alphagenome = "alphagenome", zoonomia = "zoonomia",
               gpn_star = c("gpn_star", "gpn_star_bin"), avi = "avi")
GROUPS <- c(baselineLD = "baselineLD", borzoi = "borzoi",
            alphagenome = "alphagenome", zoonomia = "zoonomia",
            gpn_star = "gpnstar", gpn_star_bin = "gpnstar", avi = "avi")
TRANSFORMS <- c(
  baselineLD = "as published in baselineLD v2.2",
  borzoi = "sign(x) log(1 + |x|), then standardized genome-wide; the panel's allele where the file lists two",
  alphagenome = "log(1 + x) (CONTACT_MAPS: log(x)), then standardized genome-wide",
  zoonomia = "as published in baselineLD v2.3",
  gpnstar = "standardized genome-wide; the panel's allele where the file lists two",
  gpnstar_bin01 = "1 for the top 0.1% of the standardized score genome-wide, else 0",
  avi = "Atlas raw_score standardized over the matched SNPs; 0 where unmatched"
)
message("=== PART 5: the release ===")
have <- names(GROUPS)[vapply(names(GROUPS), done_src, logical(1))]
if (!all(c(SOURCES, "gpn_star_bin") %in% have)) {
  stop("sources missing: ",
       paste(setdiff(c(SOURCES, "gpn_star_bin"), have), collapse = ", "))
}
models <- MODELS[map_lgl(MODELS, ~ all(.x %in% have))]
rel_file <- function(m, chr) sprintf(REL_FILE, m, chr)
src_cols <- function(src, chr = 1L) {
  sc <- names(open_dataset(src_file(src, chr))$schema)
  return(setdiff(sc, c(KEYS, "base")))
}

# the column guide, from the chromosome 1 schemas
guide <- map_dfr(names(models), function(m) {
  map_dfr(models[[m]], function(src) {
    sc <- src_cols(src)
    g <- if (src == "gpn_star_bin") {
      paste0(GROUPS[[src]], str_extract(sc, "_bin[^_]+$"))
    } else {
      rep(GROUPS[[src]], length(sc))
    }
    tibble(column = sc, file = m, group = g, transform = unname(TRANSFORMS[g]))
  })
})
if (anyDuplicated(guide$column)) {
  stop("column names shared between sources: ",
       paste(guide$column[duplicated(guide$column)], collapse = ", "))
}
message(sprintf("  %d columns: %s", nrow(guide),
                paste(sprintf("%s %d", names(table(guide$group)),
                              table(guide$group)), collapse = ", ")))

built <- vapply(1:22, function(chr) {
  all(file.exists(rel_file(names(models), chr)))
}, logical(1))
if (any(built)) {
  chr1 <- which(built)[1]
  for (m in names(models)) {
    on_disk <- setdiff(names(open_dataset(rel_file(m, chr1))$schema),
                       c(KEYS, "base"))
    if (!identical(on_disk, guide$column[guide$file == m])) {
      stop("the release on disk has other columns than the sources here; ",
           "remove ", REL, "/*.parquet to rebuild it")
    }
  }
}
write_tsv(guide, COLUMNS)
rel_snps <- vector("list", 22)
for (chr in 1:22) {
  if (built[chr]) {
    rel_snps[[chr]] <- read_parquet(rel_file(names(models)[1], chr),
                                    col_select = "SNP")$SNP
    next
  }
  srcs <- map(set_names(have), read_src, chr)
  keep <- Reduce(intersect, map(srcs, "SNP"))
  rel_snps[[chr]] <- keep
  # integer keys: PolyFun matches them to the .bim as text
  sp <- srcs$baselineLD %>% filter(SNP %in% keep) %>%
    select(all_of(KEYS)) %>%
    mutate(CHR = as.integer(CHR), BP = as.integer(BP), base = 1)
  for (m in names(models)) {
    d <- sp
    for (src in models[[m]]) {
      d <- d %>% inner_join(srcs[[src]] %>% select(SNP, all_of(src_cols(src))),
                            by = "SNP")
    }
    if (nrow(d) != nrow(sp) || anyDuplicated(d$SNP)) {
      stop("chr", chr, " ", m, ": rows do not align with the SNP set")
    }
    if (!identical(setdiff(names(d), c(KEYS, "base")),
                   guide$column[guide$file == m])) {
      stop("chr", chr, " ", m, ": columns differ from chromosome 1")
    }
    if (anyNA(d)) {
      stop("chr", chr, " ", m, ": missing values")
    }
    write_parquet(d, rel_file(m, chr))
  }
  rm(srcs)
  message(sprintf("  chr%-2d %d SNPs, %d files", chr, nrow(sp), length(models)))
}
# the SNP list is the release's rows: the overlap less the SNPs a source
# scores for neither allele of the panel's variant
rel_snps <- unlist(rel_snps)
if (!setequal(rel_snps, snps)) {
  message(sprintf("  release: %d SNPs, %d of the overlap dropped; snps.tsv rewritten",
                  length(rel_snps), length(setdiff(snps, rel_snps))))
  write_tsv(tibble(SNP = rel_snps), SNPS)
}
cat(sprintf("release: %d SNPs (%d in the overlap of the score files)\n",
            length(rel_snps), length(snps)))

if (length(allele_tally)) {
  for (src in names(allele_tally)) {
    cat(sprintf("%s: %d SNPs take the panel's allele from another row, %d scored for neither are dropped\n",
                src, allele_tally[[src]][["other_row"]], allele_tally[[src]][["none"]]))
  }
}

cat("\nnext: Rscript annotations/01-build-annotations.R from the repository root,\n",
    "then sbatch annotations/02-compute-ldscores.sh\n", sep = "")
