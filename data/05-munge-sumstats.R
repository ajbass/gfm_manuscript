##################################################
### File munges the UK Biobank summary stats ###
## Alleles matched to the 1000G EUR reference,
## INFO >= 0.8, chr6:24-45 Mb removed, on the
## analysis SNPs; then PolyFun's munge to parquet
##################################################
library(tidyverse)

SUMSTATS <- Sys.getenv("SUMSTATS",
                       "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")
POLYFUN <- Sys.getenv("POLYFUN", "/rds/user/ab3105/hpc-work/polyfun")
BIM <- file.path(DATA_GWASAI,
                 "baselineLD2.3/1000G_EUR_Phase3_plink/1000G.EUR.QC.%d.bim")
OUT_DIR <- file.path(SUMSTATS, "sumstats/UKBB-indep")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Array task: every 25th trait file
args <- commandArgs(trailingOnly = TRUE)
id2 <- as.numeric(args[1])
trait_files <- list.files(file.path(SUMSTATS,
                                    "UKBB_Price_37traits/indep_traits"),
                          full.names = TRUE)
trait_files <- trait_files[seq_along(trait_files) %% 25 == id2]

snps <- read_tsv("qc_panukk/1000G_SNPs_overlap.tsv", show_col_types = FALSE)
rsid <- map_dfr(1:22, ~ read_tsv(sprintf(BIM, .x),
                                 col_names = c("chr", "SNP", "CM", "pos", "alt",
                                               "ref"),
                                 show_col_types = FALSE)) %>%
  filter(SNP %in% snps$SNP) %>%
  select(-CM)

for (f in trait_files) {
  print(f)
  header <- colnames(read_tsv(f, n_max = 0, show_col_types = FALSE))
  cols <- which(header %in% c("SNP", "CHR", "POS", "A1", "A2", "REF", "EAF",
                              "Beta", "se", "P", "N", "INFO"))
  # P as text, then as.numeric: readr misparses values below ~1e-290
  # (1e-320 -> 1e-307); the floor at the smallest normal double follows
  data <- read_tsv(f, col_select = all_of(cols), show_col_types = FALSE,
                   col_types = cols(P = col_character(), .default = col_guess())) %>%
    mutate(P = suppressWarnings(as.numeric(P)))

  # Direct allele match
  tt0 <- rsid %>%
    inner_join(data %>% select(SNP, ref = A1, alt = A2, EAF, Beta, se, N, INFO,
                               P),
               by = c("SNP", "ref", "alt"))

  # Flipped allele match: beta and frequency follow
  tt1 <- rsid %>%
    inner_join(data %>% select(SNP, ref = A2, alt = A1, EAF, Beta, se, N, INFO,
                               P),
               by = c("SNP", "ref", "alt")) %>%
    mutate(Beta = -Beta, EAF = 1 - EAF)

  data <- bind_rows(tt0, tt1) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    filter(!is.na(P)) %>%
    mutate(EAF = pmin(EAF, 1 - EAF)) %>%
    select(SNP, chr, pos, ref, alt, Beta, se, P, EAF, INFO, N) %>%
    rename(A1 = alt, A2 = ref, BP = pos, CHR = chr,
           beta = Beta, MAF = EAF, pval = P) %>%
    filter(!(CHR == 6 & BP > 24e6 & BP < 45e6)) %>%
    filter(INFO >= 0.8)
  data$pval[data$pval < .Machine$double.xmin] <- .Machine$double.xmin

  fname <- str_remove(basename(f), "\\.gz$")
  out <- file.path(OUT_DIR, fname)
  write_tsv(data, out)
  system(paste("gzip -f", out))
  system(paste("python3", file.path(POLYFUN, "munge_polyfun_sumstats.py"),
               "--sumstats", paste0(out, ".gz"),
               "--out", paste0(out, ".parquet")))
  system(paste("rm", paste0(out, ".gz")))
}
