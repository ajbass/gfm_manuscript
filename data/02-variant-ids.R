##################################################
### File maps the analysis SNPs to UKB ids ###
## The rsIDs scored by every model
##################################################
library(tidyverse)
library(data.table)

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) to build somewhere else
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")

# Pan-UKB variant QC file
QC <- "./qc_panukk/full_variant_qc_metrics.txt.bgz"
qc_cols <- match(c("rsid", "varid"),
                 strsplit(readLines(gzfile(QC), n = 1), "\t")[[1]])
if (anyNA(qc_cols)) {
  stop("rsid and varid columns not found in ", QC)
}
qc <- fread(cmd = sprintf("zcat %s | cut -f %d,%d", QC, qc_cols[1], qc_cols[2]),
            col.names = c("rsid", "varid"), sep = "\t", showProgress = FALSE) %>%
  as_tibble()

snps <- read_tsv(file.path(DATA_GWASAI, "foundational_models/release/snps.tsv"),
                 show_col_types = FALSE)
write_tsv(snps, "./qc_panukk/1000G_SNPs_overlap.tsv")

ids <- snps %>% rename(rsid = SNP) %>%
  left_join(qc, by = "rsid") %>%
  filter(!is.na(varid)) %>%
  select(varid)
write_tsv(ids, "./qc_panukk/list_of_variant_ids.txt", col_names = FALSE)
