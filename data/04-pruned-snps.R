##################################################
### File collects the pruned SNPs as rsIDs ###
## The independent SNPs the sfFDR pi0 model is
## fitted on: qc_panukk/pruned_snps_0.3.tsv
##################################################
library(tidyverse)

qc <- read_tsv("./qc_panukk/full_variant_qc_metrics.txt.bgz",
               col_select = c("rsid", "varid"), show_col_types = FALSE)

files <- list.files("./prune_data_0.3", pattern = "\\.prune\\.in$",
                    full.names = TRUE)
if (length(files) != 22) {
  stop("expected 22 .prune.in files, found ", length(files))
}
pruned <- map_dfr(files, ~ read_tsv(.x, col_names = "varid",
                                    show_col_types = FALSE))

out <- pruned %>% left_join(qc, by = "varid") %>%
  filter(!is.na(rsid)) %>%
  distinct(rsid) %>%
  rename(SNP = rsid)
message(nrow(pruned), " pruned variants, ", nrow(out), " with an rsID")
write_tsv(out, "./qc_panukk/pruned_snps_0.3.tsv")
