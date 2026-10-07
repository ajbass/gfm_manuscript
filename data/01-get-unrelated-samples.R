##################################################
### File lists the UK Biobank samples used ###
## Individuals used in the UKB genetic PCA (field
## 22020: unrelated, all ancestries), minus
## withdrawals
##################################################
library(tidyverse)
library(data.table)

# Cluster paths this analysis was run on. Set BIOBANK to the UK Biobank
# release directory elsewhere
BIOBANK <- Sys.getenv("BIOBANK",
                      "/home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank")
pheno_file <- file.path(BIOBANK, "phenotypes/ukb672224.tab")
withdrawal_file <- file.path(BIOBANK, "phenotypes/withdraw98032_19.txt")

df <- fread(pheno_file, select = c("f.eid", "f.22020.0.0"), sep = "\t")
setDF(df)
colnames(df) <- c("IID", "used.in.pca")
df <- df[!is.na(df$used.in.pca), ]

withdrawals <- read_tsv(withdrawal_file, col_names = FALSE) %>%
  rename(IID = X1)

out <- data.frame(IID = df$IID, FID = df$IID) %>%
  filter(!(IID %in% withdrawals$IID))

write_tsv(out[, c("FID", "IID")], file = "./unrelated_individuals.tab")
