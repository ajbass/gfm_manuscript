# Code to prepare the UK Biobank inputs

Overview of files:

- `01-get-unrelated-samples.R`: lists the individuals in the UK Biobank PCA
- `02-variant-ids.R`: maps the analysis SNPs to Pan-UKB variant ids and copies the SNP list into `qc_panukk/`
- `03-prune-snps.sh`: LD-prunes the imputed genotypes (r2 < 0.3, 200 kb) on those SNPs
- `04-pruned-snps.R`: lists the pruned variants for sfFDR
- `05-munge-sumstats.R`: munges the summary statistics onto the analysis SNPs and the 1000G alleles (INFO >= 0.8, chr6:24-45 Mb removed)
