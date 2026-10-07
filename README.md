## Overview

This repository contains code to reproduce the results in

> Bass AJ, Park J, Mukherjee S, Whittaker JC, Wallace C. Genomic foundation models do not yet boost common variant discovery. *Submitted*; 2026.

Overview of directories:

- `data/`: UK Biobank inputs (samples, LD-pruned SNPs and summary statistics)
- `annotations/`: constructing the annotation sets and their LD scores
- `analysis/`: expected chi-squared fits, sfFDR, and model summaries
- `run-all.sh`: submits every cluster stage

The annotations are on Zenodo (https://doi.org/10.5281/zenodo.23063752).
The result tables of the paper's run are in a second record (https://doi.org/10.5281/zenodo.23110145). See subdirectories for additional details.

## Getting started

The scripts read paths from three environment variables: `DATA_GWASAI`
(reference data and annotations), `SUMSTATS` (summary statistics) and `ANDY`
(output). Set them first (they default to the original run's paths), and
`RBIN` and `POLYFUN` if R or PolyFun are elsewhere:

```bash
export DATA_GWASAI=/path/to/data    # annotations and reference data
export SUMSTATS=/path/to/sumstats   # summary statistics
export ANDY=/path/to/fits           # written by the pipeline
export RBIN=/path/to/R/bin          # R's bin directory (default $HOME/R-4.4.2/bin)
export POLYFUN=/path/to/polyfun     # the PolyFun checkout for the LD scores
```

The SLURM scripts use the original cluster's account (`-A CWALLACE-SL2-CPU`),
partitions and conda environment. Change these, and the `GLUE` line of
`run-all.sh`, to run elsewhere.

1. **Annotations** (https://doi.org/10.5281/zenodo.23063752): one tar per
   model, `columns.tsv` (the column guide), and the LD-pruned SNPs. From the
   repository root:

   ```bash
   mkdir -p $DATA_GWASAI/foundational_models/release qc_panukk
   for m in baselineLD borzoi alphagenome zoonomia gpn_star avi; do
     tar -xf $m.tar -C $DATA_GWASAI/foundational_models/release
   done
   cp columns.tsv $DATA_GWASAI/foundational_models/release/
   cp pruned_snps_0.3.tsv qc_panukk/
   # the analysis SNPs
   Rscript -e 'library(tidyverse); library(arrow)
     map_dfr(1:22, ~ read_parquet(sprintf("%s/foundational_models/release/baselineLD.%d.parquet",
       Sys.getenv("DATA_GWASAI"), .x), col_select = "SNP")) %>%
       write_tsv("qc_panukk/1000G_SNPs_overlap.tsv")'
   ```

2. **1000 Genomes EUR reference**: the Phase 3 PLINK files and HapMap3
   weights from the S-LDSC reference files (https://zenodo.org/records/8292725):
   `$DATA_GWASAI/baselineLD2.3/1000G_EUR_Phase3_plink/1000G.EUR.QC.<chr>.{bed,bim,fam}`
   and
   `$DATA_GWASAI/baselineLD2.3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.<chr>.l2.ldscore.gz`.
3. **UK Biobank summary statistics** (the Price lab BOLT-LMM release,
   https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/UKBB),
   one file per trait as
   `$SUMSTATS/UKBB_Price_37traits/indep_traits/<trait>.sumstats.gz`. Every
   file there is run. Munge them onto the analysis SNPs with
   `mkdir -p logs && sbatch data/05-munge-sumstats-submit.sh`. 
4. **PLINK 2** as `./plink2` in the repository root, and the software below.
5. **Run** from the repository root, once the munging jobs have finished:

   ```bash
   bash run-all.sh
   ```

   It builds the annotation sets (`annotations/01-build-annotations.R`), then
   submits their LD scores (`annotations/02-compute-ldscores.sh`), and every
   later stage with the dependencies chained.

Software: R 4.4.2 with `tidyverse` 2.0.0, `data.table` 1.16.4, `arrow`
20.0.0.2, `glmnet` 4.1-10, `glmnetUtils` 1.1.9, `digest` 0.6.37, `sffdr`
1.1.2, `topr` 2.0.2 (the versions the analysis was run with); Python 3.12 with
[PolyFun](https://github.com/omerwe/polyfun), and PLINK 2.0. 
