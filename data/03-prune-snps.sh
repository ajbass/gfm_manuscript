#!/bin/bash
#SBATCH -J get-indp-snps
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-22%2
#SBATCH -p icelake-himem
#SBATCH -o logs/prune_%A_%a.out
#SBATCH -e logs/prune_%A_%a.err
##################################################
### File LD-prunes the UKB imputed genotypes ###
## r2 < 0.3 in 200 kb windows among the unrelated
## individuals, on the analysis SNPs
##################################################

chr=${SLURM_ARRAY_TASK_ID}
BIOBANK=${BIOBANK:-/home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank}
mkdir -p logs prune_data_0.3

# variant ids are set to CHR:BP_REF_ALT to match the Pan-UKB ids in
# qc_panukk/list_of_variant_ids.txt
./plink2 \
  --memory 1250000 \
  --bgen   ${BIOBANK}/genotypes-imputed/ukb22828_c${chr}_b0_v3.bgen ref-first \
  --sample ${BIOBANK}/genotypes-imputed/ukb22828_c${chr}_b0_v3_s487160.sample \
  --keep   ./unrelated_individuals.tab \
  --set-all-var-ids @:#_\$r_\$a --new-id-max-allele-len 113 missing \
  --extract ./qc_panukk/list_of_variant_ids.txt \
  --indep-pairwise 200kb 1 0.3 \
  --out ./prune_data_0.3/chr_prune_genotyped_200kb_0.3_${chr}
