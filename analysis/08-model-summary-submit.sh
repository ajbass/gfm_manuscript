#!/bin/bash

#SBATCH -J model-summary
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --time=20:00:00
#SBATCH --mem=120000mb
#SBATCH --cpus-per-task=12
#SBATCH -p icelake-himem
#SBATCH -o logs/summary_%x_%A.out
#SBATCH -e logs/summary_%x_%A.err

##################################################
### File submits the cross-trait aggregation ###
## sbatch 08-model-summary-submit.sh          joint
## sbatch 08-model-summary-submit.sh subset   subset
## The baseline reference is echisq/baselineLD2.2
##################################################

RBIN=${RBIN:-$HOME/R-4.4.2/bin}
. /etc/profile.d/modules.sh

module load anaconda
source activate /home/ab3105/miniforge3/envs/polyfun

mkdir -p logs

VARIANT=${1:-joint}
if [[ "$VARIANT" != "joint" && "$VARIANT" != "subset" ]]; then
  echo "Error: variant must be 'joint' or 'subset', got '$VARIANT'"
  exit 1
fi

echo "Variant: $VARIANT"
echo "Started: $(date)   Node: $(hostname)"

VARIANT=$VARIANT $RBIN/Rscript analysis/08-model-summary.R

echo "Finished: $(date)"
