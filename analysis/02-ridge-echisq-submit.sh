#!/bin/bash
#SBATCH -J ridge-echisq
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=96000mb
#SBATCH --array=0-12%13
#SBATCH -p icelake-himem
#SBATCH -o logs/ridge_%A_%a.out
#SBATCH -e logs/ridge_%A_%a.err

##################################################
### File submits the echisq regression step ###
## sbatch 02-ridge-echisq-submit.sh <model>
## VARIANT=joint sbatch ... <model> for the joint
## fit (baselineLD2.2 and baselineLD2.2_<gfm>)
##################################################
# 13 shards, matching the %% 13 in the R script
RBIN=${RBIN:-$HOME/R-4.4.2/bin}
MODEL=${1:?"usage: sbatch 02-ridge-echisq-submit.sh <model_name>"}

. /etc/profile.d/modules.sh
mkdir -p logs
echo "Time: `date`   Node: `hostname`   Task: ${SLURM_ARRAY_TASK_ID}   Model: ${MODEL}   Variant: ${VARIANT:-subset}"
$RBIN/R -f analysis/02-ridge-echisq.R --args ${SLURM_ARRAY_TASK_ID} ${MODEL}
echo "Finished: `date`"
