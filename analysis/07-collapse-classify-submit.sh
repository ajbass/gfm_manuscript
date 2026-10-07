#!/bin/bash
#SBATCH -J collapse-cls
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16000mb
#SBATCH --array=0-13%14
#SBATCH -p icelake
#SBATCH -o logs/collapse_%A_%a.out
#SBATCH -e logs/collapse_%A_%a.err

##################################################
### File submits the locus classification array ###
## One task per group of traits; combine with
## Rscript 07-collapse-classify.R merge afterward
##################################################

RBIN=${RBIN:-$HOME/R-4.4.2/bin}
. /etc/profile.d/modules.sh
mkdir -p logs
echo "Time: `date`   Task: ${SLURM_ARRAY_TASK_ID}"
$RBIN/R -f analysis/07-collapse-classify.R --args ${SLURM_ARRAY_TASK_ID} 14
echo "Finished: `date`"
