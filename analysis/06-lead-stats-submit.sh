#!/bin/bash
#SBATCH -J lead-stats
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=24000mb
#SBATCH --array=0-29%15
#SBATCH -p icelake
#SBATCH -o logs/leadstats_%A_%a.out
#SBATCH -e logs/leadstats_%A_%a.err

##################################################
### File submits the lead-statistics array ###
## Each task writes a shard; combine them with
## Rscript 06-lead-stats.R merge when all are done
##################################################

RBIN=${RBIN:-$HOME/R-4.4.2/bin}
. /etc/profile.d/modules.sh
mkdir -p logs
echo "Time: `date`   Task: ${SLURM_ARRAY_TASK_ID}"
$RBIN/R -f analysis/06-lead-stats.R --args ${SLURM_ARRAY_TASK_ID} 30
echo "Finished: `date`"
