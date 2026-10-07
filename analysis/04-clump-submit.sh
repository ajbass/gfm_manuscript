#!/bin/bash
#SBATCH -J clump-1000G
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=16000mb
#SBATCH --array=0-29%15
#SBATCH -p icelake
#SBATCH -o logs/clump1000G_%A_%a.out
#SBATCH -e logs/clump1000G_%A_%a.err

##################################################
### File submits the 1000G clumping array ###
## 1000G EUR is 489 samples, so memory and
## walltime are modest
##################################################

RBIN=${RBIN:-$HOME/R-4.4.2/bin}
. /etc/profile.d/modules.sh
mkdir -p logs

echo "Time: `date`"
echo "Node: `hostname`   Task: ${SLURM_ARRAY_TASK_ID}   Dir: `pwd`"

$RBIN/R -f analysis/04-clump.R --args ${SLURM_ARRAY_TASK_ID} 30 || exit 1

echo "Finished: `date`"
