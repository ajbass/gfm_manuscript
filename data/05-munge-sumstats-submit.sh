#!/bin/bash
#SBATCH -J munge-sumstats
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --time=04:00:00
#SBATCH --mem=32000mb
#SBATCH --array=0-24
#SBATCH -p icelake
#SBATCH -o logs/munge_%A_%a.out
#SBATCH -e logs/munge_%A_%a.err

##################################################
### File submits the summary statistic munging ###
## One array task per 25th trait file; submit
## from the repository root
##################################################
module load conda
module load anaconda
source activate /home/ab3105/miniforge3/envs/polyfun

RBIN=${RBIN:-$HOME/R-4.4.2/bin}
$RBIN/R -f data/05-munge-sumstats.R --args ${SLURM_ARRAY_TASK_ID}
