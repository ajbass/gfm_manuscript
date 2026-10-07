#!/bin/bash
#SBATCH -J ldsc-h2
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=32000mb
#SBATCH --array=0-24%13
#SBATCH -p icelake-himem
#SBATCH -o logs/ldsch2_%A_%a.out
#SBATCH -e logs/ldsch2_%A_%a.err

##################################################
### File submits the S-LDSC heritability step ###
## One array task per trait shard; completed
## trait x model combinations are skipped
##################################################
RBIN=${RBIN:-$HOME/R-4.4.2/bin}
. /etc/profile.d/modules.sh
mkdir -p logs
module load conda
module load anaconda
source activate /home/ab3105/miniforge3/envs/polyfun

echo "Time: `date`   Node: `hostname`   Task: ${SLURM_ARRAY_TASK_ID}"
$RBIN/R -f analysis/01-ldsc-h2.R --args ${SLURM_ARRAY_TASK_ID}
echo "Finished: `date`"
