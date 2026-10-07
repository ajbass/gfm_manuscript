#!/bin/bash
#SBATCH -J sffdr
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=48000mb
#SBATCH --array=0-12%13
#SBATCH -p icelake-himem
#SBATCH -o logs/sffdr_%A_%a.out
#SBATCH -e logs/sffdr_%A_%a.err

##################################################
### File submits the sfFDR step ###
## sbatch 03-sffdr-submit.sh            all models
## sbatch 03-sffdr-submit.sh <model>    one model
## VARIANT=joint sbatch ...             joint fits
##################################################
RBIN=${RBIN:-$HOME/R-4.4.2/bin}
MODEL=${1:-}

. /etc/profile.d/modules.sh
mkdir -p logs
echo "Time: `date`   Node: `hostname`   Task: ${SLURM_ARRAY_TASK_ID}   Model: ${MODEL:-ALL}   Variant: ${VARIANT:-subset}"
$RBIN/R -f analysis/03-sffdr.R --args ${SLURM_ARRAY_TASK_ID} ${MODEL} || exit 1
echo "Finished: `date`"
