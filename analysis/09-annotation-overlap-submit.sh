#!/bin/bash
#SBATCH -J annotation-overlap
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=24000mb
#SBATCH -p icelake
#SBATCH -o logs/overlap_%j.out
#SBATCH -e logs/overlap_%j.err

##################################################
### File submits the overlap measurement ###
## sbatch 09-annotation-overlap-submit.sh <model> [joint|subset]
## One arm carrying baselineLD per job; run-all.sh
## submits one for every baselineLD+<model> arm
##################################################
RBIN=${RBIN:-$HOME/R-4.4.2/bin}
MODEL=${1:?model name required}
VARIANT=${2:-joint}

. /etc/profile.d/modules.sh
mkdir -p logs
echo "Time: `date`   Node: `hostname`   Model: ${MODEL}   Variant: ${VARIANT}"
VARIANT=$VARIANT $RBIN/Rscript analysis/09-annotation-overlap.R ${MODEL}
echo "Finished: `date`"
