#!/bin/bash
##################################################
### File submits the whole pipeline in order ###
## Run from the repository root.
##################################################

#   bash run-all.sh               # the paper's arms, all but the h2 step
#   ARMS=all bash run-all.sh      # every arm run (annotations/arms.tsv)
#   RUN_H2=1 bash run-all.sh      # add S-LDSC h2
#   SKIP_BUILD=1 bash run-all.sh  # skip the annotation step
#
set -euo pipefail

# Cluster paths this analysis was run on. Export DATA_GWASAI (reference
# data) or ANDY (fit stores) to rebuild somewhere else
export DATA_GWASAI=${DATA_GWASAI:-/home/ab3105/rds/hpc-work/data_gwasAI}
export ANDY=${ANDY:-/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY}
export SUMSTATS=${SUMSTATS:-/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY}
export POLYFUN=${POLYFUN:-/rds/user/ab3105/hpc-work/polyfun}
export RBIN=${RBIN:-$HOME/R-4.4.2/bin}
# The arm table and which of its arms to run
export ARMS_TSV=${ARMS_TSV:-annotations/arms.tsv}
export ARMS=${ARMS:-paper}
DATA_DIR=$DATA_GWASAI
# The run tree is wherever this script lives
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_TREE=${SOURCE_TREE:-/home/ab3105/rds/hpc-work/gwasAI}
RUN_H2=${RUN_H2:-0}
SKIP_BUILD=${SKIP_BUILD:-0}
# Seconds between the per-model submissions (two per arm), to be gentle
# on the scheduler; SUBMIT_PAUSE=0 submits in one burst
SUBMIT_PAUSE=${SUBMIT_PAUSE:-2}

# sbatch flags for the small glue jobs
GLUE="-A CWALLACE-SL2-CPU -p icelake --time=04:00:00 --mem=32000mb -o logs/%x_%A.out -e logs/%x_%A.err"

cd "$WORK_DIR"
mkdir -p logs

# Model lists from the arm table
# baselineLD2.2 has no subset form: the joint fit is the reference for both
# variants
arm_models () {  # $1 = variant
  awk -F'\t' -v v="$1" -v want="$ARMS" \
      'NR > 1 && $2 == v && (want == "all" || $4 == 1) { print $1 }' "$ARMS_TSV"
}
JOINT_MODELS=$(arm_models joint)
SUBSET_MODELS=$(arm_models subset)
if [[ -z $JOINT_MODELS || -z $SUBSET_MODELS ]]; then
  echo "ERROR: no arms read from $ARMS_TSV (ARMS=$ARMS)" >&2
  exit 1
fi

# Run local inputs, linked from SOURCE_TREE
if [[ ! -e qc_panukk/pruned_snps_0.3.tsv ]]; then
  if [[ -e "$SOURCE_TREE/qc_panukk/pruned_snps_0.3.tsv" ]]; then
    ln -sfn "$SOURCE_TREE/qc_panukk" qc_panukk
    echo "linked qc_panukk -> $SOURCE_TREE/qc_panukk"
  else
    echo "ERROR: qc_panukk/pruned_snps_0.3.tsv not found here or in $SOURCE_TREE" >&2
    exit 1
  fi
fi
if [[ ! -x plink2 ]]; then
  if [[ -x "$SOURCE_TREE/plink2" ]]; then
    ln -sfn "$SOURCE_TREE/plink2" plink2
    echo "linked plink2 -> $SOURCE_TREE/plink2"
  else
    echo "ERROR: plink2 not found here or in $SOURCE_TREE" >&2
    exit 1
  fi
fi

# Annotations, on the login node and blocking: the arms are selections
# from the per-variant release, which 00-prepare-scores.R builds first
REL=$DATA_DIR/foundational_models/release
if [[ "$SKIP_BUILD" != "1" ]]; then
  if [[ ! -e $REL/baselineLD.1.parquet ]]; then
    echo "ERROR: no release under $REL" >&2
    echo "  run annotations/00-prepare-scores.R first (see annotations/README.md)" >&2
    exit 1
  fi
  echo "=== annotations: building the $ARMS arms from the release ==="
  "$RBIN"/Rscript annotations/01-build-annotations.R
else
  echo "=== annotations: SKIP_BUILD=1, assuming they are on disk ==="
fi

# A model with no annotation files would strand the chain
MISS=""
shopt -s nullglob
while read -r M V D; do
  f=("$DATA_DIR/foundational_models/$M/$D"/"$M".*.annot.parquet)
  [[ ${#f[@]} -lt 22 ]] && MISS="$MISS $V/$M(${#f[@]}/22)"
done < <(awk -F'\t' -v want="$ARMS" \
             'NR > 1 && (want == "all" || $4 == 1) { print $1, $2, $3 }' "$ARMS_TSV")
shopt -u nullglob
if [[ -n $MISS ]]; then
  echo "ERROR: annotation files missing for:$MISS" >&2
  echo "  build them first -- refusing to submit a doomed chain" >&2
  exit 1
fi

# LD scores
LD=$(sbatch --parsable annotations/02-compute-ldscores.sh)
echo "LD scores               -> job $LD"

# Optional S-LDSC heritability
if [[ "$RUN_H2" == "1" ]]; then
  H2=$(sbatch --parsable --dependency=afterok:$LD analysis/01-ldsc-h2-submit.sh)
  echo "S-LDSC h2 (optional)    -> job $H2"
fi

# Overlap of each model's LD scores with baselineLD, from the arms that
# carry baselineLD beside one model (the joint arm where there is one).
# Reads the LD scores only, so it runs beside the fits
overlap_arms () {
  awk -F'\t' -v want="$ARMS" \
      'NR > 1 && (want == "all" || $4 == 1) && $5 ~ /^baselineLD\+[^+]+$/ { print $1, $2 }' \
      "$ARMS_TSV" | sort -k1,1 -k2,2 | awk '!seen[$1]++'
}
while read -r M V; do
  O=$(sbatch --parsable --dependency=afterok:$LD analysis/09-annotation-overlap-submit.sh "$M" "$V")
  echo "overlap $V $M -> job $O"
done < <(overlap_arms)

# Ridge then sfFDR
SFF_IDS=""
RIDGE_IDS=""
submit_fit () {  # $1 = model, $2 = variant
  local RID SID
  RID=$(VARIANT=$2 sbatch --parsable --dependency=afterok:$LD analysis/02-ridge-echisq-submit.sh "$1")
  SID=$(VARIANT=$2 sbatch --parsable --dependency=afterok:$RID analysis/03-sffdr-submit.sh "$1")
  RIDGE_IDS="${RIDGE_IDS:+$RIDGE_IDS:}$RID"
  SFF_IDS="${SFF_IDS:+$SFF_IDS:}$SID"
  echo "ridge/sffdr $2 $1 -> jobs $RID -> $SID"
  sleep "$SUBMIT_PAUSE"
}
for M in $JOINT_MODELS;  do submit_fit "$M" joint;  done
for M in $SUBSET_MODELS; do submit_fit "$M" subset; done

# Clumping into independent loci
CL=$(sbatch --parsable --dependency=afterok:$SFF_IDS analysis/04-clump-submit.sh)
echo "clumping                -> job $CL"

# Locus counts, the discovery result
CO=$(sbatch --parsable $GLUE -J collect-clumps --dependency=afterok:$CL \
     --wrap "$RBIN/Rscript analysis/05-collect-clumps.R")
echo "locus counts            -> job $CO"

# Lead statistics. Shards are keyed by task number, so stale ones are
# cleared first
CLR8=$(sbatch --parsable $GLUE -J clear-leadshards --dependency=afterok:$CO \
      --wrap "rm -f indep-discoveries-1000G/_leadstats_shards/*.tsv")
L8=$(sbatch --parsable --dependency=afterok:$CLR8 analysis/06-lead-stats-submit.sh)
M8=$(sbatch --parsable $GLUE -J leadstats-merge --dependency=afterok:$L8 \
     --wrap "$RBIN/Rscript analysis/06-lead-stats.R merge")
echo "lead stats              -> jobs $CLR8 -> $L8 -> $M8"

# Shared / distinct classification
CLR9=$(sbatch --parsable $GLUE -J clear-collapseshards --dependency=afterok:$M8 \
      --wrap "rm -f indep-discoveries-1000G/_collapse_shards/*.tsv")
L9=$(sbatch --parsable --dependency=afterok:$CLR9 analysis/07-collapse-classify-submit.sh)
M9=$(sbatch --parsable $GLUE -J collapse-merge --dependency=afterok:$L9 \
     --wrap "$RBIN/Rscript analysis/07-collapse-classify.R merge")
echo "classification          -> jobs $CLR9 -> $L9 -> $M9"

# R2 metrics and annotation statistics
S10=$(sbatch --parsable --dependency=afterok:$RIDGE_IDS analysis/08-model-summary-submit.sh subset)
J10=$(sbatch --parsable --dependency=afterok:$RIDGE_IDS analysis/08-model-summary-submit.sh joint)
echo "model summary           -> jobs $S10 (subset), $J10 (joint)"

# Consistency of the stores and tables, last
V=$(sbatch --parsable $GLUE -J verify-run --dependency=afterok:$M9:$S10:$J10 \
    --wrap "$RBIN/Rscript verify-run.R")
echo "verify                  -> job $V (logs/verify-run_$V.out)"

echo
echo "All stages submitted. Monitor with:  squeue -u \$USER"
echo "When everything is COMPLETED, copy the result tables down and run the"
echo "figure scripts locally (see figures/README.md)."
