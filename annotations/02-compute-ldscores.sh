#!/bin/bash
#SBATCH -J ld-combos
#SBATCH -A CWALLACE-SL2-CPU
#SBATCH --time=4:00:00
#SBATCH --mem=32000mb
#SBATCH --array=1-22%10
#SBATCH -p icelake-himem
#SBATCH -o logs/ldcombo_%A_%a.out
#SBATCH -e logs/ldcombo_%A_%a.err
##################################################
### File computes LD scores for every arm  ###
## One chromosome per task, all arms within
## Flags: --allow-missing, no window flag
##################################################

# Submit from the repository root. The arms come from annotations/arms.tsv,
# the paper's unless ARMS=all; existing outputs are skipped
chr=${SLURM_ARRAY_TASK_ID}
mkdir -p logs

module load conda
module load anaconda
source activate /home/ab3105/miniforge3/envs/polyfun

DATA_GWASAI=${DATA_GWASAI:-/home/ab3105/rds/hpc-work/data_gwasAI}
POLYFUN=${POLYFUN:-/rds/user/ab3105/hpc-work/polyfun}
ARMS_TSV=${ARMS_TSV:-annotations/arms.tsv}
ARMS=${ARMS:-paper}
FM=${DATA_GWASAI}/foundational_models
BF=${DATA_GWASAI}/baselineLD2.3/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}

# model and directory of each arm
awk -F'\t' -v want="$ARMS" 'NR > 1 && (want == "all" || $4 == 1) { print $1, $3 }' \
    "$ARMS_TSV" | while read -r M D; do
  ANNOT=${FM}/${M}/${D}/${M}.${chr}.annot.parquet
  OUT=${FM}/${M}/${D}/${M}.${chr}.l2.ldscore.parquet
  if [ -f "$OUT" ]; then
    echo "skip (exists): $OUT"
    continue
  fi
  if [ ! -f "$ANNOT" ]; then
    echo "skip (no annotations): $ANNOT"
    continue
  fi
  echo "=== ${M} (${D}) chr${chr} ==="
  python3 ${POLYFUN}/compute_ldscores.py \
    --allow-missing \
    --bfile ${BF} \
    --annot ${ANNOT} \
    --out   ${OUT}
  # --allow-missing lets panel SNPs without an annotation row through, and
  # silently drops an annotation row whose CHR:BP:A1:A2 is not in the
  # panel; every annotation row must come back with an LD score
  python3 - "$ANNOT" "$OUT" <<'PY'
import sys
import pyarrow.parquet as pq
a, o = (pq.ParquetFile(f).metadata.num_rows for f in sys.argv[1:3])
if a != o:
    print(f"!! LD scores hold {o} of {a} annotation rows: rows not matched to the panel were dropped")
PY
done

echo "Finished: `date`"
