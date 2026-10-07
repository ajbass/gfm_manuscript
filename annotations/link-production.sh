#!/bin/bash
##################################################
### File links the paper's annotation files ###
## Into the layout the scripts read, so the
## analysis runs on the sets the paper used
## without rebuilding them or their LD scores
##################################################

# bash annotations/link-production.sh   (from the repository root)
# The production column of annotations/arms.tsv says where the paper's run
# kept each arm under REF_GWASAI (default: DATA_GWASAI itself); each of its
# annotation, l2.M and LD score files is linked as
#   foundational_models/<model>/<dir>/<model>.<chr>.<suffix>
# Existing files are left alone, so a rebuilt arm is never replaced
set -euo pipefail
DATA_GWASAI=${DATA_GWASAI:-/home/ab3105/rds/hpc-work/data_gwasAI}
REF_GWASAI=${REF_GWASAI:-$DATA_GWASAI}
ARMS_TSV=${ARMS_TSV:-annotations/arms.tsv}
FM=$DATA_GWASAI/foundational_models
RFM=$REF_GWASAI/foundational_models

n_link=0
n_have=0
n_miss=0
while read -r M D P; do
  mkdir -p "$FM/$M/$D"
  for chr in $(seq 1 22); do
    for sfx in annot.parquet l2.M l2.ldscore.parquet; do
      src=$RFM/$P$chr.$sfx
      dst=$FM/$M/$D/$M.$chr.$sfx
      if [[ -e $dst ]]; then
        n_have=$((n_have + 1))
      elif [[ -e $src ]]; then
        ln -s "$src" "$dst"
        n_link=$((n_link + 1))
      else
        n_miss=$((n_miss + 1))
        [[ $chr -eq 1 ]] && echo "not in the reference: $src"
      fi
    done
  done
done < <(awk -F'\t' 'NR > 1 { print $1, $3, $6 }' "$ARMS_TSV")
echo "linked $n_link files, $n_have already in place, $n_miss absent from $RFM"
