#!/bin/bash
##################################################
### File links the inputs into a fresh tree ###
## So a rebuild writes beside nothing it reads:
## the downloads and reference data are linked
## from REF_GWASAI, everything derived is built
## anew under DATA_GWASAI
##################################################

# REF_GWASAI=<tree the paper was run on> DATA_GWASAI=<new tree> \
#   bash annotations/link-inputs.sh
# Links, under DATA_GWASAI: baselineLD2.3/ (the 1000G plink files, HM3
# weights and 37->38 map, all read-only), the AlphaGenome and GPN-Star
# score files, the Borzoi list and raw components, the baselineLD v2.2 and
# v2.3 annotation downloads and the Atlas AVI file. Nothing already present
# is replaced
set -euo pipefail
DATA_GWASAI=${DATA_GWASAI:?set DATA_GWASAI to the new tree}
REF_GWASAI=${REF_GWASAI:?set REF_GWASAI to the tree the paper was run on}
if [[ $(readlink -f "$DATA_GWASAI") == $(readlink -f "$REF_GWASAI") ]]; then
  echo "ERROR: DATA_GWASAI and REF_GWASAI are the same tree" >&2
  exit 1
fi
FM=foundational_models
BZ=$FM/borzoi/borzoi_102_annotation_set

link () {  # $1 = path relative to both trees
  local src=$REF_GWASAI/$1 dst=$DATA_GWASAI/$1
  if [[ -e $dst ]]; then
    echo "in place: $1"
  elif [[ -e $src ]]; then
    mkdir -p "$(dirname "$dst")"
    ln -s "$src" "$dst"
    echo "linked:   $1"
  else
    echo "MISSING:  $src"
  fi
}

link baselineLD2.3
# the API export, or a gzipped copy of it
if [[ -e $REF_GWASAI/$FM/alphagenome/scores.csv ]]; then
  link $FM/alphagenome/scores.csv
else
  link $FM/alphagenome/scores.csv.gz
fi
link $FM/gpn-star/scores.parquet
link $BZ/SNP_LIST_1000G.tsv
link $BZ/annot_raw
for chr in $(seq 1 22); do
  link $FM/baselineLD2.2/annot/baselineLD.$chr.annot.gz
  link $FM/zoonomia/annot/baselineLD.$chr.annot.gz
done
# the Atlas file and its index, not the directory, so the extraction cache
# is written here
for f in "$REF_GWASAI"/$FM/avi/raw/*.gz "$REF_GWASAI"/$FM/avi/raw/*.bgz \
         "$REF_GWASAI"/$FM/avi/raw/*.tbi; do
  [[ -e $f ]] && link "$FM/avi/raw/$(basename "$f")"
done
echo "done: build with DATA_GWASAI=$DATA_GWASAI"
