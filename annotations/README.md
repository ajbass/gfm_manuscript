# Code to build the annotation sets

Overview of files:

- `aa_vcf.ipynb`: the notebook used to get the AlphaGenome scores
- `00-prepare-scores.R`: builds the per-variant annotation release from the downloads and the AlphaGenome scores
- `01-build-annotations.R`: writes every annotation set (arm) as a selection of the release's columns
- `02-compute-ldscores.sh`: computes LD scores for each set against the 1000G EUR reference
- `03-package-release.R`: bundles the release for distribution (one tar per model, `columns.tsv`, the LD-pruned SNP list)
- `arms.tsv`: the table of arms every script reads
- `link-inputs.sh`: links the downloads into a fresh `DATA_GWASAI`
- `link-production.sh`: links the files of the original run into the layout the scripts read

The scores come from:

- baselineLD v2.2 and v2.3 (with Zoonomia): https://zenodo.org/records/8292725
- Borzoi: https://storage.googleapis.com/seqnn-share/sniff/borzoi_102_annotation_set
- AlphaGenome: the API (https://www.alphagenomedocs.com/)
- AVI: the AlphaGenome Atlas (https://deepmind.google.com/science/alphagenome)
- GPN-Star: https://huggingface.co/datasets/songlab/ldsc

The output is available for download on Zenodo (https://doi.org/10.5281/zenodo.23063752).
