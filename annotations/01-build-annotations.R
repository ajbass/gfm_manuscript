##################################################
### File builds the annotation sets ###
##################################################
library(tidyverse)
library(arrow)

# The arm table: model, variant, dir, paper, groups
ARMS_TSV <- Sys.getenv("ARMS_TSV", "annotations/arms.tsv")
ARM_SET <- Sys.getenv("ARMS", "paper")
if (!ARM_SET %in% c("paper", "all")) {
  stop("ARMS must be 'paper' or 'all'")
}
arms <- read_tsv(ARMS_TSV, show_col_types = FALSE)
if (ARM_SET == "paper") {
  arms <- arms %>% filter(paper == 1)
}

# Cluster paths this analysis was run on.
setwd(Sys.getenv("DATA_GWASAI", "/home/ab3105/rds/hpc-work/data_gwasAI"))
REL <- "./foundational_models/release"
REL_FILE <- file.path(REL, "%s.%d.parquet")
COLUMNS <- file.path(REL, "columns.tsv")
OUT <- "./foundational_models/%s/%s/%s.%d.%s"

KEYS <- c("SNP", "CHR", "BP", "A1", "A2")

# The baselineLD classes
GRP_MAFLD <- c("MAFbin1", "MAFbin2", "MAFbin3", "MAFbin4", "MAFbin5",
               "MAFbin6", "MAFbin7", "MAFbin8", "MAFbin9", "MAFbin10",
               "MAF_Adj_Predicted_Allele_Age", "MAF_Adj_LLD_AFR",
               "Recomb_Rate_10kb", "Nucleotide_Diversity_10kb",
               "Backgrd_Selection_Stat", "CpG_Content_50kb", "MAF_Adj_ASMC")

GRP_CODING <- c("Coding_UCSC", "Coding_UCSC.flanking.500",
                "Intron_UCSC", "Intron_UCSC.flanking.500",
                "UTR_3_UCSC", "UTR_3_UCSC.flanking.500",
                "UTR_5_UCSC", "UTR_5_UCSC.flanking.500",
                "synonymous", "non_synonymous")

GRP_CONSERVATION <- c("Conserved_LindbladToh",
                      "Conserved_LindbladToh.flanking.500",
                      "GERP.NS", "GERP.RSsup4",
                      "Conserved_Vertebrate_phastCons46way",
                      "Conserved_Vertebrate_phastCons46way.flanking.500",
                      "Conserved_Mammal_phastCons46way",
                      "Conserved_Mammal_phastCons46way.flanking.500",
                      "Conserved_Primate_phastCons46way",
                      "Conserved_Primate_phastCons46way.flanking.500")

GRP_REGULATORY <- c("CTCF_Hoffman", "CTCF_Hoffman.flanking.500",
                    "DGF_ENCODE", "DGF_ENCODE.flanking.500",
                    "DHS_peaks_Trynka", "DHS_Trynka", "DHS_Trynka.flanking.500",
                    "Enhancer_Andersson", "Enhancer_Andersson.flanking.500",
                    "Enhancer_Hoffman", "Enhancer_Hoffman.flanking.500",
                    "FetalDHS_Trynka", "FetalDHS_Trynka.flanking.500",
                    "H3K27ac_Hnisz", "H3K27ac_Hnisz.flanking.500",
                    "H3K27ac_PGC2", "H3K27ac_PGC2.flanking.500",
                    "H3K4me1_peaks_Trynka", "H3K4me1_Trynka",
                    "H3K4me1_Trynka.flanking.500",
                    "H3K4me3_peaks_Trynka", "H3K4me3_Trynka",
                    "H3K4me3_Trynka.flanking.500",
                    "H3K9ac_peaks_Trynka", "H3K9ac_Trynka",
                    "H3K9ac_Trynka.flanking.500",
                    "PromoterFlanking_Hoffman",
                    "PromoterFlanking_Hoffman.flanking.500",
                    "Promoter_UCSC", "Promoter_UCSC.flanking.500",
                    "Repressed_Hoffman", "Repressed_Hoffman.flanking.500",
                    "SuperEnhancer_Hnisz", "SuperEnhancer_Hnisz.flanking.500",
                    "TFBS_ENCODE", "TFBS_ENCODE.flanking.500",
                    "Transcr_Hoffman", "Transcr_Hoffman.flanking.500",
                    "TSS_Hoffman", "TSS_Hoffman.flanking.500",
                    "WeakEnhancer_Hoffman", "WeakEnhancer_Hoffman.flanking.500",
                    "BivFlnk", "BivFlnk.flanking.500")

GRP_QTL <- c("GTEx_eQTL_MaxCPP", "BLUEPRINT_H3K27acQTL_MaxCPP",
             "BLUEPRINT_H3K4me1QTL_MaxCPP", "BLUEPRINT_DNA_methylation_MaxCPP")

GRP_REGEVO <- c("Human_Promoter_Villar", "Human_Promoter_Villar.flanking.500",
                "Human_Enhancer_Villar", "Human_Enhancer_Villar.flanking.500",
                "Ancient_Sequence_Age_Human_Promoter",
                "Ancient_Sequence_Age_Human_Promoter.flanking.500",
                "Ancient_Sequence_Age_Human_Enhancer",
                "Ancient_Sequence_Age_Human_Enhancer.flanking.500",
                "Human_Enhancer_Villar_Species_Enhancer_Count",
                "Human_Promoter_Villar_ExAC",
                "Human_Promoter_Villar_ExAC.flanking.500")

## The columns each group names
if (!file.exists(COLUMNS) ||
    !file.exists(sprintf(REL_FILE, "baselineLD", 1L))) {
  stop("no release under ", REL, ": run annotations/00-prepare-scores.R first")
}
guide <- read_tsv(COLUMNS, show_col_types = FALSE)
bl_ann <- guide$column[guide$group == "baselineLD"]
CLASSES <- list(mafld = intersect(GRP_MAFLD, bl_ann),
                regulatory = intersect(GRP_REGULATORY, bl_ann),
                conservation = intersect(GRP_CONSERVATION, bl_ann),
                coding = intersect(GRP_CODING, bl_ann),
                qtlregevo = intersect(c(GRP_QTL, GRP_REGEVO), bl_ann),
                functional = setdiff(bl_ann, GRP_MAFLD))
n_class <- sum(lengths(CLASSES[names(CLASSES) != "functional"]))
message(sprintf("baselineLD classes: %s (%d of %d)",
                paste(sprintf("%s %d", names(CLASSES), lengths(CLASSES)),
                      collapse = ", "),
                n_class, length(bl_ann)))
if (n_class != length(bl_ann)) {
  warning("class membership does not partition baselineLD: ", n_class,
          " assigned of ", length(bl_ann))
}
group_cols <- function(g) {
  if (g %in% names(CLASSES)) {
    return(CLASSES[[g]])
  }
  return(guide$column[guide$group == g])
}

# Each arm's columns
arms <- arms %>%
  mutate(cols = map(str_split(groups, fixed("+")),
                    ~ unlist(map(.x, group_cols))),
         files = map(cols, ~ unique(guide$file[match(.x, guide$column)])),
         missing = map_chr(str_split(groups, fixed("+")),
                           ~ paste(.x[lengths(map(.x, group_cols)) == 0],
                                   collapse = ", ")))
for (i in which(arms$missing != "")) {
  message(sprintf("  %s/%s SKIPPED: no columns for %s", arms$variant[i],
                  arms$model[i], arms$missing[i]))
}
arms <- arms %>% filter(missing == "")
dup <- arms$model[map_lgl(arms$cols, ~ anyDuplicated(.x) > 0)]
if (length(dup)) {
  stop("arms with a column twice: ", paste(dup, collapse = ", "))
}
arm_file <- function(a, chr, sfx) {
  return(sprintf(OUT, a$model, a$dir, a$model, chr, sfx))
}
arms <- arms %>%
  mutate(done = map_lgl(seq_len(n()), function(i) {
    all(file.exists(arm_file(arms[i, ], 1:22, "annot.parquet")))
  }))
for (i in which(arms$done)) {
  message(sprintf("  %s/%s exists", arms$variant[i], arms$model[i]))
}
arms <- arms %>% filter(!done)
message(sprintf("building %d arms (%s set)", nrow(arms), ARM_SET))
walk(unique(sprintf("./foundational_models/%s/%s", arms$model, arms$dir)),
     dir.create, showWarnings = FALSE, recursive = TRUE)

files <- unique(unlist(arms$files))
for (chr in 1:22) {
  if (!nrow(arms)) {
    break
  }
  d <- NULL
  for (f in files) {
    r <- read_parquet(sprintf(REL_FILE, f, chr))
    if (!identical(setdiff(names(r), c(KEYS, "base")),
                   guide$column[guide$file == f])) {
      stop("chr", chr, " ", f, ": release columns do not match columns.tsv")
    }
    d <- if (is.null(d)) {
      r
    } else {
      inner_join(d, r %>% select(-all_of(KEYS[-1]), -base), by = "SNP")
    }
    if (nrow(d) != nrow(r)) {
      stop("chr", chr, " ", f, ": release files do not share one SNP set")
    }
    rm(r)
  }
  for (i in seq_len(nrow(arms))) {
    a <- arms[i, ]
    x <- d %>% select(all_of(KEYS), base, all_of(a$cols[[1]]))
    write_parquet(x, arm_file(a, chr, "annot.parquet"))
    # .l2.M: column sums, base included, no header
    m <- x %>% select(-all_of(KEYS)) %>% summarise(across(everything(), sum))
    write_tsv(m, arm_file(a, chr, "l2.M"), col_names = FALSE)
    if (chr == 1) {
      message(sprintf("  %-30s %-6s %d SNPs, base + %d annotations", a$model,
                      a$variant, nrow(x), ncol(x) - length(KEYS) - 1L))
    }
  }
  message("chr", chr, " done")
}

cat("\nnext: sbatch annotations/02-compute-ldscores.sh from the repository root\n",
    "(same arms; existing LD scores are skipped), then analysis/ in order or\n",
    "run-all.sh for the whole chain\n", sep = "")
