##################################################
### File measures overlap with baselineLD ###
##################################################
library(tidyverse)
library(arrow)

# How to run, from the repository root:
#   VARIANT=joint  Rscript analysis/09-annotation-overlap.R baselineLD2.2_borzoi
#   VARIANT=subset Rscript analysis/09-annotation-overlap.R baselineLD2.2_avi
# The arm must carry baselineLD beside the model's columns (the LD scores
# are the design matrix the ridge sees, so overlap is measured there).
# Writes, under annotation_overlap/,
#   <model>.columns.tsv     one row per model column, every statistic
#   <model>.summary.tsv     medians by statistic, baselineLD group, assay
#   <model>.cor.tsv         model column x baselineLD column correlations
#   <model>.corfull.tsv     the square matrix over both and base, which the
#                           partial correlations of the figure need
#   <model>.baselineld.tsv  each baselineLD column on the model's block

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) to read another tree.
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")
FM <- file.path(DATA_GWASAI, "foundational_models")
COLUMNS <- file.path(FM, "release", "columns.tsv")
OUT_DIR <- "annotation_overlap"
META <- c("CHR", "SNP", "BP", "CM", "A1", "A2", "MAF")

VARIANT <- Sys.getenv("VARIANT", "joint")
if (!VARIANT %in% c("subset", "joint")) {
  stop("VARIANT must be 'subset' or 'joint'")
}
ARMS_TSV <- Sys.getenv("ARMS_TSV", "annotations/arms.tsv")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: VARIANT=joint|subset ",
       "Rscript analysis/09-annotation-overlap.R <model>")
}
model_name <- args[1]
arm <- read_tsv(ARMS_TSV, show_col_types = FALSE) %>%
  filter(model == model_name, variant == VARIANT)
if (nrow(arm) != 1) {
  stop("no ", VARIANT, " arm named ", model_name, " in ", ARMS_TSV)
}
PREFIX <- file.path(FM, model_name, arm$dir, paste0(model_name, "."))
message("model : ", model_name, " (", VARIANT, ")")
message("files : ", PREFIX, "<chr>.l2.ldscore.parquet")

# The baselineLD classes, evaluated out of the build script so the two
# cannot drift. QTL and regulatory evolution stay apart
for (x in parse("annotations/01-build-annotations.R")) {
  if (is.call(x) && identical(x[[1]], as.name("<-")) &&
      startsWith(as.character(x[[2]]), "GRP_")) {
    eval(x)
  }
}
GROUPS <- list(base = "base",
               coding = GRP_CODING,
               conservation = GRP_CONSERVATION,
               regulatory = GRP_REGULATORY,
               qtl = GRP_QTL,
               regulatory_evolution = GRP_REGEVO,
               mafld = GRP_MAFLD)

## The columns: baselineLD from the release guide, the model's the rest
guide <- read_tsv(COLUMNS, show_col_types = FALSE)
bl_cols <- guide$column[guide$group == "baselineLD"]
ann <- names(open_dataset(paste0(PREFIX, 1L, ".l2.ldscore.parquet"))$schema) %>%
  setdiff(META)
if (!all(c("base", bl_cols) %in% ann)) {
  stop("the arm does not carry base and every baselineLD column")
}
model_cols <- setdiff(ann, c("base", bl_cols))
if (!length(model_cols)) {
  stop("the arm has no columns beyond baselineLD")
}
grp <- map(GROUPS, ~ intersect(.x, ann))
unassigned <- setdiff(bl_cols, unlist(grp))
if (length(unassigned)) {
  warning("baselineLD columns in no class: ",
          paste(unassigned, collapse = ", "))
}
message(sprintf("columns: %d baselineLD, %d %s", length(bl_cols),
                length(model_cols), model_name))

## Cross products over the genome, one chromosome in memory at a time
p <- length(ann)
XtX <- matrix(0, p, p, dimnames = list(ann, ann))
sums <- setNames(numeric(p), ann)
n <- 0
for (chr in 1:22) {
  f <- paste0(PREFIX, chr, ".l2.ldscore.parquet")
  if (!file.exists(f)) {
    stop("LD score file not found: ", f)
  }
  m <- read_parquet(f, col_select = ann) %>% as.matrix()
  storage.mode(m) <- "double"
  m <- m[complete.cases(m), , drop = FALSE]
  XtX <- XtX + crossprod(m)
  sums <- sums + colSums(m)
  n <- n + nrow(m)
  message(sprintf("chr %2d  %7d SNPs  (%d so far)", chr, nrow(m), n))
  rm(m)
}
mu <- sums / n
cv <- (XtX - n * outer(mu, mu)) / (n - 1)
sd <- sqrt(diag(cv))
keep <- is.finite(sd) & sd > 0
if (any(!keep)) {
  message("constant, dropped: ", paste(ann[!keep], collapse = ", "))
}
cr <- cv[keep, keep] / outer(sd[keep], sd[keep])
grp <- map(grp, ~ intersect(.x, colnames(cr)))
all_bl <- intersect(bl_cols, colnames(cr))
tg <- intersect(model_cols, colnames(cr))

# R2 of each target column on a block of predictors, from the correlation
# matrix; a small ridge keeps a collinear block invertible
block_r2 <- function(targets, preds, cr, lambda = 1e-4) {
  if (!length(preds)) {
    return(setNames(rep(NA_real_, length(targets)), targets))
  }
  inv <- solve(cr[preds, preds, drop = FALSE] + diag(lambda, length(preds)))
  r2 <- vapply(targets, function(t) {
    r <- cr[preds, t]
    return(min(max(as.numeric(t(r) %*% inv %*% r), 0), 1))
  }, numeric(1))
  return(r2)
}

## Per model column: marginal R2 per class, the increment each class adds
## over regulatory, the variance unique to each class among the six, and
## the marginal R2 with the total LD score (base) partialled out
nb <- setdiff(names(GROUPS), "base")
stats <- tibble(column = tg,
                assay = if_else(str_detect(tg, ":"),
                                str_extract(tg, "^[^:]+"), tg),
                r2_all = block_r2(tg, c("base", all_bl), cr),
                r2_functional = block_r2(tg, all_bl, cr))
for (g in names(GROUPS)) {
  stats[[paste0("marginal_", g)]] <- block_r2(tg, grp[[g]], cr)
}
r2_reg <- block_r2(tg, grp$regulatory, cr)
# base is a class of one, and the figure draws it beside the others, so it
# takes an increment over regulatory like every other class
for (g in setdiff(names(GROUPS), "regulatory")) {
  stats[[paste0("added_", g)]] <-
    block_r2(tg, union(grp$regulatory, grp[[g]]), cr) - r2_reg
}
for (g in nb) {
  stats[[paste0("unique_", g)]] <- stats$r2_functional -
    block_r2(tg, setdiff(all_bl, grp[[g]]), cr)
}
# every LD score carries the total LD score, which inflates every
# correlation. the partial correlations remove that shared part
rb <- cr[, "base"]
pc <- (cr - outer(rb, rb)) / sqrt(outer(1 - rb^2, 1 - rb^2))
diag(pc) <- 1
pc <- pc[colnames(cr) != "base", colnames(cr) != "base"]
for (g in nb) {
  stats[[paste0("partial_", g)]] <- block_r2(tg, grp[[g]], pc)
}

## Medians over the model's columns, for every assay and for all of them
long <- stats %>%
  pivot_longer(-c(column, assay), names_to = "statistic", values_to = "r2") %>%
  separate(statistic, c("statistic", "group"), sep = "_", extra = "merge")
summ <- bind_rows(long %>% mutate(assay = "all"), long) %>%
  group_by(statistic, group, assay) %>%
  summarise(n = n(), median_r2 = median(r2, na.rm = TRUE), .groups = "drop")

## The other direction, and the raw correlations
group_of <- function(col) {
  g <- names(grp)[map_lgl(grp, ~ col %in% .x)]
  return(if (length(g)) g[1] else NA_character_)
}
bl_on_model <- tibble(column = all_bl,
                      group = map_chr(all_bl, group_of),
                      r2_on_model = block_r2(all_bl, tg, cr)) %>%
  arrange(desc(r2_on_model))
cor_long <- as_tibble(cr[tg, all_bl, drop = FALSE], rownames = "column") %>%
  pivot_longer(-column, names_to = "baselineld", values_to = "r")

keep <- c(tg, all_bl, "base")
cor_full <- as_tibble(cr[keep, keep, drop = FALSE], rownames = "column")

dir.create(OUT_DIR, showWarnings = FALSE)
out <- function(what) file.path(OUT_DIR, paste0(model_name, ".", what, ".tsv"))
write_tsv(stats, out("columns"))
write_tsv(summ, out("summary"))
write_tsv(cor_long, out("cor"))
write_tsv(cor_full, out("corfull"))
write_tsv(bl_on_model, out("baselineld"))

cat(sprintf("\n%s (%s): %d columns, %d SNPs\n", model_name, VARIANT,
            length(tg), n))
cat(sprintf("median R2 on baselineLD: %.3f with base, %.3f without\n",
            median(stats$r2_all), median(stats$r2_functional)))
cat("median R2 by class (marginal / added over regulatory / unique /",
    "base partialled):\n")
for (g in nb) {
  cat(sprintf("  %-22s %.3f  %s  %.3f  %.3f\n", g,
              median(stats[[paste0("marginal_", g)]]),
              if (g == "regulatory") "  -  " else
                sprintf("%.3f", median(stats[[paste0("added_", g)]])),
              median(stats[[paste0("unique_", g)]]),
              median(stats[[paste0("partial_", g)]])))
}
cat(sprintf("baselineLD on the model's block: median R2 %.3f\n",
            median(bl_on_model$r2_on_model)))
cat("ALL DONE\n")
