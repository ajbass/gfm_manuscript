##################################################
### File computes the expected chi-squared ###
##################################################
library(tidyverse)
library(arrow)
library(topr)
library(sffdr)
library(digest)
library(glmnet)
library(glmnetUtils)

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) or ANDY (fit stores) to rebuild somewhere else.
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")
ANDY <- Sys.getenv("ANDY", "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
# inputs, read only
SUMSTATS <- Sys.getenv("SUMSTATS",
                       "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
FM <- file.path(DATA_GWASAI, "foundational_models")

# VARIANT switches input and output together:
#   subset = sets without full baselineLD
#   joint  = full joint LD scores; baselineLD2.2
#            and baselineLD2.2_<gfm> only
VARIANT <- Sys.getenv("VARIANT", "subset")
if (!VARIANT %in% c("subset", "joint")) {
  stop("VARIANT must be 'subset' or 'joint'")
}
OUT_BASE <- if (VARIANT == "subset") {
  file.path(ANDY, "sffdr/echisq_subset")
} else {
  file.path(ANDY, "sffdr/echisq")
}

# Model -> LD score prefix, from the arm table for this variant; only the
# baselineLD2.2 and baselineLD2.2_<gfm> models have a joint form, and
# mafld_only is the architecture control
ARMS_TSV <- Sys.getenv("ARMS_TSV", "annotations/arms.tsv")
arms <- read_tsv(ARMS_TSV, show_col_types = FALSE) %>%
  filter(variant == VARIANT)
MODEL_PREFIX <- setNames(file.path(FM, arms$model, arms$dir,
                                   paste0(arms$model, ".")),
                         arms$model)

# Arguments: array id and one model name
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript <script> <array_id> <model_name>\n  model_name one of:\n    ",
       paste(names(MODEL_PREFIX), collapse = "\n    "))
}
id2 <- as.numeric(args[1])
model_name <- args[2]
if (!model_name %in% names(MODEL_PREFIX)) {
  stop("Unknown model name: '", model_name, "'\n  allowed: ",
       paste(names(MODEL_PREFIX), collapse = ", "))
}

message("Running covariate-adjusted regression for model: ", model_name)
message("Variant: ", VARIANT, "   writing under: ", OUT_BASE)
message("Job ID index: ", id2)

# LD score merge, allowing an allele flip
get_model_path <- function(model, chr)
  paste0(MODEL_PREFIX[[model]], chr, ".l2.ldscore.parquet")

get_ldscores <- function(sumstats, model_name) {
  target_snps <- sumstats$SNP
  ld_list <- lapply(1:22, function(chr) {
    path <- get_model_path(model_name, chr)
    if (!file.exists(path)) {
      stop(paste("LD score file not found:", path))
    }
    arrow::open_dataset(path) %>% filter(SNP %in% target_snps) %>% collect()
  })
  comb <- bind_rows(ld_list)
  p1 <- sumstats %>%
    inner_join(comb, by = c("SNP", "BP", "A1", "A2", "CHR")) %>%
    filter(!is.na(base))
  # a flipped match relabels the alleles, so the effect flips with them
  p2 <- sumstats %>% mutate(A1_tmp = A1, A1 = A2, A2 = A1_tmp) %>%
    select(-A1_tmp) %>%
    inner_join(comb, by = c("SNP", "BP", "A1", "A2", "CHR")) %>%
    filter(!is.na(base)) %>%
    mutate(Beta = -Beta, EAF = 1 - EAF)
  bind_rows(p1, p2) %>% distinct(SNP, .keep_all = TRUE)
}

# fail now, not once per trait
miss <- (1:22)[!file.exists(vapply(1:22,
                                   function(c) get_model_path(model_name, c),
                                   ""))]
if (length(miss)) {
  stop("LD scores missing for ", model_name, ", chr: ",
       paste(miss, collapse = ","),
       "\n  compute LD scores for this annotation set first")
}
message("preflight OK: LD scores present for all 22 chromosomes")

# Traits come from sumstats/, not S-LDSC
sumstats_dir <- file.path(SUMSTATS, "sumstats/UKBB-indep")
if (!dir.exists(sumstats_dir)) {
  stop("sumstats directory not found: ", sumstats_dir)
}
all_files <- list.files(sumstats_dir, pattern = "\\.sumstats\\.parquet$",
                        full.names = TRUE)
if (!length(all_files)) {
  stop("no .sumstats.parquet files in ", sumstats_dir)
}
id <- (seq_along(all_files) %% 13) == id2
filename_subset <- all_files[id]
message("Processing ", length(filename_subset), " traits out of ",
        length(all_files))

out_dir <- file.path(OUT_BASE, model_name)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message("writing to: ", out_dir)

# HM3 regression weights
weights_path <- file.path(DATA_GWASAI, "baselineLD2.3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.")
hm3_list <- lapply(1:22, function(i) {
  read_tsv(paste0(weights_path, i, ".l2.ldscore.gz"), col_types = cols(),
           lazy = FALSE) %>%
    select(SNP, L2)
})
hm3 <- bind_rows(hm3_list)

exact_p <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[!is.na(x) & x < .Machine$double.xmin] <- .Machine$double.xmin
  return(x)
}
p_cols <- function(f, which = c("P", "fP", "fpvalues")) {
  hdr <- names(read_tsv(f, n_max = 0, show_col_types = FALSE))
  want <- intersect(which, hdr)
  spec <- replicate(length(want), col_character(), simplify = FALSE)
  names(spec) <- want
  return(do.call(cols, c(spec, list(.default = col_guess()))))
}

# Fit one trait at a time
for (file in filename_subset) {
  message("Processing: ", file)
  tryCatch({
    sid <- str_remove(basename(file), "\\.sumstats\\.parquet$")

    # resume: skip completed traits
    if (file.exists(file.path(out_dir, paste0(sid, ".sffdr.gz"))) &&
        file.exists(file.path(out_dir, paste0(sid, ".coefs.gz")))) {
      message("  skip (already done): ", sid)
      next
    }

    set.seed(readBin(digest(sid, raw = TRUE), "integer"))
    print(sid)

    sumstats <- arrow::read_parquet(file) # `file` IS the sumstats parquet

    meta_path <- file.path(SUMSTATS, "UKBB_Price_37traits/indep_traits",
                           paste0(sid, ".sumstats.gz"))
    if (!file.exists(meta_path)) {
      message("Meta missing: ", meta_path)
      next
    }
    getP <- read_tsv(meta_path, col_types = p_cols(meta_path, "P"), lazy = FALSE) %>%
      mutate(P = exact_p(P))

    # allele matching, allowing a flip
    ref_data <- getP %>% select(SNP, POS, A1, A2, CHR, EAF, Beta, se, P) %>%
      rename(BP = POS, A1_ref = A1, A2_ref = A2)
    match_direct <- sumstats %>%
      inner_join(ref_data, by = c("SNP", "CHR", "BP")) %>%
      filter(A1 == A1_ref & A2 == A2_ref)
    match_flipped <- sumstats %>%
      inner_join(ref_data, by = c("SNP", "CHR", "BP")) %>%
      filter(A1 == A2_ref & A2 == A1_ref) %>%
      mutate(Beta = -Beta, EAF = 1 - EAF) # EAF is the A1 frequency
    sumstats <- bind_rows(match_direct, match_flipped) %>%
      select(-A1_ref, -A2_ref) %>% distinct(SNP, .keep_all = TRUE)
    rm(ref_data, getP, match_direct, match_flipped)

    data <- get_ldscores(sumstats, model_name)
    data <- data %>% left_join(hm3, by = "SNP") %>% arrange(SNP) %>%
      distinct(SNP, .keep_all = TRUE)

    data$Z <- abs(qnorm(data$P / 2)) * sign(data$Beta)
    N_val <- median(data$N, na.rm = TRUE)
    outlier_idx <- which(data$Z^2 > pmax(0.001 * N_val, 80))
    if (length(outlier_idx) > 0) data$L2[outlier_idx] <- NA
    colnames(data) <- str_replace_all(colnames(data), "[-:]", "_")

    # layout below: [12 meta][annots][L2]
    meta_expect <- c("SNP", "CHR", "BP", "A1", "A2", "EAF", "N", "Beta", "se",
                     "P")
    if (!all(meta_expect %in% names(data)[1:12])) {
      stop("unexpected column layout; first 12 are: ",
           paste(names(data)[1:12], collapse = ", "))
    }
    if (names(data)[ncol(data)] != "L2") {
      stop("expected L2 as the final column, found ", names(data)[ncol(data)])
    }

    ldscores <- data[, -c(1:12)]
    stat <- data$Z^2
    z00 <- ldscores[, -ncol(ldscores)]
    weight_col <- ldscores[[ncol(ldscores)]]
    valid_rows <- !is.na(weight_col)

    score_col <- grepl("_PC_[0-9]+$", colnames(z00)) |
      colnames(z00) %in% c("ATAC", "CAGE", "CHIP_HISTONE", "CHIP_TF", "CONTACT_MAPS",
                           "DNASE", "PROCAP", "RNA_SEQ", "AVI",
                           "GPN_STAR_P", "GPN_STAR_M", "GPN_STAR_V")
    # A set that contains the baselineLD model keeps PolyFun's sum over every
    # column
    if (all(c("Coding_UCSC", "MAFbin1", "Conserved_LindbladToh") %in% colnames(z00))) {
      score_col[] <- FALSE
    }
    lds <- rowSums(z00[, !score_col, drop = FALSE])
    # On record in the log: what the weight's sum ran over, and what it left out
    message(sprintf("  weights: LD-score sum over %d of %d columns: %s", sum(!score_col), ncol(z00),
                    paste(colnames(z00)[!score_col], collapse = " ")))
    message(sprintf("  weights: %d model-score columns left out%s", sum(score_col),
                    if (any(score_col)) paste0(": ", paste(colnames(z00)[score_col], collapse = " ")) else ""))
    w_ld <- pmax(weight_col, 1)
    coef <- (mean(stat[valid_rows]) - 1) / (data$N[1] * mean(lds[valid_rows]))
    coef <- max(coef, 0)
    lds <- pmax(lds, 1)
    w_vec <- 1 / (2 * (1 + data$N[1] * coef * lds)^2 * w_ld)

    model_df <- z00
    model_df$chisq <- stat
    fmod <- as.formula(paste("chisq~", paste(colnames(z00), collapse = "+")))
    echisq_ridge <- rep(NA, nrow(model_df))
    CHR0 <- data$CHR
    coef_list <- list()

    # Ridge when regularization matters, else OLS.
    n_predictors <- ncol(z00)
    RIDGE_MIN_PREDICTORS <- as.integer(Sys.getenv("RIDGE_MIN_PREDICTORS", "5"))
    use_ridge <- n_predictors >= RIDGE_MIN_PREDICTORS
    message(sprintf("Found %d predictors (threshold %d). Using %s.",
                    n_predictors, RIDGE_MIN_PREDICTORS,
                    ifelse(use_ridge, "Ridge (cv.glmnet)",
                           "OLS (glmnet lambda=1e-10)")))

    for (chr in 1:22) {
      train_idx <- which(CHR0 != chr & valid_rows)
      test_idx <- which(CHR0 == chr)
      if (length(train_idx) == 0) {
        next
      }

      if (use_ridge) {
        rfit <- glmnetUtils::cv.glmnet(fmod, weights = w_vec[train_idx],
                                       alpha = 0,
                                       foldid = as.integer(as.factor(CHR0[train_idx])),
                                       data = model_df[train_idx, ])
        best_lambda_val <- rfit$lambda.min
        echisq_ridge[test_idx] <- predict(rfit, newdata = model_df[test_idx, ],
                                          s = "lambda.min", type = "response")
        raw_coefs <- as.matrix(coef(rfit, s = "lambda.min"))
      } else {
        rfit <- glmnetUtils::glmnet(fmod, data = model_df[train_idx, ],
                                    weights = w_vec[train_idx], alpha = 0,
                                    lambda = 1e-10)
        best_lambda_val <- 1e-10
        echisq_ridge[test_idx] <- predict(rfit, newdata = model_df[test_idx, ],
                                          s = 1e-10, type = "response")
        raw_coefs <- as.matrix(coef(rfit, s = 1e-10))
      }

      pred_names <- rownames(raw_coefs)[rownames(raw_coefs) != "(Intercept)"]
      train_sds <- apply(z00[train_idx, pred_names, drop = FALSE], 2, sd,
                         na.rm = TRUE)
      std_coefs_vec <- raw_coefs[pred_names, 1] * train_sds

      coef_list[[chr]] <- tibble(
        Annotation = pred_names,
        Beta_Raw = raw_coefs[pred_names, 1],
        Beta_Std = std_coefs_vec,
        Best_Lambda = best_lambda_val,
        LOCO_HeldOut_CHR = chr
      )
    }

    # a subset of results, to save space
    df_out <- data %>% select(SNP, CHR, BP, A1, A2, EAF, N, Beta, se, Z, P) %>%
      mutate(echisq = echisq_ridge)
    write_tsv(df_out, file.path(out_dir, paste0(sid, ".sffdr.gz")))

    if (length(coef_list) > 0) {
      df_coefs <- bind_rows(coef_list)
      df_coefs$Trait <- sid
      df_coefs$Model <- model_name
      write_tsv(df_coefs, file.path(out_dir, paste0(sid, ".coefs.gz")))
    }

    rm(df_out, data, model_df, coef_list)
    gc()

  }, error = function(e) {
    message("ERROR processing file: ", file, "\n", e$message)
  })
}

message("done: ", model_name)
