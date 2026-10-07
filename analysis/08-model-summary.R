##################################################
### File summarizes models across traits ###
##################################################
library(tidyverse)
library(data.table)

# Configuration
VARIANT <- Sys.getenv("VARIANT", "joint")
if (!VARIANT %in% c("joint", "subset")) {
  stop("VARIANT must be 'joint' or 'subset'")
}

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) or ANDY (fit stores) to rebuild somewhere else.
DATA_GWASAI <- Sys.getenv("DATA_GWASAI",
                          "/home/ab3105/rds/hpc-work/data_gwasAI")
SFFDR <- file.path(Sys.getenv("ANDY", "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY"), "sffdr")
base_dir2 <- file.path(SFFDR, "echisq") # baseline reference
base_dir <- if (VARIANT == "subset") {
  file.path(SFFDR, "echisq_subset")
} else {
  file.path(SFFDR, "echisq")
}
baseline_model <- "baselineLD2.2"

SUFFIX <- if (VARIANT == "subset") "_subset" else ""
OUT_R2 <- paste0("trait_level_r2_metrics_ridge", SUFFIX, ".tsv")
OUT_ANN <- paste0("summary_annotation_stats", SUFFIX, ".tsv")

# The arms of this variant, the paper's unless ARMS=all
ARMS_TSV <- Sys.getenv("ARMS_TSV", "annotations/arms.tsv")
ARM_SET <- Sys.getenv("ARMS", "paper")
arms <- read_tsv(ARMS_TSV, show_col_types = FALSE) %>%
  filter(variant == VARIANT, model != baseline_model)
if (ARM_SET == "paper") {
  arms <- arms %>% filter(paper == 1)
}
foundational_models <- arms$model
# Drop models with no directory here
have <- foundational_models[dir.exists(file.path(base_dir,
                                                 foundational_models))]
absent <- setdiff(foundational_models, have)
if (length(absent)) {
  message("!! no ", VARIANT, " directory (skipped): ",
          paste(absent, collapse = ", "))
}
foundational_models <- have
if (!length(foundational_models)) {
  stop("no model directories found under ", base_dir)
}

message("variant : ", VARIANT)
message("baseline: ", file.path(base_dir2, baseline_model))
message("models  : ", base_dir)
message("outputs : ", OUT_R2, ", ", OUT_ANN)

# Traits from the baseline coefs files. '.full.' fits have g = 1, so
# their p-values are NA and inflate the denominator
baseline_files <- list.files(file.path(base_dir2, baseline_model),
                             pattern = "\\.coefs\\.gz$")
n_all <- length(baseline_files)
baseline_files <- baseline_files[!str_detect(baseline_files, "\\.full\\.")]
traits <- str_replace(baseline_files, "\\.coefs\\.gz$", "")

message(sprintf("coefs files in %s: %d total, %d after dropping '.full.' variants",
                baseline_model, n_all, length(baseline_files)))
if (n_all != length(baseline_files)) {
  message(sprintf("  -> excluded %d '.full.' files (whole-data fits, 1 fold)",
                  n_all - length(baseline_files)))
}
if (any(str_detect(traits, "\\.full"))) {
  stop("trait list still contains '.full' entries -- check the filter")
}

message(paste("Found", length(traits), "traits to process."))
message(paste("Models to process:",
              paste(foundational_models, collapse = ", ")))

# Per-model coverage, so gaps show up front
for (m in c(baseline_model, foundational_models)) {
  d <- file.path(base_dir, m)
  cf <- list.files(d, pattern = "\\.coefs\\.gz$")
  cf <- cf[!str_detect(cf, "\\.full\\.")]
  sf <- list.files(d, pattern = "\\.sffdr\\.gz$")
  sf <- sf[!str_detect(sf, "\\.full\\.")]
  have <- str_replace(cf, "\\.coefs\\.gz$", "")
  miss <- setdiff(traits, have)
  message(sprintf("  %-32s coefs %2d, sffdr %2d%s", m, length(cf), length(sf),
                  if (length(miss)) paste0("   MISSING: ",
                                           paste(miss,
                                                 collapse = ", ")) else ""))
}

# HM3 weights for the weighted R2
message("Loading HM3 weights for weighted R2 calculation...")
weights_path <- file.path(DATA_GWASAI, "baselineLD2.3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.")
hm3_list <- lapply(1:22, function(i) {
  read_tsv(paste0(weights_path, i, ".l2.ldscore.gz"), col_types = cols()) %>%
    select(SNP, L2)
})
hm3_weights <- bind_rows(hm3_list)

# Raw ridge coefficients; folds kept for jackknifing
get_trait_coefficients <- function(model_path, trait, model_name) {
  # never accept a ".full" trait name
  if (str_detect(trait, "\\.full")) {
    return(NULL)
  }
  file_path <- file.path(model_path, paste0(trait, ".coefs.gz"))
  if (!file.exists(file_path)) {
    return(NULL)
  }

  df <- fread(file_path)
  out <- df %>%
    mutate(Trait = trait, Model = model_name) %>%
    select(Trait, Model, Annotation, Beta_Std, Beta_Raw, Best_Lambda)
  return(out)
}

# R2 and enrichment, baseline versus joint
get_r2_metrics <- function(baseline_dir, joint_dir, trait, joint_model_name,
                           hm3_weights) {
  if (str_detect(trait, "\\.full")) {
    return(NULL)
  }

  f_base <- file.path(baseline_dir, paste0(trait, ".sffdr.gz"))
  if (!file.exists(f_base)) {
    return(NULL)
  }

  df_base <- fread(f_base, select = c("SNP", "Z", "N", "echisq")) %>%
    filter(!is.na(Z), !is.infinite(Z)) %>%
    arrange(SNP) %>% distinct(SNP, .keep_all = TRUE)

  f_joint <- file.path(joint_dir, paste0(trait, ".sffdr.gz"))
  if (!file.exists(f_joint)) {
    return(NULL)
  }

  df_joint <- fread(f_joint, select = c("SNP", "echisq")) %>%
    rename(echisq_joint = echisq) %>%
    arrange(SNP) %>% distinct(SNP, .keep_all = TRUE)

  df <- inner_join(df_base, df_joint, by = "SNP") %>%
    inner_join(hm3_weights, by = "SNP") %>%
    filter(!is.na(L2))

  N_val <- median(df$N, na.rm = TRUE)
  max_chisq <- pmax(0.001 * N_val, 80)

  df <- df %>% mutate(obs_chisq = Z^2) %>% filter(obs_chisq <= max_chisq)

  obs_chisq <- df$obs_chisq
  w <- 1 / pmax(df$L2, 1)

  pred_base <- pmax(df$echisq, 1.0)
  pred_joint <- pmax(df$echisq_joint, 1.0)

  r2_pearson_base <- cor(pred_base, obs_chisq, use = "complete.obs")^2
  r2_pearson_joint <- cor(pred_joint, obs_chisq, use = "complete.obs")^2
  rho2_spearman_base <- cor(pred_base, obs_chisq, method = "spearman",
                            use = "complete.obs")^2
  rho2_spearman_joint <- cor(pred_joint, obs_chisq, method = "spearman",
                             use = "complete.obs")^2

  mean_obs_w <- sum(w * obs_chisq, na.rm = TRUE) / sum(w, na.rm = TRUE)
  sst_w <- sum(w * (obs_chisq - mean_obs_w)^2, na.rm = TRUE)
  sse_base_w <- sum(w * (obs_chisq - pred_base)^2, na.rm = TRUE)
  sse_joint_w <- sum(w * (obs_chisq - pred_joint)^2, na.rm = TRUE)
  strict_oos_R2_base <- 1 - (sse_base_w / sst_w)
  strict_oos_R2_joint <- 1 - (sse_joint_w / sst_w)

  wmae_base <- sum(w * abs(obs_chisq - pred_base), na.rm = TRUE) / sum(w, na.rm = TRUE)
  wmae_joint <- sum(w * abs(obs_chisq - pred_joint), na.rm = TRUE) / sum(w, na.rm = TRUE)

  dev_base <- 2 * sum(w * (obs_chisq * log((obs_chisq + 1e-9) / pred_base) - (obs_chisq - pred_base)), na.rm = TRUE) / sum(w, na.rm = TRUE)
  dev_joint <- 2 * sum(w * (obs_chisq * log((obs_chisq + 1e-9) / pred_joint) - (obs_chisq - pred_joint)), na.rm = TRUE) / sum(w, na.rm = TRUE)

  # Calibration over 100 bins of each model's own prediction: each model
  # is binned on its own ranking, so the two are compared on equal terms
  bin_stats <- function(pred) {
    b <- tibble(obs = obs_chisq, pred = pred, w = w) %>%
      mutate(Pred_Bin = ntile(pred, 100)) %>%
      group_by(Pred_Bin) %>%
      summarise(Mean_Obs = sum(w * obs) / sum(w),
                Mean_Pred = sum(w * pred) / sum(w), .groups = "drop")
    return(c(pearson = cor(b$Mean_Pred, b$Mean_Obs)^2,
             strict = 1 - sum((b$Mean_Obs - b$Mean_Pred)^2) /
               sum((b$Mean_Obs - mean(b$Mean_Obs))^2)))
  }
  bins_base <- bin_stats(pred_base)
  bins_joint <- bin_stats(pred_joint)
  bin_r2_pearson_base <- bins_base[["pearson"]]
  bin_r2_pearson_joint <- bins_joint[["pearson"]]
  bin_strict_R2_base <- bins_base[["strict"]]
  bin_strict_R2_joint <- bins_joint[["strict"]]

  df_ranked <- tibble(obs = obs_chisq, pred_b = pred_base, pred_j = pred_joint,
                      w = w) %>%
    mutate(Bin_Base = ntile(pred_b, 100), Bin_Joint = ntile(pred_j, 100))

  top10_obs_base <- sum(df_ranked$w[df_ranked$Bin_Base >= 91] * df_ranked$obs[df_ranked$Bin_Base >= 91]) / sum(df_ranked$w[df_ranked$Bin_Base >= 91])
  top10_obs_joint <- sum(df_ranked$w[df_ranked$Bin_Joint >= 91] * df_ranked$obs[df_ranked$Bin_Joint >= 91]) / sum(df_ranked$w[df_ranked$Bin_Joint >= 91])
  top5_obs_base <- sum(df_ranked$w[df_ranked$Bin_Base >= 96] * df_ranked$obs[df_ranked$Bin_Base >= 96]) / sum(df_ranked$w[df_ranked$Bin_Base >= 96])
  top5_obs_joint <- sum(df_ranked$w[df_ranked$Bin_Joint >= 96] * df_ranked$obs[df_ranked$Bin_Joint >= 96]) / sum(df_ranked$w[df_ranked$Bin_Joint >= 96])
  top1_obs_base <- sum(df_ranked$w[df_ranked$Bin_Base == 100] * df_ranked$obs[df_ranked$Bin_Base == 100]) / sum(df_ranked$w[df_ranked$Bin_Base == 100])
  top1_obs_joint <- sum(df_ranked$w[df_ranked$Bin_Joint == 100] * df_ranked$obs[df_ranked$Bin_Joint == 100]) / sum(df_ranked$w[df_ranked$Bin_Joint == 100])

  sum_w <- sum(w, na.rm = TRUE)
  wmse_base <- sse_base_w / sum_w
  wmse_joint <- sse_joint_w / sum_w

  out <- tibble(
    Trait = trait,
    Comparison = paste0("Baseline_vs_", joint_model_name),
    Pearson_r2_Base = r2_pearson_base, Pearson_r2_Joint = r2_pearson_joint,
    Delta_Pearson_r2 = r2_pearson_joint - r2_pearson_base,
    Strict_OOS_R2_Base = strict_oos_R2_base,
    Strict_OOS_R2_Joint = strict_oos_R2_joint,
    Delta_Strict_OOS_R2 = strict_oos_R2_joint - strict_oos_R2_base,
    Bin_Pearson_R2_Base = bin_r2_pearson_base,
    Bin_Pearson_R2_Joint = bin_r2_pearson_joint,
    Delta_Bin_Pearson_R2 = bin_r2_pearson_joint - bin_r2_pearson_base,
    Top10_Obs_Base = top10_obs_base, Top10_Obs_Joint = top10_obs_joint,
    Delta_Top10_Obs = top10_obs_joint - top10_obs_base,
    Top5_Obs_Base = top5_obs_base, Top5_Obs_Joint = top5_obs_joint,
    Delta_Top5_Obs = top5_obs_joint - top5_obs_base,
    Top1_Obs_Base = top1_obs_base, Top1_Obs_Joint = top1_obs_joint,
    Delta_Top1_Obs = top1_obs_joint - top1_obs_base,
    Bin_Strict_R2_Base = bin_strict_R2_base,
    Bin_Strict_R2_Joint = bin_strict_R2_joint,
    Delta_Bin_Strict_R2 = bin_strict_R2_joint - bin_strict_R2_base,
    WMSE_Base = wmse_base, WMSE_Joint = wmse_joint,
    Delta_WMSE = wmse_joint - wmse_base,
    WMAE_Base = wmae_base, WMAE_Joint = wmae_joint,
    Delta_WMAE = wmae_joint - wmae_base,
    Deviance_Base = dev_base, Deviance_Joint = dev_joint,
    Delta_Deviance = dev_joint - dev_base,
    Spearman_rho2_Base = rho2_spearman_base,
    Spearman_rho2_Joint = rho2_spearman_joint,
    Delta_Spearman_rho2 = rho2_spearman_joint - rho2_spearman_base
  )
  return(out)
}

# Main loop over traits and models
all_coefs_list <- list()
all_r2_list <- list()
skipped <- list()

message("Starting loop over traits...")

for (trait in traits) {
  base_coefs <- get_trait_coefficients(file.path(base_dir2, baseline_model),
                                       trait, baseline_model)
  if (!is.null(base_coefs)) all_coefs_list[[length(all_coefs_list) + 1]] <- base_coefs
  print(trait)

  for (f_model in foundational_models) {
    joint_coefs <- get_trait_coefficients(file.path(base_dir, f_model), trait,
                                          f_model)
    if (!is.null(joint_coefs)) {
      all_coefs_list[[length(all_coefs_list) + 1]] <- joint_coefs
    } else {
      skipped[[length(skipped) + 1]] <- tibble(Trait = trait, Model = f_model,
                                               what = "coefs")
    }

    r2_metrics <- get_r2_metrics(file.path(base_dir2, baseline_model),
                                 file.path(base_dir, f_model), trait, f_model,
                                 hm3_weights)
    if (!is.null(r2_metrics)) {
      all_r2_list[[length(all_r2_list) + 1]] <- r2_metrics
    } else {
      skipped[[length(skipped) + 1]] <- tibble(Trait = trait, Model = f_model,
                                               what = "r2")
    }
  }
}

trait_level_coefs <- bind_rows(all_coefs_list)
trait_level_r2 <- bind_rows(all_r2_list)

# Report anything skipped, never drop silently
if (length(skipped)) {
  sk <- bind_rows(skipped)
  message("\n!! SKIPPED (input file absent) -- these become missing traits downstream:")
  print(sk %>% count(Model, what), n = 40)
  print(sk, n = 60)
  write_tsv(sk, "aggregation_skipped.tsv")
} else message("\nno traits skipped -- all models have every input file.")

# Cross-trait jackknife summary statistics
message("Generating Jackknife summary statistics across all traits...")

summary_annotation_stats <- trait_level_coefs %>%
  group_by(Model, Trait, Annotation) %>%
  summarise(
    g = n(),
    Mean_Beta_Std = mean(Beta_Std, na.rm = TRUE),
    Median_Beta_Std = median(Beta_Std, na.rm = TRUE),
    SD_Beta_Std = sd(Beta_Std, na.rm = TRUE),
    Max_Beta_Std = max(abs(Beta_Std), na.rm = TRUE),
    Avg_Lambda = mean(Best_Lambda, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SE_Jackknife = SD_Beta_Std * ((g - 1) / sqrt(g)),
    Z_score = Mean_Beta_Std / SE_Jackknife,
    P_value = 2 * pt(-abs(Z_score), df = pmax(g - 1, 1))
  ) %>%
  group_by(Model, Trait) %>%
  mutate(FDR_Q_value = p.adjust(P_value, method = "fdr")) %>%
  ungroup() %>%
  group_by(Model, Trait) %>%
  mutate(Rank_Median = rank(-abs(Median_Beta_Std)),
         Rank_Max = rank(-abs(Max_Beta_Std))) %>%
  arrange(Model, Trait, Rank_Median)

# Every row must come from a 22-fold LOCO fit
bad_g <- summary_annotation_stats %>% filter(g != 22) %>% count(Model, Trait, g)
if (nrow(bad_g)) {
  message("\n!! rows where the number of folds is not 22 (a '.full' whole-data fit\n",
          "   has g = 1, which yields NA p-values and never reaches significance):")
  print(bad_g, n = 40)
} else message("\nvalidation: all rows have g = 22 folds. No '.full' contamination.")

full_rows <- summary_annotation_stats %>% filter(str_detect(Trait, "\\.full"))
if (nrow(full_rows)) {
  message("!! '.full' trait names still present -- the filter did not work:")
  print(full_rows %>% count(Model, Trait), n = 40)
}

cat("\ntrait counts per model (should be identical across models):\n")
print(summary_annotation_stats %>% distinct(Model, Trait) %>% count(Model),
      n = 40)
cat("\nR2 metric rows per model:\n")
print(trait_level_r2 %>% count(Comparison), n = 40)

# Save outputs
message("Saving files...")
write_tsv(trait_level_r2, OUT_R2)
write_tsv(summary_annotation_stats, OUT_ANN)
message("Done! Files created successfully.")
