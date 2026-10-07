##################################################
### File applies sfFDR to the GWAS p-values ###
##################################################
library(tidyverse)
library(sffdr)
library(digest)
options(warn = 1)

# Custom pi0: rank-transform covariates, anchor knots at anchor_value's
# rank, prune until each bin has enough SNPs
pi0_model_anchor <- function(
    z,
    indep_snps,
    anchor_value = 1.0,
    knots = c(0.50, 0.75, 0.90, 0.95, 0.99, 0.999),
    min_indep_per_bin = 100,
    verbose = TRUE
) {

  # rank transform
  z <- as.matrix(z)
  z_raw <- z

  z_ranked <- apply(z, 2, function(col) {
    n_valid <- sum(!is.na(col))
    if (n_valid == 0) {
      return(col)
    }
    (rank(col, ties.method = "random", na.last = "keep")) / n_valid
  })

  if (is.null(colnames(z))) colnames(z) <- paste0("z", seq_len(ncol(z)))
  colnames(z_ranked) <- make.names(colnames(z), unique = TRUE)
  selected_terms <- character(0)

  if (verbose) {
    message("==================================================")
    message(sprintf("Building Pi0 Model (Anchor Filtering)"))
    message(sprintf("  Anchor Value   : %.2f", anchor_value))
    message("--------------------------------------------------")
  }

  # iterate covariates
  for (i in seq_len(ncol(z))) {
    var_name <- colnames(z_ranked)[i]
    raw_vals <- z_raw[, i]
    rank_vals <- z_ranked[, i]

    valid_idx <- !is.na(rank_vals)
    valid_indices_map <- which(valid_idx)

    if (is.logical(indep_snps)) {
      curr_is_indep <- indep_snps[valid_indices_map]
    } else {
      curr_is_indep <- valid_indices_map %in% indep_snps
    }
    indep_ranks <- rank_vals[valid_idx][curr_is_indep]

    # anchor, then filter candidate knots
    anchor_rank <- mean(raw_vals[valid_idx] <= anchor_value, na.rm = TRUE)
    valid_candidates <- knots[knots > (anchor_rank + 0.01)]
    use_anchor <- (anchor_rank > 0.01 && anchor_rank < 0.99)

    if (use_anchor) {
      if (verbose) {
        message(sprintf("  [%s]: Anchor at Rank %.4f. Dropped candidates below anchor.",
                        var_name, anchor_rank))
      }
      current_knots <- sort(c(anchor_rank, valid_candidates))
    } else {
      current_knots <- sort(knots)
    }

    # prune until every bin has enough SNPs
    final_knots <- current_knots

    repeat {
      if (length(final_knots) == 0) {
        break
      }
      bins <- c(0, final_knots, 1)
      bin_cuts <- cut(indep_ranks, breaks = bins, include.lowest = TRUE)
      counts <- as.numeric(table(bin_cuts))

      if (all(counts >= min_indep_per_bin)) {
        break
      }

      bad_bin_idx <- which.min(counts)
      knot_to_remove <- NULL
      if (use_anchor) {
        anchor_idx <- which(final_knots == anchor_rank)
        if (length(anchor_idx) > 0 && bad_bin_idx == anchor_idx && anchor_idx > 1) {
          knot_to_remove <- anchor_idx - 1
        } else if (length(anchor_idx) > 0 && bad_bin_idx == (anchor_idx + 1) &&
                   anchor_idx < length(final_knots)) {
          knot_to_remove <- anchor_idx + 1
        }
      }

      if (is.null(knot_to_remove)) {
        if (bad_bin_idx > length(final_knots)) {
          knot_to_remove <- length(final_knots)
        } else {
          knot_to_remove <- bad_bin_idx
        }
      }

      if (final_knots[knot_to_remove] == anchor_rank) use_anchor <- FALSE
      final_knots <- final_knots[-knot_to_remove]
    }

    # spline, or linear if no knots survive
    if (length(final_knots) > 0) {
      knot_str <- paste(round(final_knots, 5), collapse = ", ")
      term <- sprintf("splines::ns(%s, knots = c(%s), Boundary.knots = c(0, 1))",
                      var_name, knot_str)
      if (verbose) message(sprintf("  [%s]: Final Model -> %s", var_name,
                                   knot_str))
    } else {
      term <- var_name
      if (verbose) message(sprintf("  [%s]: Final Model -> Linear", var_name))
    }
    selected_terms <- c(selected_terms, term)
  }

  if (length(selected_terms) == 0) {
    fmod <- formula("~ 1")
  } else {
    fmod <- formula(paste("~", paste(selected_terms, collapse = " + ")))
  }

  return(list(fmod = fmod, zt = as.data.frame(z_ranked)))
}

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) or ANDY (fit stores) to rebuild somewhere else.
ANDY <- Sys.getenv("ANDY", "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")

# Reads <trait>.sffdr.gz from IN_BASE, writes <trait>.sffdr_final.gz to
# OUT_BASE
VARIANT <- Sys.getenv("VARIANT", "subset")
if (!VARIANT %in% c("subset", "joint")) {
  stop("VARIANT must be 'subset' or 'joint'")
}
# The arms of this variant, from the arm table
ARMS_TSV <- Sys.getenv("ARMS_TSV", "annotations/arms.tsv")
MODELS <- read_tsv(ARMS_TSV, show_col_types = FALSE) %>%
  filter(variant == VARIANT) %>% pull(model)
IN_BASE <- if (VARIANT == "subset") {
  file.path(ANDY, "sffdr/echisq_subset")
} else {
  file.path(ANDY, "sffdr/echisq")
}
OUT_BASE <- if (VARIANT == "subset") {
  file.path(ANDY, "sffdr/final_results/subset")
} else {
  file.path(ANDY, "sffdr/final_results")
}

PRUNED_SNPS <- "./qc_panukk/pruned_snps_0.3.tsv"
TOTAL_JOBS <- 13

# Arguments: array id, optionally one model
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  id2 <- 1
  message("No args provided, defaulting to ID 1")
} else {
  id2 <- as.numeric(args[1])
}
if (length(args) >= 2) {
  if (!args[2] %in% MODELS) {
    stop("Unknown model '", args[2], "'. Allowed:\n  ",
         paste(MODELS, collapse = "\n  "))
  }
  MODELS <- args[2]
}

message("variant : ", VARIANT)
message("reading : ", IN_BASE)
message("writing : ", OUT_BASE)
message("models  : ", paste(MODELS, collapse = ", "))

# Gather files; report models with no input
present <- MODELS[dir.exists(file.path(IN_BASE, MODELS))]
absent <- setdiff(MODELS, present)
if (length(absent)) {
  message("!! no input directory (skipped): ", paste(absent, collapse = ", "))
}
if (!length(present)) {
  stop("no input directories found under ", IN_BASE)
}

all_files <- unlist(lapply(present, function(m) {
  f <- list.files(file.path(IN_BASE, m), pattern = "\\.sffdr\\.gz$",
                  full.names = TRUE)
  f <- f[!str_detect(basename(f), "\\.full\\.")] # skip whole-data fits
  message(sprintf("  %-28s %d traits", m, length(f)))
  f
}))
if (!length(all_files)) {
  stop("no .sffdr.gz files found under ", IN_BASE)
}

my_indices <- which((seq_along(all_files) %% TOTAL_JOBS) == id2)
my_files <- all_files[my_indices]
message(sprintf("This job (ID %s) will process %d of %d files.",
                id2, length(my_files), length(all_files)))

# Independent SNPs used to fit pi0
if (!file.exists(PRUNED_SNPS)) {
  stop("pruned SNP list not found at '", PRUNED_SNPS, "'\n",
       "  (path is relative to the working directory: ", getwd(), ")")
}

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

indep_snps_ref <- read_tsv(PRUNED_SNPS, show_col_types = FALSE)

# apply sffdr
n_done <- 0L
n_skip <- 0L
n_err <- 0L
for (file in my_files) {

  # skip before tryCatch: next() cannot cross it
  current_model <- basename(dirname(file))
  trait <- str_remove(basename(file), "\\.sffdr\\.gz$")
  output_dir <- file.path(OUT_BASE, current_model)
  out_file <- file.path(output_dir, paste0(trait, ".sffdr_final.gz"))
  if (file.exists(out_file)) {
    message("skip (exists): ", current_model, " / ", trait)
    n_skip <- n_skip + 1L
    next
  }

  tryCatch({
    message(paste("Processing:", trait, "Model:", current_model))
    # the rank transform breaks ties at random; seed per model and trait so
    # a rerun reproduces the fit
    set.seed(readBin(digest(paste(current_model, trait), raw = TRUE),
                     "integer"))

    df <- read_tsv(file, show_col_types = FALSE, col_types = p_cols(file, "P")) %>%
      mutate(P = exact_p(P))

    echisq_val <- df$echisq
    p_val <- df$P
    is_indep <- df$SNP %in% indep_snps_ref$SNP

    fpi0m_tail <- pi0_model_anchor(
      z = as.matrix(echisq_val),
      indep_snps = is_indep,
      anchor_value = 1.0,
      knots = c(0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999),
      min_indep_per_bin = 100
    )

    fpi0_ridge <- fpi0est(
      p = p_val,
      z = fpi0m_tail$zt,
      pi0_model = fpi0m_tail$fmod,
      indep_snps = is_indep
    )$fpi0

    out <- sffdr(p_val, surrogate = echisq_val, fpi0 = fpi0_ridge)

    df$fpi0 <- out$fpi0
    df$fpvalues <- out$fpvalues
    df$flfdr <- out$flfdr
    df$fqvalues <- out$fqvalues

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write_tsv(df, out_file)
    n_done <- n_done + 1L

    rm(df, out, fpi0m_tail, fpi0_ridge)
    gc()

  }, error = function(e) {
    n_err <<- n_err + 1L
    message("Error processing file: ", file)
    message("  message: ", conditionMessage(e))
    cl <- conditionCall(e)
    if (!is.null(cl)) message("  call   : ", paste(deparse(cl), collapse = " "))
  })
}

message(sprintf("\nfinished: %d written, %d skipped, %d errors", n_done, n_skip,
                n_err))
