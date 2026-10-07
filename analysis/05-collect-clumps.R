##################################################
### File assembles the clumps into tables ###
##################################################
library(tidyverse)

# Cluster paths this analysis was run on. Set DATA_GWASAI (reference
# data) or ANDY (fit stores) to rebuild somewhere else.
BASE <- Sys.getenv("ANDY", "/home/ab3105/rds/rds-basis-YRWZsjDGyaU/People/ANDY")
# Tree-relative: this run's clump tree
OUT_ROOT <- file.path(normalizePath("."), "indep-discoveries-1000G")
AH_WINDOW <- 100000

collapse_ah <- function(chr, bp) {
  n <- length(chr)
  if (!n) {
    return(0L)
  }
  o <- order(chr, bp)
  chr <- chr[o]
  bp <- bp[o]
  if (n == 1L) {
    return(1L)
  }
  return(sum(c(TRUE, chr[-1] != chr[-n] | diff(bp) > AH_WINDOW)))
}

# Load the stage-1 clump files
failed <- list.files(OUT_ROOT, pattern = "\\.failed$", recursive = TRUE)
if (length(failed)) {
  stop(length(failed), " plink calls failed (rerun 04-clump.R); first: ",
       failed[1])
}
files <- list.files(OUT_ROOT, pattern = "\\.clumps$", recursive = TRUE,
                    full.names = TRUE)
files <- files[!str_detect(files, "/_paired/|/_union/")]
files <- files[!str_detect(basename(files), "\\.full-")]
message("clump files found: ", length(files))
if (!length(files)) {
  stop("nothing to collect under ", OUT_ROOT)
}

parse_one <- function(f) {
  rel <- str_remove(f, fixed(paste0(OUT_ROOT, "/")))
  # variant / model / <trait>-<chr>-<type>.clumps
  parts <- str_split(rel, "/")[[1]]
  stem <- str_remove(basename(f), fixed(".clumps"))
  bits <- str_split(stem, "-")[[1]]
  type <- bits[length(bits)]
  chr <- suppressWarnings(as.integer(bits[length(bits) - 1]))
  trait <- paste(bits[seq_len(length(bits) - 2)], collapse = "-")

  d <- suppressWarnings(suppressMessages(read_table(f, show_col_types = FALSE)))
  if (!nrow(d)) {
    return(NULL)
  }
  nm <- names(d)
  return(tibble(variant = parts[1], model = parts[2], trait = trait,
                type = type,
                CHR = d[[intersect(c("#CHROM", "CHR"), nm)[1]]],
                POS = d[[intersect(c("POS", "BP"), nm)[1]]],
                SNP = d[[intersect(c("ID", "SNP"), nm)[1]]],
                SP2 = if ("SP2" %in% nm) d$SP2 else NA_character_))
}

all_clumps <- map_dfr(files, parse_one) %>%
  mutate(CHR = as.integer(CHR), POS = as.numeric(POS))
write_tsv(all_clumps, file.path(OUT_ROOT, "all_clumps_1000G.tsv"))

# Combined region tables, downstream format
walk(unique(paste(all_clumps$variant, all_clumps$model, all_clumps$type)),
     function(k) {
  p <- str_split(k, " ")[[1]]
  d <- all_clumps %>% filter(variant == p[1], model == p[2], type == p[3]) %>%
    transmute(`disease type` = str_remove(trait, "^UKBB-"), `#CHROM` = CHR, POS,
              SNP, SP2)
  f <- file.path(OUT_ROOT,
                 sprintf("%s_%s_%s_combined_region_1000G.txt", p[1], p[2],
                         p[3]))
  write_tsv(d, f)
  message("wrote ", basename(f), "  (", nrow(d), " clumps)")
})

counts <- all_clumps %>% group_by(variant, model, trait, type) %>%
  summarise(n_lead = n(), n_loci = collapse_ah(CHR, POS), .groups = "drop")
# raw lead counts kept beside the collapsed counts, for sensitivity
write_tsv(counts, file.path(OUT_ROOT, "lead_counts_1000G.tsv"))

wide <- counts %>% select(-n_lead) %>%
  pivot_wider(names_from = type, values_from = n_loci, names_prefix = "n_",
              values_fill = 0) %>% # a type with no loci anywhere is 0
  mutate(pct_change = round(100 * (n_fP - n_P) / n_P, 2))
write_tsv(wide, file.path(OUT_ROOT, "locus_counts_1000G.tsv"))

empt <- list.files(OUT_ROOT, pattern = "\\.empty$", recursive = TRUE)
if (length(empt)) {
  cat(sprintf("\n%d chromosome-jobs produced no clumps (expected for chromosomes\n", length(empt)))
  cat("with no significant SNPs; investigate only if a large trait is affected).\n")
}
err <- list.files(OUT_ROOT, pattern = "\\.stderr$", recursive = TRUE,
                  full.names = TRUE)
err <- err[file.size(err) > 0]
if (length(err)) {
  cat(sprintf("\n!! %d plink calls left non-empty stderr; first few:\n",
              length(err)))
  print(head(err, 5))
}

cat("\n=== traits per variant/model (expect the full set) ===\n")
print(counts %>% distinct(variant, model, trait) %>% count(variant, model,
                                                           name = "n_traits"),
      n = 40)
# QC: P is the same in every model, so its locus count must agree across
# them
cat("\n=== QC: standard-P locus counts should be identical across models ===\n")
qc <- counts %>% filter(type == "P") %>%
  group_by(trait) %>%
  summarise(n_models = n(), n_distinct_counts = n_distinct(n_loci),
            counts = paste(sort(unique(n_loci)), collapse = ","),
            .groups = "drop")
print(qc %>% filter(n_distinct_counts > 1), n = 30)
if (all(qc$n_distinct_counts == 1)) {
  cat("PASS: all traits agree across models\n")
} else {
  cat("NOTE: traits above disagree. Small off-by-one differences usually mean the\n")
  cat("      input SNP sets differ slightly between models (annotation LD-score\n")
  cat("      merges do not all retain exactly the same SNPs). Check with the\n")
  cat("      input-count diagnostic before worrying.\n")
}

cat("\n=== % change in loci, functional vs standard p-values ===\n")
print(wide %>% group_by(variant, model) %>%
        summarise(n = n(),
                  median_pct = round(median(pct_change, na.rm = TRUE), 2),
                  q25 = round(quantile(pct_change, .25, na.rm = TRUE), 2),
                  q75 = round(quantile(pct_change, .75, na.rm = TRUE), 2),
                  .groups = "drop"), n = 40)
