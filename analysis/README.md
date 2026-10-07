# Code to reproduce the GWAS analysis

Overview of files:

- `01-ldsc-h2.R`: runs stratified LD score regression per trait and model 
- `02-ridge-echisq.R`: fits the expected chi-squared statistic
- `03-sffdr.R`: applies sfFDR
- `04-clump.R`: clumps significant SNPs into independent loci against the 1000G EUR panel
- `05-collect-clumps.R`: merges the per-chromosome clumps and collapses signals within 100 kb
- `06-lead-stats.R`: adds Z, N, alleles and both p-values to each clump lead
- `07-collapse-classify.R`: classifies each locus of the joint baselineLD fit as standard only, functional only, or shared
- `08-model-summary.R`: computes out-of-sample R², signal enrichment, and annotation coefficient statistics
- `09-annotation-overlap.R`: measures how much of a model's LD scores the baselineLD classes explain
