#!/usr/bin/env Rscript

# Calculate OrA overlap enrichment from mapped background and DEG flags.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
x <- fread(args[1])
x[, DEG := is.finite(log2FoldChange) & is.finite(pvalue) &
           abs(log2FoldChange) > 1 & pvalue < 0.05]

ans <- rbindlist(lapply(unique(x$condition), function(cc) {
  rbindlist(lapply(c("body", "flank2k"), function(mode) {
    d <- x[condition == cc & get(paste0(mode, "_mapped"))]
    hit <- d[[paste0(mode, "_hit")]]
    N <- nrow(d)
    K <- sum(hit)
    n <- sum(d$DEG)
    k <- sum(d$DEG & hit)
    data.table(
      condition = cc,
      overlap_mode = mode,
      background_N = N,
      background_overlap_K = K,
      DEG_n = n,
      DEG_overlap_k = k,
      fold_enrichment = (k / n) / (K / N),
      P = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    )
  }))
}))

ans[, FDR := p.adjust(P, method = "BH")]
fwrite(ans, args[2], sep = "\t")
