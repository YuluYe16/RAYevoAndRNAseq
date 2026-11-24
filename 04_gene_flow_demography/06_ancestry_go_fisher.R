#!/usr/bin/env Rscript
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 2L, !file.exists(a[2]))
d <- read.delim(a[1], check.names = FALSE)
stopifnot(nrow(d) == 3146L, !anyDuplicated(d$GO_ID),
          all(d$OrA_BP_annotated_genes == 2213L), all(d$TRJ_BP_annotated_genes == 1835L))
tests <- lapply(seq_len(nrow(d)), function(i) {
    x <- d$OrA_count[i]; y <- d$TRJ_count[i]
    z <- rbind(c(x, 2213L-x), c(y, 1835L-y))
    f <- fisher.test(z, alternative = "two.sided", conf.level = 0.95)
    data.frame(p = f$p.value, OR = unname(f$estimate),
               OR_CI_low = f$conf.int[1], OR_CI_high = f$conf.int[2])
})
ans <- cbind(d, do.call(rbind, tests))
ans$BH_FDR_all3146 <- p.adjust(ans$p, method = "BH")
ans$OrA_vs_TRJ_difference_pp <- ans$OrA_percent - ans$TRJ_percent
write.table(ans, a[2], sep = "\t", row.names = FALSE, quote = FALSE)
