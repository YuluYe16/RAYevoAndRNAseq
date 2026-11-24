#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(DESeq2); library(clusterProfiler)})
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a)==4L, !file.exists(a[4]), dir.exists(dirname(a[4])))
x <- fread(a[1]); meta <- fread(a[2]); go <- fread(a[3])
stopifnot(all(c("sample","material","condition") %in% names(meta)),
          all(c("gene_id","GO_ID","GO_description") %in% names(go)),
          "gene_id" %in% names(x), !anyDuplicated(x$gene_id), !anyDuplicated(meta$sample))
if ("included" %in% names(meta)) meta <- meta[included == TRUE]
meta[, condition := factor(condition, levels=c("Control","Submergence"))]
samples <- meta$sample
stopifnot(nrow(meta)==12L, setequal(meta$material,c("RAY01","RAY02")),
          all(table(meta$material,meta$condition)==3L), !anyNA(meta$condition),
          setequal(setdiff(names(x),"gene_id"), samples))
m <- as.matrix(x[, ..samples]); rownames(m) <- x$gene_id
stopifnot(all(is.finite(m)), all(m >= 0), all(m == round(m)))
storage.mode(m) <- "integer"
stopifnot(dir.create(a[4]))
out <- function(x, f) fwrite(x, file.path(a[4], f), sep="\t", na="NA")
filter_counts <- function(y) y[rowSums(y >= 10L) >= 3L, , drop=FALSE]

cd <- as.data.frame(meta); rownames(cd) <- cd$sample
pc <- prcomp(t(assay(vst(DESeqDataSetFromMatrix(filter_counts(m),cd,~1),blind=TRUE))),scale.=FALSE)
out(as.data.table(pc$x[,1:3,drop=FALSE],keep.rownames="sample"),"PCA_coordinates.tsv")
out(data.table(PC=paste0("PC",seq_along(pc$sdev)),variance_fraction=pc$sdev^2/sum(pc$sdev^2)),
    "PCA_variance.tsv")

de <- lapply(c("RAY01","RAY02"), function(mat) {
    ids <- meta[material==mat, sample]
    dds <- DESeq(DESeqDataSetFromMatrix(filter_counts(m[,ids,drop=FALSE]),cd[ids,,drop=FALSE],~condition),quiet=TRUE)
    z <- as.data.table(as.data.frame(results(dds,contrast=c("condition","Submergence","Control"),alpha=.05)),
                       keep.rownames="gene_id")
    z[, status := "NS"]
    z[!is.na(padj) & padj < .05 & log2FoldChange >= 1, status := "UP"]
    z[!is.na(padj) & padj < .05 & log2FoldChange <= -1, status := "DOWN"]
    out(z,paste0(mat,"_DESeq2_all_genes.tsv"))
    out(z[status!="NS"],paste0(mat,"_DEG.tsv"))
    z
})
go <- unique(go[grepl("^GO:[0-9]{7}$",GO_ID) & !is.na(GO_description),
                .(gene_id,GO_ID,GO_description)])
bg <- Reduce(intersect,list(de[[1]][!is.na(padj),gene_id],
                           de[[2]][!is.na(padj),gene_id], go$gene_id))
for (direction in c("UP","DOWN")) {
    g <- intersect(de[[1]][status==direction,gene_id],de[[2]][status==direction,gene_id])
    fg <- intersect(g,bg)
    out(data.table(gene_id=g),paste0("Shared_",direction,"_input_genes.tsv"))
    out(data.table(gene_id=bg),paste0("Shared_",direction,"_GO_background.tsv"))
    stopifnot(length(fg)>0)
    er <- enricher(fg,universe=bg,TERM2GENE=go[,.(GO_ID,gene_id)],
                   TERM2NAME=unique(go[,.(GO_ID,GO_description)]),pvalueCutoff=1,qvalueCutoff=1,
                   pAdjustMethod="BH",minGSSize=3,maxGSSize=5000)
    out(as.data.table(as.data.frame(er)),paste0("Shared_",direction,"_GO_all_terms.tsv"))
}
writeLines(capture.output(sessionInfo()),file.path(a[4],"sessionInfo.txt"))
