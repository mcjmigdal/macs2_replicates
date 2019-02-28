#!/usr/bin/env R
# Migdal Feb 2019
# Combine p-values bedgraphs using Fisher method
# outputs combined bedgraph with p-values corrected using BH method
# Usage: Rscript combine_pvalues.R out.bedGraph in_1.bedGraph in_2.bedGraph ...
library(rtracklayer)
library(GenomicRanges)

args <- commandArgs(TRUE)
out <- args[1] # Output file name
input <- args[-1] # Input bedgraphs

stopifnot(length(input) > 2)

# constants
e <- exp(1)
log10e <- log10(e) # const for changing log10 to ln
K <- 2 * length(input) # degrees of freedom

# sum input bedGraphs
bdg_sum <- lapply(input, import, format = "bedGRaph")
bdg_sum <- GRangesList(bdg_sum)
bdg_sum <- coverage(bdg_sum, weight = "score")
bdg_sum <- GRanges(bdg_sum)

# calculate Fisher's method p-value
# chi2 = 2 * sum(ln(pval_1) + ln(pval_2) + ...)
bdg_sum$score <- 2 * (bdg_sum$score / log10e)
bdg_sum$score <- pchisq(bdg_sum$score, df = K, lower.tail=FALSE)
bdg_sum$score[is.infinite(bdg_sum$score)] <- 324 # -log10(1e-324) and below return Inf

# adjust p-value
idxes <- 1:length(bdg_sum) + c(0, cumsum(width(bdg_sum)[1:(length(bdg_sum) - 1)]) - 1:(length(bdg_sum) - 1))
pvalues <- rep(bdg_sum$score, times = width(bdg_sum))
qvalues <- p.adjust(pvalues, method = "BH")
qvalues <- qvalues[idxes]
qvalues <- -log10(qvalues)
qvalues[is.infinite(qvalues)] <- 324
bdg_sum$score <- qvalues

export(bdg_sum, con = out, format = "bedGRaph")
