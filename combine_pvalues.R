library(rtracklayer)
library(GenomicRanges)

args <- commandArgs(TRUE)
out <- args[1]
input <- args[-1]

# constants
e <- exp(1)
log10e <- log10(e)

# sum replicates p-values
bdg_chi2 <- lapply(input, import, format = "bedGRaph")
bdg_chi2 <- GRangesList(bdg_chi2)
bdg_chi2 <- coverage(bdg_chi2, weight = "score")
bdg_chi2 <- GRanges(bdg_chi2)

# calculate Fisher's method p-value
bdg_chi2$score <- 2 * (bdg_chi2$score / log10(e))
bdg_chi2$score <- pchisq(bdg_chi2$score, df=2 * length(input), lower.tail=FALSE)
bdg_chi2$score[is.infinite(bdg_chi2$score)] <- 324 # -log10(1e-324) and below return Inf

# adjust p-value
idxes <- 1:length(bdg_chi2) + c(0, cumsum(width(bdg_chi2)[1:(length(bdg_chi2) - 1)]) - 1:(length(bdg_chi2) - 1))
pvalues <- rep(bdg_chi2$score, times = width(bdg_chi2))
qvalues <- p.adjust(pvalues, method = "BH")
qvalues <- qvalues[idxes]
qvalues <- -log10(qvalues)
qvalues[is.infinite(qvalues)] <- 324
bdg_chi2$score <- qvalues

export(bdg_chi2, con = out, format = "bedGRaph")
