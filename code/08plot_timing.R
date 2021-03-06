
rm(list=ls())

library(ggplot2)
library(plyr)

load("timing_tmp.RData")

################################################################################
# Ugly code to load the results from the commandline flashpca
d1 <- read.table(pipe("grep real scca_expression_standardised\\.txt*.log"), header=FALSE, sep="",
   stringsAsFactors=FALSE)
d1$label <- gsub(
   "scca_expression_standardised\\.txt_chr|:real|\\.log|\\_rep[[:digit:]]+",
   "", d1$V1)
d1$time <- t(sapply(strsplit(d1$V2, "m|s"), as.numeric)) %*% c(60, 1)
d3 <- aggregate(time ~ label, FUN=median, data=d1)
colnames(d3)[2] <- "median"
d3$lq <- aggregate(time ~ label, FUN=quantile, probs=0.25, data=d1)[,2]
d3$uq <- aggregate(time ~ label, FUN=quantile, probs=0.75, data=d1)[,2]

# Find out how many SNPs in each subset
d2 <- read.table(pipe("wc -l *chr*.bim"), header=FALSE, sep="",
   stringsAsFactors=FALSE)
d2 <- d2[grepl("\\.bim", d2$V2), ]
d2$label <- gsub(
   "hapmap3_r3_b36_fwd\\.qc\\.poly\\.filtered.wexpr_chr|\\.bim",
   "", d2$V2)
colnames(d2)[1:2] <- c("nsnps", "file")
colnames(d1)[1:2] <- c("logfile", "time")
d <- merge(d2, d3, by="label")
d$neval <- 30
d$expr <- "SCCAcmdline_rand"
d$npheno <- ncol(E)

# Read timing results for the experiments with the 10k genes
d1 <- read.table(
   pipe("grep real scca_expression_standardised_10000\\.txt*.log"),
   header=FALSE, sep="", stringsAsFactors=FALSE)
d1$label <- gsub(
   "scca_expression_standardised_10000\\.txt_chr|:real|\\.log|\\_rep[[:digit:]]+",
   "", d1$V1)
d1$time <- t(sapply(strsplit(d1$V2, "m|s"), as.numeric)) %*% c(60, 1)
d3 <- aggregate(time ~ label, FUN=median, data=d1)
colnames(d3)[2] <- "median"
d3$lq <- aggregate(time ~ label, FUN=quantile, probs=0.25, data=d1)[,2]
d3$uq <- aggregate(time ~ label, FUN=quantile, probs=0.75, data=d1)[,2]

# Find out how many SNPs in each subset, again... [TODO: duplicated]
d2 <- read.table(pipe("wc -l *chr*.bim"), header=FALSE, sep="",
   stringsAsFactors=FALSE)
d2 <- d2[grepl("\\.bim", d2$V2), ]
d2$label <- gsub(
   "hapmap3_r3_b36_fwd\\.qc\\.poly\\.filtered.wexpr_chr|\\.bim",
   "", d2$V2)
colnames(d2)[1:2] <- c("nsnps", "file")
colnames(d1)[1:2] <- c("logfile", "time")
d10k <- merge(d2, d3, by="label")

d10k$neval <- 30
d10k$expr <- "SCCAcmdline_rand"
d10k$npheno <- 10000

res3 <- rbind(res2, d[, vn], d10k[, vn])

res3$method <- sapply(strsplit(as.character(res3$expr), "_"), head, n=1)
res3$init <- sapply(strsplit(as.character(res3$expr), "_"), tail, n=1)
res3$expr <- factor(as.character(res3$expr))

################################################################################
# Plotting
res4 <- subset(res3, npheno %in% c(1000, 10000, ncol(E)) & expr != "PMA")

# http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/#a-colorblind-friendly-palette
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73",
   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

annot <- data.frame(
   init=rep(c("mean", "rand", "svd"), each=3),
   label=rep(c("flashpcaR\n(R package)", "flashpca\n(commandline)",
      "PMA\n(R package)"), 3),
   x=rep(c(2.3e4, 3.7e5, 2500), 3),
   y=rep(c(0.3, 6.1, 35), 3),
   colour="black",
   npheno=1000
)
annot$x[annot$label == "flashpca\n(commandline)" & annot$init != "rand"] <- NA

scale_y_log2 <- function(...)
{
   scale_y_continuous(..., trans=scales::log_trans(2))
}



g <- ggplot(res4, aes(x=nsnps, y=median, colour=method, shape=factor(npheno)))
g <- g + geom_line(size=1.5)
g <- g + geom_point(size=3)
g <- g + facet_wrap(~ init)
g <- g + theme_bw()
g <- g + scale_x_log10("# of SNPs", breaks=10^c(2:6))
g <- g + scale_y_log10("Median wall time (seconds)", breaks=10^c(-2:10))
#g <- g + scale_y_log2("Median wall time (seconds)", breaks=2^c(-4:10))
g <- g + scale_colour_manual(values=cbbPalette, guide=FALSE)
g <- g + scale_shape_discrete("Number of\nphenotypes")
g <- g + theme(legend.position=c(0.94, 0.17),
   plot.margin=unit(c(1, 5, 0, 0), "mm"))
#g <- g + annotate("text", label="flashpcaR\n(R package)",
#   x=2.3e4, y=0.3, colour="black")
#g <- g + annotate("text", label="flashpca\n(commandline)",
#   x=4.1e5, y=10, colour="black")
#g <- g + annotate("text", label="PMA\n(R package)", x=2500, y=35)
g <- g + geom_text(aes(x=x, y=y, label=label), colour="black", data=annot)

pdf("scca_timing_full.pdf", width=12, height=4.5)
print(g)
dev.off()

pdf("scca_timing_full_linearscale.pdf", width=10, height=4.5)
print(g + scale_x_continuous() + scale_y_continuous())
dev.off()

res5 <- subset(res3, npheno %in% c(1000, 10000, ncol(E)) & init == "rand")

g <- ggplot(res5, aes(x=nsnps, y=median, colour=expr, shape=factor(npheno)))
g <- g + geom_line(size=1.5)
g <- g + geom_point(size=3)
g <- g + theme_bw()
g <- g + scale_x_log10("# of SNPs", breaks=10^c(2:6))
g <- g + scale_y_log10("Median wall time (seconds)", breaks=10^c(-2:10))
#g <- g + scale_y_log2("Median wall time (seconds)", breaks=2^c(-4:10))
g <- g + scale_colour_manual(values=cbbPalette, guide=FALSE)
g <- g + scale_shape_discrete("Number of\nphenotypes")
g <- g + theme(legend.position=c(0.85, 0.15),
   plot.margin=unit(c(1, 5, 0, 0), "mm"))
g <- g + annotate("text", label="flashpcaR\n(R package)",
   x=3e4, y=0.8, colour="black")
g <- g + annotate("text", label="flashpca\n(commandline)",
   x=4.5e5, y=10, colour="black")
g <- g + annotate("text", label="PMA\n(R package)", x=3500, y=35)

pdf("scca_timing.pdf", width=5, height=5)
print(g)
dev.off()

save.image(file="timing.RData")

