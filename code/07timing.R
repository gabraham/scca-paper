
rm(list=ls())
graphics.off()

# Ensure the correct flashpcaR package is used
library(devtools)
install_url("file://flashpcaR_1.2.6.tar.gz", subdir=".")

library(flashpcaR)
library(PMA)
library(ggplot2)
library(microbenchmark)
library(plyr)

# v1.2.6 is first version to support initialisation of the V matrix in scca
stopifnot(packageVersion("flashpcaR") >= "1.2.6")

load("data.RData")

# faster crossprod using RcppEigen. Must be loaded after the data image.
source("crossprod.R")

################################################################################
# The actual timing experiments for subsets of the data

nsnps <- c(1000, 5000, 10000, 20000, 50000, 75000)
npheno <- c(100, 1000, 5000, 10000, 15000)

lambda1 <- 1e-2
lambda2 <- 1e-2

nreps <- 30
n1 <- 5000
n2 <- 50000

for(n in nsnps) {
   write.table(colnames(G)[1:n], file=paste0("snps_", n, ".txt"),
      quote=FALSE, col.names=FALSE, row.names=FALSE)
}

res <- lapply(nsnps, function(n) {
   Gs <- G[, 1:n]
   r <- lapply(npheno, function(p) {
      Es <- E[, 1:p]

      cat(n, p, "\n")

      K <- fprod(Es, Gs)
      f <- flashpca(K, ndim=1, verbose=FALSE, stand="none", check_geno=FALSE)
      Vinit <- f$vectors[,1]
      #Vinit <- rnorm(ncol(Es))
      # default time unit is milliseconds
      # For PMA::CCA, allow more iterations than the default 15 (scca
      # allows up to 1000 by default)
      bn <- microbenchmark(
         SCCA=scca(Gs, Es, ndim=1, lambda1=lambda1, lambda2=lambda2,
            verbose=FALSE, stand="none", check_geno=FALSE,
            V=Vinit),
         PMA_no_svd=CCA(Gs, Es, K=1, penaltyx=0.1, penaltyz=0.1,
            v=Vinit, niter=1000),
         times=nreps)
      print(bn)
      rm(Es)
      gc()
      data.frame(summary(bn, unit="s"), nsnps=n, npheno=p)
   })
   rm(Gs)
   gc()
   do.call(rbind, r)
})

vn <- c("expr", "median", "lq", "uq", "neval", "nsnps", "npheno")
res2 <- do.call(rbind, lapply(res, function(x) x[,vn]))

################################################################################
# Load the results from the commandline flashpca
d1 <- read.table(pipe("grep real scca*.log"), header=FALSE, sep="",
   stringsAsFactors=FALSE)
d1$label <- gsub("scca_chr|:real|\\.log|\\_rep[[:digit:]]+", "", d1$V1)
d1$time <- t(sapply(strsplit(d1$V2, "m|s"), as.numeric)) %*% c(60, 1)
d3 <- aggregate(time ~ label, FUN=median, data=d1)
colnames(d3)[2] <- "median"
d3$lq <- aggregate(time ~ label, FUN=quantile, probs=0.25, data=d1)[,2]
d3$uq <- aggregate(time ~ label, FUN=quantile, probs=0.75, data=d1)[,2]

d2 <- read.table(pipe("wc -l *chr*.bim"), header=FALSE, sep="",
   stringsAsFactors=FALSE)
d2 <- d2[grepl("\\.bim", d2$V2), ]
d2$label <- gsub(
   "hapmap3_r2_b36_fwd\\.consensus\\.qc\\.poly_founders_filtered_chr|\\.bim",
   "", d2$V2)
colnames(d2)[1:2] <- c("nsnps", "file")
colnames(d1)[1:2] <- c("logfile", "time")
d <- merge(d2, d3, by="label")

d$neval <- 30
d$expr <- "SCCA_cmdline"
d$npheno <- ncol(E)

res3 <- rbind(res2, d[, vn])

################################################################################
# Plotting
res4 <- subset(res3, npheno %in% c(1000, 10000, ncol(E)) & expr != "PMA")

# http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/#a-colorblind-friendly-palette
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73",
   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g <- ggplot(res4, aes(x=nsnps, y=median,
   #colour=factor(npheno), shape=expr, linetype=expr))
   colour=expr, shape=factor(npheno)))
g <- g + geom_line(size=1.5)
g <- g + geom_point(size=3)
#g <- g + geom_errorbar(aes(ymax=uq, ymin=lq))
g <- g + theme_bw()
g <- g + scale_x_log10("# of SNPs", breaks=10^c(2:6))
g <- g + scale_y_log10("Median wall time (seconds)", breaks=10^c(0:10))
#g <- g + scale_linetype_discrete(breaks=NULL)
#g <- g + scale_colour_discrete("Number of\nphenotypes")
g <- g + scale_colour_manual(values=cbbPalette, guide=FALSE)
g <- g + scale_shape_discrete("Number of\nphenotypes")
#g <- g + scale_shape_discrete(breaks=NULL)
g <- g + theme(legend.position=c(0.85, 0.15),
   plot.margin=unit(c(1, 5, 0, 0), "mm"))
# flashpcaR
#g <- g + annotate("segment", x=2e4, y=5e-2, xend=1e4, yend=0.4, colour="black",
#   size=0.5, arrow=arrow(length=unit(0.03, "npc")), alpha=0.3)
#g <- g + annotate("segment", x=2e4, y=5e-2, xend=5e3, yend=1.2, colour="black",
#   size=0.5, arrow=arrow(length=unit(0.03, "npc")), alpha=0.3)
g <- g + annotate("text", label="flashpcaR\n(R package)",
   x=3e4, y=0.3, colour="black")
#g <- g + annotate("point", x=1.3e4, y=0.36, colour=cbbPalette[1],
#   fill=cbbPalette[1], shape=19, size=2.5)
# flashpca
#g <- g + annotate("segment", x=3e5, y=2, xend=2e5, yend=15, colour="black",
#   size=0.5, arrow=arrow(length=unit(0.03, "npc")))
g <- g + annotate("text", label="flashpca\n(command line)",
   x=4.5e5, y=10, colour="black")
#g <- g + annotate("point", x=2e5, y=12, fill=cbbPalette[3],
#   colour=cbbPalette[3], shape=22, size=2.5)
# PMA
#g <- g + annotate("segment", x=1500, y=25, xend=3000, yend=4, colour="black",
#   size=0.5, arrow=arrow(length=unit(0.03, "npc")))
g <- g + annotate("text", label="PMA\n(R package)", x=9000, y=35)
#g <- g + annotate("point", x=5000, y=41, fill=cbbPalette[2],
#   colour=cbbPalette[2], shape=24, size=2.5)

pdf("scca_timing.pdf", width=5, height=5)
print(g)
dev.off()

save.image(file="timing.RData")

