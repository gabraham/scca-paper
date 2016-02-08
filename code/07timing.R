
rm(list=ls())
graphics.off()

library(flashpcaR)
library(PMA)
library(ggplot2)
library(microbenchmark)
library(plyr)


load("data.RData")

# faster crossprod using RcppEigen. Must be loaded after the data image.
source("crossprod.R")

################################################################################
# The actual timing experiments for subsets of the data

nsnps <- c(1000, 5000, 10000, 20000, 50000)
npheno <- c(100, 5000, 10000)

lambda1 <- 1e-3
lambda2 <- 1e-3

nreps <- 30
n1 <- 5000
n2 <- 20000

for(n in nsnps) {
   write.table(colnames(G)[1:n], file=paste0("snps_", n, ".txt"),
      quote=FALSE, col.names=FALSE, row.names=FALSE)
}

res <- lapply(nsnps, function(n) {
   Gs <- G[, 1:n]
   r <- lapply(npheno, function(p) {
      Es <- E[, 1:p]

      cat(n, p, "\n")

      # Only run the PMA without-pre-SVD on small data, otherwise very slow
      if(n <= n1) {
	 K <- fprod(Es, Gs)
	 f <- flashpca(K, ndim=1, verbose=FALSE, stand="none", check_geno=FALSE)
      	 # default time unit is milliseconds
      	 bn <- microbenchmark(
      	    SCCA=scca(Gs, Es, ndim=1, lambda1=lambda1, lambda2=lambda2,
      	 	    verbose=FALSE, stand="none", check_geno=FALSE),
      	    PMA=CCA(Gs, Es, K=1, penaltyx=0.9, penaltyz=0.9),
      	    PMA_no_svd=CCA(Gs, Es, K=1, penaltyx=0.9, penaltyz=0.9,
	       v=f$vectors[,1]),
      	    times=nreps)
      } else if(n <= n2) {
	 K <- fprod(Es, Gs)
	 f <- flashpca(K, ndim=1, verbose=FALSE, stand="none", check_geno=FALSE)
      	 # default time unit is milliseconds
      	 bn <- microbenchmark(
      	    SCCA=scca(Gs, Es, ndim=1, lambda1=lambda1, lambda2=lambda2,
      	 	    verbose=FALSE, stand="none", check_geno=FALSE),
      	    PMA_no_svd=CCA(Gs, Es, K=1, penaltyx=0.9, penaltyz=0.9,
	       v=f$vectors[,1]),
      	    times=nreps)
      } else {
      	 bn <- microbenchmark(
      	    SCCA=scca(Gs, Es, ndim=1, lambda1=lambda1, lambda2=lambda2,
      	 	    verbose=FALSE, stand="none", check_geno=FALSE),
      	    times=nreps
	 )
      }
      print(bn)
      rm(Es)
      gc()
      if(n <= n2) {
	 data.frame(summary(bn, unit="s"), nsnps=n, npheno=p)
      } else {
	 data.frame(summary(bn, unit="s"), nsnps=n, npheno=p)
      }
   })
   rm(Gs)
   gc()
   do.call(rbind, r)
})

vn <- c("expr", "median", "neval", "nsnps", "npheno")
res2 <- do.call(rbind, lapply(res, function(x) x[,vn]))

################################################################################
# Load the results from the commandline flashpca
d1 <- read.table(pipe("grep real scca*.log"), header=FALSE, sep="",
   stringsAsFactors=FALSE)
d1$label <- gsub("scca_chr|:real|\\.log|\\_rep[[:digit:]]+", "", d1$V1)
d1$time <- t(sapply(strsplit(d1$V2, "m|s"), as.numeric)) %*% c(60, 1)
d3 <- aggregate(time ~ label, FUN=median, data=d1)
colnames(d3)[2] <- "median"

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
res4 <- subset(res3, npheno %in% c(100, 10000, ncol(E)) & expr != "PMA")

g <- ggplot(res4, aes(x=nsnps, y=median,
   colour=factor(npheno), shape=expr, linetype=expr))
g <- g + geom_line(size=1)
g <- g + geom_point(size=3)
g <- g + theme_bw()
g <- g + scale_x_log10("# of SNPs", breaks=10^c(2:6))
g <- g + scale_y_log10("Median wall time (seconds)",
   breaks=10^c(0:10))
g <- g + scale_linetype_discrete(breaks=NULL)
g <- g + scale_colour_discrete("Number of\nphenotypes")
g <- g + scale_shape_discrete(breaks=NULL)
g <- g + theme(legend.position=c(0.8, 0.2),
   plot.margin=unit(c(1, 5, 0, 0), "mm"))
g <- g + annotate("text", label="flashpcaR\n(R package)", x=2e4, y=3e-2, colour="black")
g <- g + annotate("text", label="flashpca\n(command line)", x=4.5e5, y=1.2,
   colour="black")
g <- g + annotate("text", label="PMA\n(R package)", x=4000, y=62)
g <- g + annotate("segment", x=2e4, y=5e-2, xend=1.5e4, yend=0.19, colour="black",
   size=0.5, arrow=arrow(length=unit(0.03, "npc")))
g <- g + annotate("segment", x=5e5, y=2, xend=5e5, yend=15, colour="black",
   size=0.5, arrow=arrow(length=unit(0.03, "npc")))
g <- g + annotate("segment", x=5000, y=30, xend=7000, yend=9, colour="black",
   size=0.5, arrow=arrow(length=unit(0.03, "npc")))
g <- g + annotate("point", x=8.4e3, y=3.6e-2, colour="black", shape=19,
   size=2.5, fill="black")
g <- g + annotate("point", x=1.7e5, y=1.4, colour="black", shape=24, size=2.5,
   fill="black")
g <- g + annotate("point", x=2000, y=75, colour="black", shape=22, size=2.5,
   fill="black")

pdf("scca_timing.pdf", width=5, height=5)
print(g)
dev.off()

save.image(file="timing.RData")

