
rm(list=ls())
graphics.off()

library(plink2R)
library(flashpcaR)
library(PMA)
library(fields)
library(doMC)

source("crossprod.R")

registerDoMC(cores=10)



################################################################################
# Load expression data
lf <- list.files(pattern="_p3_expression\\.txt", path="../data", full.name=TRUE)

res <- vector("list", length(lf))
names(res) <- lf
for(f in lf) {
   h <- read.table(f, nrows=1, stringsAsFactors=FALSE, sep="")
   x <- read.table(f, skip=2, stringsAsFactors=FALSE, sep="", row.names=1)
   colnames(x) <- h[-(1:2)]
   res[[f]] <- x
}

expr <- do.call(cbind, res)
expr <- t(expr)

################################################################################
# Load genotype data
dat <- read_plink(
   "../data/hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered_chr1",
   impute="random")

rownames(dat$bed) <- dat$fam[,2]
rownames(dat$fam) <- dat$fam[,2]
pop <- gsub("\\.\\./data/", "", sapply(strsplit(rownames(expr), "_"), head, n=1))
rownames(expr) <- sapply(strsplit(rownames(expr), "\\."), tail, n=1)
names(pop) <- rownames(expr)
pop <- factor(pop)

# Get the individuals common to the expression and genotype data (may not
# necessarily be in same order as PLINK may re-order samples)
nm <- rownames(expr)[rownames(expr) %in% dat$fam[,2]]
expr <- expr[nm, ]
pop <- pop[nm]
dat$fam <- dat$fam[nm, ]
dat$bed <- dat$bed[nm, ]

for(p in levels(pop)) {
   write.table(dat$fam[pop == p, 1:2], file=paste0("samples_", p, ".txt"),
      col.names=FALSE, row.names=FALSE, quote=FALSE)
}

# shorter names
E <- expr
G <- dat$bed

stopifnot(all(names(pop) == rownames(G)))
stopifnot(all(rownames(E) == rownames(G)))

################################################################################
# PCA of the gene expression data
pca.E <- flashpca(E, stand="none", ndim=10, nextra=50)

png("HM3_expression_PCA.png", width=1200, height=1200, res=200)
colnames(pca.E$projection) <- paste0("PC", 1:ncol(pca.E$projection))
pairs(pca.E$projection[, 1:5], gap=0, cex=1, pch=20, col=pop)
dev.off()

fam <- dat$fam

save(G, E, pop, pca.E, fam, nm, file="data_stage1.RData")

#
## now run flashpca commandline PCA
#
## Following Stranger2012, adjust certain populations for admixture using PCs
## of the genotypes 
#
#pn <- c("GIH", "LWK", "MEX", "MKK")
#E.pn <- foreach(n=levels(pop), .combine="rbind") %dopar% {
#
#   En <- E[pop == n, ]
#   if(n %in% pn) {
#      Gp <- as.matrix(read.table(paste0("pcs_", n, ".txt"), header=FALSE,
#         sep=""))
#      fam <- read.table(
#         paste0("hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered_maf0.05_",
#            n, "_thinned.fam"), header=FALSE, sep="", stringsAsFactors=FALSE)
#      rownames(Gp) <- fam[,2]
#      nm <- names(pop[pop == n])
#      Gp2 <- Gp[fam[,2] %in% nm, ]
#      Gp2 <- Gp2[nm, ]
#      En <- regress.out(cbind(1, Gp2), En)
#   }
#   scale(En)
#   #En
#}
#E.resid <- E.pn[rownames(E), ]
#
#stopifnot(all(rownames(E.resid) == names(pop)))
#
#pca.E.resid <- flashpca(E.resid, stand="center", ndim=10, nextra=50)
#png("HM3_expression_PCA_residuals.png", width=2400, height=1200, res=200)
#colnames(pca.E.resid$projection) <- paste0("PC", 1:ncol(pca.E.resid$projection))
#par(mfrow=c(2, 4))
##pairs(pca.E.resid$projection[, 1:5], gap=0, cex=1, pch=20, col=pop)
#for(i in 2:9) {
#   plot(pca.E.resid$proj[, 1], pca.E.resid$proj[, i], col=pop,
#      xlab="PC1", ylab=paste0("PC", i))
#   legend("topright", fill=seq(levels(pop)), legend=levels(pop))
#}
#dev.off()
#
#pc <- as.matrix(
#   read.table("pcs_hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered_maf0.05_thinned.txt",
#   header=FALSE, sep="", stringsAsFactors=FALSE))
#fam <- read.table("hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered_maf0.05_thinned.fam",
#   header=FALSE, sep="", stringsAsFactors=FALSE)
#rownames(pc) <- fam[,2]
#pc <- pc[rownames(E), ]
#
#M <- model.matrix(~ pop + pc[, 1:10])
#E.resid2 <- regress.out(M, E)
#
#pca.E.resid2 <- flashpca(E.resid2, stand="center", ndim=10, nextra=50)
#colnames(pca.E.resid2$projection) <- paste0("PC", 1:ncol(pca.E.resid2$projection))
#
#png("HM3_expression_PCA_residuals2.png", width=2400, height=1200, res=200)
#par(mfrow=c(2, 4))
##pairs(pca.E.resid2$projection[, 1:5], gap=0, cex=1, pch=20, col=pop)
#for(i in 2:9) {
#   plot(pca.E.resid2$proj[, 1], pca.E.resid2$proj[, i], col=pop,
#      xlab="PC1", ylab=paste0("PC", i))
#   legend("topright", fill=seq(levels(pop)), legend=levels(pop))
#}
#dev.off()
#
#
#################################################################################
## Set expression as the residuals
#E <- E.resid
#
#################################################################################
#
#ve <- apply(E, 2, sd) > 0.1
#table(ve)
#E <- E[, ve]
#
#dim(G)
#dim(E)
#
## We scale both G and E
#E <- scale(E)
#G <- scale(G)
#
#write.table(data.frame(dat$fam[nm, 1:2], format(E, digits=5)),
#   file="expression_standardised.txt", row.names=FALSE, col.names=FALSE,
#   quote=FALSE)
#write.table(dat$fam[nm, 1:2], file="common_samples.txt",
#   row.names=FALSE, col.names=FALSE, quote=FALSE)
#
## Regress out 10 PCs (of expression) from the expression data
#E.resid <- regress.out(cbind(1, pca.E$projection), E)
#
#save.image(file="data.RData")
#
