
rm(list=ls())
graphics.off()

library(plink2R)
library(flashpcaR)

################################################################################
# Load expression data
lf <- list.files(pattern="^(CEU|CHB|GIH|JPT|LWK|MEX|MKK|YRI)_p3_expression\\.txt",
   path="../data", full.name=TRUE)

res <- vector("list", length(lf))
names(res) <- lf
for(f in lf) {
   cat("reading", f, "\n")
   h <- read.table(f, nrows=1, stringsAsFactors=FALSE, sep="")
   x <- read.table(f, skip=2, stringsAsFactors=FALSE, sep="", row.names=1)
   colnames(x) <- h[-(1:2)]
   res[[f]] <- x
}

expr <- do.call(cbind, res)
expr <- t(expr)

################################################################################
# Load genotype data for chr 1
dat <- read_plink(
   "../data/hapmap3_r2_b36_fwd.qc.poly/hapmap3_r3_b36_fwd.qc.poly.filtered.wexpr.chr1",
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

# shorter names
E <- expr
G <- dat$bed
fam <- dat$fam

stopifnot(all(names(pop) == rownames(G)))
stopifnot(all(rownames(E) == rownames(G)))

################################################################################
# PCA of the gene expression data
pca.E <- flashpca(E, stand="none", ndim=10, nextra=50)

png("HM3_expression_PCA.png", width=1200, height=1200, res=200)
colnames(pca.E$projection) <- paste0("PC", 1:ncol(pca.E$projection))
pairs(pca.E$projection[, 1:5], gap=0, cex=1, pch=20, col=pop)
dev.off()


save(G, E, pop, pca.E, fam, nm, file="data_stage1.RData")

