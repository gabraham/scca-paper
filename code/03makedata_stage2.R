

rm(list=ls())
graphics.off()

library(flashpcaR)
library(doMC)

registerDoMC(cores=10)

load("data_stage1.RData")

source("crossprod.R")
source("regress_out.R")

# Following Stranger2012, adjust certain populations for admixture using PCs
# of the genotypes 
pn <- c("GIH", "LWK", "MEX", "MKK")
E.pn <- foreach(n=levels(pop), .combine="rbind") %dopar% {

   En <- E[pop == n, ]
   if(n %in% pn) {
      Pc <- as.matrix(read.table(paste0("eigenvectors_", n, ".txt"), header=FALSE,
         sep=""))
      fam <- read.table(
         paste0("hapmap3_r2_b36_fwd.consensus.qc.poly_founders_filtered_",
            n, "_thinned.fam"), header=FALSE, sep="", stringsAsFactors=FALSE)
      rownames(Pc) <- fam[,2]
      nm <- names(pop[pop == n])
      Pc2 <- Pc[fam[,2] %in% nm, ]
      Pc2 <- Pc2[nm, ]
      En <- regress.out(cbind(1, Pc2), En)
   }
   En
}
E.resid <- E.pn[rownames(E), ]


stopifnot(all(rownames(E.resid) == names(pop)))

# Check for structure in the expression data
pca.E.resid <- flashpca(E.resid, stand="center", ndim=10, nextra=50)
png("HM3_expression_PCA_residuals.png", width=2400, height=1200, res=200)
colnames(pca.E.resid$projection) <- paste0("PC", 1:ncol(pca.E.resid$projection))
par(mfrow=c(2, 4))
for(i in 2:9) {
   plot(pca.E.resid$proj[, 1], pca.E.resid$proj[, i], col=pop,
      xlab="PC1", ylab=paste0("PC", i))
   legend("topright", fill=seq(levels(pop)), legend=levels(pop))
}
dev.off()

# Regress expression data on population indicator variable and 10 PCs of the
# expression data, take residuals
M <- model.matrix(~ pop + pca.E.resid$vectors[, 1:10])
E.resid.2 <- regress.out(M, E.resid)
pca.E.resid.2 <- flashpca(E.resid.2, stand="center", ndim=10, nextra=50)
png("HM3_expression_PCA_residuals2.png", width=2400, height=1200, res=200)
colnames(pca.E.resid.2$projection) <- paste0("PC", 1:ncol(pca.E.resid.2$projection))
par(mfrow=c(2, 4))
for(i in 2:9) {
   plot(pca.E.resid.2$proj[, 1], pca.E.resid$proj[, i], col=pop,
      xlab="PC1", ylab=paste0("PC", i))
   legend("topright", fill=seq(levels(pop)), legend=levels(pop))
}
dev.off()

################################################################################
# Set expression as the residuals
E <- E.resid.2

################################################################################

ve <- apply(E, 2, sd) > 0.1
table(ve)

E <- E[, ve]

dim(G)
dim(E)

# We scale both G and E
E <- scale(E)
G <- scale(G)

# Must write the phenotypes in the same order as the original fam file, as
# flashpca doesn't sort phenotype data
write.table(data.frame(fam[rownames(fam) %in% nm, 1:2], format(E, digits=5)),
   file="expression_standardised.txt", row.names=FALSE, col.names=FALSE,
   quote=FALSE)
write.table(fam[nm, 1:2], file="common_samples.txt",
   row.names=FALSE, col.names=FALSE, quote=FALSE)

save(G, E, pop, fam, file="data.RData")

