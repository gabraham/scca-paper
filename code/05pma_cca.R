

set.seed(1238)

rm(list=ls())
graphics.off()

library(plink2R)
library(flashpcaR)
library(PMA)
library(fields)
library(doMC)
library(abind)


registerDoMC(cores=10)

load("data.RData")

source("crossprod.R")

################################################################################
# Grid search for optimal penalties in cross-validation
nfolds <- 5
folds <- sample(nfolds, replace=TRUE, size=nrow(G))

# PMA
cons1 <- exp(seq(log(1e-3), log(1 - 1e-3), length=30))
cons2 <- exp(seq(log(1e-3), log(1 - 1e-3), length=25))

system.time({
   res2 <- foreach(fold=1:nfolds) %dopar% {
      trn <- folds != fold
      Gtrn <- scale(G[trn, ])
      Etrn <- scale(E[trn, ])
   
      Gtst <- scale(G[!trn, ],
         center=attr(Gtrn, "scaled:center"),
         scale=attr(Gtrn, "scaled:scale"))
      Etst <- scale(E[!trn, ],
         center=attr(Etrn, "scaled:center"),
         scale=attr(Etrn, "scaled:scale"))
   
      # SVD of X^T T, used to initialise CCA
      K <- fprod(Etrn, Gtrn)
      f <- flashpca(K, ndim=1, verbose=FALSE, stand="none", check_geno=FALSE)
   
      lapply(seq(along=cons1), function(i) {
         lapply(seq(along=cons2), function(j) {
	    s <- CCA(Gtrn, Etrn, K=1, penaltyx=cons1[i], penaltyz=cons2[j],
         	    standardize=FALSE, v=f$vectors[,1])
	    rownames(s$u) <- colnames(G)
	    rownames(s$v) <- colnames(E)
   
	    ptrn <- diag(cor(Gtrn %*% s$u, Etrn %*% s$v))
	    ptst <- diag(cor(Gtst %*% s$u, Etst %*% s$v))
	    nzx <- colSums(s$u != 0)
	    nzy <- colSums(s$v != 0)
         
	    cat("fold:", fold, "i:", i, "j:", j,
	       "nzx:", nzx, "nzy:", nzy, "ptst:", ptst, "\n")
   	 
            list(ptrn=ptrn, ptst=ptst, nzx=nzx, nzy=nzy)
         })
      })
   }
})
saveRDS(res2, file="res2.rds")

nzx2 <- abind(lapply(res2, sapply, sapply, function(x) x$nzx[1]), along=3)
nzy2 <- abind(lapply(res2, sapply, sapply, function(x) x$nzy[1]), along=3)
ptrn2 <- abind(lapply(res2, sapply, sapply, function(x) x$ptrn[1]), along=3)
ptst2 <- abind(lapply(res2, sapply, sapply, function(x) x$ptst[1]), along=3)
ptst2.avg <- apply(ptst2, 1:2, mean, na.rm=TRUE)
ptst2.avg.max <- max(ptst2.avg, na.rm=TRUE)
w2 <- which(ptst2.avg == ptst2.avg.max, arr.ind=TRUE)

################################################################################
# Compare with CCA permutation strategy 
pen <- seq(1e-3, 0.5, length=25)
r.prm <- CCA.permute(G, E, penaltyx=pen, penaltyz=pen, nperms=100)
saveRDS(r.prm, file="PMA_CCA_permutations.rds")


