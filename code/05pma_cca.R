

set.seed(1238)

rm(list=ls())
graphics.off()

library(plink2R)
library(flashpcaR)
library(PMA)
library(fields)
library(doMC)
library(abind)



load("data.RData")

source("crossprod.R")
source("makesim.R")

################################################################################
# Grid search for optimal penalties in cross-validation
nfolds <- 3
registerDoMC(cores=nfolds)
folds <- sample(nfolds, replace=TRUE, size=nrow(G))

cons1 <- exp(seq(log(1e-3), log(1 - 1e-3), length=30))
cons2 <- exp(seq(log(1e-3), log(1 - 1e-3), length=25))
v0 <- rnorm(ncol(E))

tm2 <- system.time({
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
   
      lapply(seq(along=cons1), function(i) {
         lapply(seq(along=cons2), function(j) {
	    s <- CCA(Gtrn, Etrn, K=1, penaltyx=cons1[i], penaltyz=cons2[j],
	       standardize=FALSE, v=v0, niter=1e3)
	    rownames(s$u) <- colnames(G)
	    rownames(s$v) <- colnames(E)
   
	    ptrn <- diag(cor(Gtrn %*% s$u, Etrn %*% s$v))
	    Pxtst <- Gtst %*% s$u
	    Pytst <- Etst %*% s$v
            ptst <- diag(cor(Pxtst, Pytst))
	    nzx <- colSums(s$u != 0)
	    nzy <- colSums(s$v != 0)
         
	    cat("fold:", fold, "i:", i, "j:", j,
	       "nzx:", nzx, "nzy:", nzy, "ptst:", ptst, "\n")
   	 
            list(ptrn=ptrn, ptst=ptst, nzx=nzx, nzy=nzy,
	       Pxtst=Pxtst, Pytst=Pytst, model=s,
	       cor.test=cor.test(Pxtst, Pytst))
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
(ptst2.avg.max <- max(ptst2.avg, na.rm=TRUE))
w2 <- which(ptst2.avg == ptst2.avg.max, arr.ind=TRUE)

stop()

################################################################################
# Compare with CCA permutation strategy 
pen <- seq(1e-3, 0.5, length=25)
system.time({
   r.prm <- CCA.permute(G, E, penaltyx=pen, penaltyz=pen, nperms=100)
})
saveRDS(r.prm, file="PMA_CCA_permutations.rds")


