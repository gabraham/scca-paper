

set.seed(1238)

rm(list=ls())
graphics.off()

library(plink2R)
library(flashpcaR)
library(PMA)
library(fields)
library(doMC)
library(abind)

# v1.2.6 is first version to support initialisation of the V matrix in scca
stopifnot(packageVersion("flashpcaR") >= "1.2.6")

load("data.RData")

source("crossprod.R")
source("makesim.R")

################################################################################
# Grid search for optimal penalties in cross-validation
nfolds <- 3
registerDoMC(cores=nfolds)
folds <- sample(nfolds, replace=TRUE, size=nrow(G))

lambda1 <- 0.02 * (1 - exp(seq(log(0.01), log(0.99), length=30)))
lambda2 <- 0.04 * (1 - exp(seq(log(0.01), log(0.99), length=25)))
v0 <- rnorm(ncol(E))

tm1 <- system.time({
   res1 <- foreach(fold=1:nfolds) %dopar% {
      trn <- folds != fold
      Gtrn <- scale(G[trn, ])
      Etrn <- scale(E[trn, ])
   
      Gtst <- scale(G[!trn, ],
         center=attr(Gtrn, "scaled:center"),
         scale=attr(Gtrn, "scaled:scale"))
      Etst <- scale(E[!trn, ],
         center=attr(Etrn, "scaled:center"),
         scale=attr(Etrn, "scaled:scale"))
   
      lapply(seq(along=lambda1), function(i) {
         lapply(seq(along=lambda2), function(j) {
	    s <- scca(Gtrn, Etrn, ndim=1,
	       lambda1=lambda1[i], lambda2=lambda2[j],
	       verbose=FALSE, stand="none", check_geno=FALSE, V=v0)
	    rownames(s$U) <- colnames(G)
            rownames(s$V) <- colnames(E)
   
            ptrn <- diag(cor(s$Px, s$Py))
	    Pxtst <- Gtst %*% s$U
	    Pytst <- Etst %*% s$V
            ptst <- diag(cor(Pxtst, Pytst))
   	    nzx <- colSums(s$U != 0)
   	    nzy <- colSums(s$V != 0)
            
   	    cat("fold:", fold, "i:", i, "j:", j,
   	       "nzx:", nzx, "nzy:", nzy, "ptst:", ptst, "\n")
   	    
	    list(ptrn=ptrn, ptst=ptst, nzx=nzx, nzy=nzy,
	       Pxtst=Pxtst, Pytst=Pytst, model=s,
	       cor.test=cor.test(Pxtst, Pytst))
         })
      })
   }
})
saveRDS(res1, file="res1.rds")

nzx1 <- abind(lapply(res1, sapply, sapply, function(x) x$nzx[1]), along=3)
nzy1 <- abind(lapply(res1, sapply, sapply, function(x) x$nzy[1]), along=3)
ptrn1 <- abind(lapply(res1, sapply, sapply, function(x) x$ptrn[1]), along=3)
ptst1 <- abind(lapply(res1, sapply, sapply, function(x) x$ptst[1]), along=3)
ptst1.avg <- apply(ptst1, 1:2, mean, na.rm=TRUE)
(ptst1.avg.max <- max(ptst1.avg, na.rm=TRUE))
w1 <- which(ptst1.avg == ptst1.avg.max, arr.ind=TRUE)

