
load("data.RData")

source("crossprod.R")

trn <- sample(1:5, nrow(E), replace=TRUE)

Gtrn <- scale(G[trn != 1, ])
Etrn <- scale(E[trn != 1, ])
Gtst <- scale(G[trn == 1, ])
Etst <- scale(E[trn == 1, ])

l1 <- 1.4e-2
l2 <- 1.83e-2

n <- nrow(Gtrn)

v <- rnorm(ncol(E))

for(i in 1:25) {
   #u <- fprod(Gtrn, Etrn %*% v) / n
   u <- fprod(Gtrn, Etrn %*% v)
   if(any(abs(u) > 0)) {
      u <- u / sqrt(sum(u^2))
   }
   u <- sign(u) * pmax(abs(u) - l1, 0)
   if(all(u == 0)) {
      cat("u is all zero, stopping\n")
      break
   }
   if(any(abs(u) > 0)) {
      u <- u / sqrt(sum(u^2))
   }
   #v <- fprod(Etrn, Gtrn %*% u) / n
   v <- fprod(Etrn, Gtrn %*% u)
   if(any(abs(v) > 0)) {
      v <- v / sqrt(sum(v^2))
   }
   v <- sign(v) * pmax(abs(v) - l2, 0)
   if(any(abs(v) > 0)) {
      v <- v / sqrt(sum(v^2))
   }
   Xtrn <- Gtrn %*% u
   Ytrn <- Etrn %*% v
   ptrn <- cor(Xtrn, Ytrn)
   Xtst <- Gtst %*% u
   Ytst <- Etst %*% v
   ptst <- cor(Xtst, Ytst)
   cat(i, ptrn, ptst, sum(u != 0), sum(v != 0),"\n")
}

