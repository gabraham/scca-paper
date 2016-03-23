
set.seed(3448973)

n <- nrow(G)
ngene <- 1000
nsnps <- 10000

G <- G[, 1:nsnps]
h2 <- 0.1

# Can't make R exactly zero because that would lead to zero variance for M,
# and hence undefined heritability
R <- matrix(rnorm(ncol(G) * ngene), ncol=ngene) * sample(
   c(1e-4, 1), size=ngene * ncol(G), replace=TRUE, prob=c(0.9999, 0.0001))
M <- G %*% R
v <- apply(M, 2, var)
N <- matrix(rnorm(n * ngene), nrow=n)
V <- sqrt((1 - h2) / h2 * v)
E <- M + sweep(N, MARGIN=2, STATS=V, FUN="*")


