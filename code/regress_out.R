
# Ridge regression, regressing out X from Y.
# Intercept is not added here, but should be part of X.
regress.out <- function(X, Y, lambda=1e-6)
{
   Y - X %*% solve(crossprod2(X, X) + lambda * diag(ncol(X)), crossprod2(X, Y))
}


