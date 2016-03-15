
#
#library('Rcpp')
#library('inline')
#
#rcpp_inc <- '
#using namespace Rcpp;
#using namespace arma;
#'
#
#src <- '
#mat m1 = as<mat>(m1in);
#mat m2 = as<mat>(m2in);
#mat cp = trans(m1) * m2;
#return(wrap(cp));
#'
#fcrossprod <- cxxfunction(signature(m1in="numeric", m2in="numeric"), src,
#   plugin='RcppArmadillo', rcpp_inc)
#

library(Rcpp)
library(inline)

prodCpp <- '
typedef Eigen::Map<Eigen::MatrixXd>   MapMatd;
const MapMatd    B(as<MapMatd>(BB));
const MapMatd    C(as<MapMatd>(CC));
return wrap(B * C);
'

crossprodCpp <- '
typedef Eigen::Map<Eigen::MatrixXd>   MapMatd;
const MapMatd    B(as<MapMatd>(BB));
const MapMatd    C(as<MapMatd>(CC));
return wrap(B.transpose() * C);
'

prod2 <- cxxfunction(signature(BB="matrix", CC="matrix"),
   prodCpp, "RcppEigen")
crossprod2 <- cxxfunction(signature(BB="matrix", CC="matrix"),
   crossprodCpp, "RcppEigen")

