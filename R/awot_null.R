#' Prepare null model for awot.
#' @param y - Continuous response variable.
#' @param tr - Binary treatment variable, 1 for treated, 0 for placebo.
#' @param x - Covariates in addition to treatment variable.
#' @return A list of objects needed for awot.
#' @references Hong Zhang, Devan Mehrotra and Judong Shen, "AWOT and CWOT for Genotype and Genotype by Treatment Interaction Joint Analysis in Pharmacogenetics GWAS".
#' @examples
#' n = 100
#' y = rnorm(n)
#' x = data.frame(x1=rnorm(n))
#' tr = rbinom(n, 1, 0.5)
#' nullmod = awot_null(y, tr, x)
#' @export
#' @import stats

awot_null = function(y, tr, x=NULL){
  n = length(y)
  if(is.null(x)){
    nullcovar = data.frame(tr=tr)
  }else{
    nullcovar = data.frame(tr=tr,x=x)
  }
  X1 = as.matrix(cbind(1, nullcovar))
  mod = glm(y~., data=nullcovar, family = "gaussian")
  s2 = sigma(mod)^2
  res = resid(mod)
  P0 = (diag(n) - X1%*%solve(t(X1)%*%X1)%*%t(X1))*s2
  return(list(res=as.matrix(res, ncol=1), P0=P0, s2=s2, tr=tr))
}
