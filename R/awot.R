#' Adaptively weighted joint test for main effect and genotype-by-treatment interaction effect for continuous endpoints.
#' @param nullmod - The null model object from the output of awot_null.
#' @param g - The variable of interest, e.g., the genotype.
#' @param weights - The pre-specified weights. The default choice is a vector of -1, -0.9,..., 0.9, 1.
#' @return The p-value of AWOT and the individual p-values of the composite genotypes.
#' @references Hong Zhang, Qing Li, Devan Mehrotra and Judong Shen. "CauchyCP: a powerful test under non-proportional hazards using Cauchy combination of change-point Cox regressions", arXiv:2101.00059.
#' @examples
#' n = 100
#' y = rbinom(n, 1, 0.3)
#' x = data.frame(x1=rnorm(n))
#' tr = rbinom(n, 1, 0.5)
#' g = rbinom(n, 2, 0.1)
#' nullmod = awot_null(y, tr, x)
#' awot(nullmod, g, weights=seq(-1,1,0.1))
#' @export
#' @import stats mvtnorm

awot = function(nullmod, g, weights=seq(-1,1,0.1)){
  y = nullmod$y
  tr = nullmod$tr
  gt = g*tr
  u = nullmod$res; P0 = nullmod$P0
  comp_covar = sapply(weights, function(s)s*g+(1-abs(s))*gt)
  mi = min(comp_covar)
  ma = max(comp_covar)
  comp_covar = (comp_covar - mi)*2/(ma - mi)
  COV = getCOV(tr, g, P0, weights)
  STAT = c(t(u)%*%comp_covar)/sqrt(diag(COV))
  PVAL = 2*pnorm(abs(STAT), lower.tail = FALSE)
  stat = min(PVAL)
  COR = cov2cor(COV)
  bound = rep(qnorm(1-stat/2), length(STAT))
  p_awot = 1 - c(pmvnorm(lower=-bound, upper=bound, mean=rep(0, length(bound)), corr=COR))
  return(list(p_awot=p_awot, PVAL_score_grids = PVAL))
}

