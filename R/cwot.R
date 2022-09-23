#' Cauchy weighted joint test for main effect and genotype-by-treatment interaction effect for binary endpoints.
#' @param nullmod - The null model object from the output of cwot_null.
#' @param g - The variable of interest, e.g., the genotype.
#' @param weights - The pre-specified weights. The default choice is a vector of -1, -0.9,..., 0.9, 1.
#' @return The p-values of CWOT, CWOT_Score, CWOT_LRT and the individual p-values of the composite genotypes.
#' @references Hong Zhang, Qing Li, Devan Mehrotra and Judong Shen. "CauchyCP: a powerful test under non-proportional hazards using Cauchy combination of change-point Cox regressions", arXiv:2101.00059.
#' @examples
#' n = 100
#' y = rbinom(n, 1, 0.3)
#' x = data.frame(x1=rnorm(n))
#' tr = rbinom(n, 1, 0.5)
#' g = rbinom(n, 2, 0.1)
#' nullmod = cwot_null(y, tr, x)
#' cwot(nullmod, g, weights=seq(-1,1,0.1))
#' @export
#' @import stats SPAtest

cwot = function(nullmod, g, weights=seq(-1,1,0.1)){
  y = nullmod$y
  tr = nullmod$nullcovar$tr
  nullcovar = as.matrix(nullmod$nullcovar)
  gt = g*tr
  comp_covar = sapply(weights, function(s)s*g+(1-abs(s))*gt)
  mi = min(comp_covar)
  ma = max(comp_covar)
  comp_covar = (comp_covar - mi)*2/(ma - mi)
  PVAL_score_grids = apply(comp_covar, 2,function(z)ScoreTest_SPA(genos=z,obj.null=nullmod$nullmod_spa)$p.value)
  PVAL_lrt_grids = apply(comp_covar, 2,function(z)anova(nullmod$nullmod_lrt , glm(y~nullcovar+ z, family="binomial"), test="LRT")[2,5])
  PVAL_grids = c(PVAL_score_grids, PVAL_lrt_grids)
  cauchy_stat = sum(tan(pi*(0.5-PVAL_grids)), na.rm=T)/(sum(!is.na(PVAL_grids)))
  cauchy_stat_lrt = sum(tan(pi*(0.5-PVAL_lrt_grids)), na.rm=T)/(sum(!is.na(PVAL_lrt_grids)))
  cauchy_stat_score = sum(tan(pi*(0.5-PVAL_score_grids)), na.rm=T)/(sum(!is.na(PVAL_score_grids)))
  p_cwot = pcauchy(cauchy_stat, lower.tail = FALSE)
  p_cwot_lrt = pcauchy(cauchy_stat_lrt, lower.tail = FALSE)
  p_cwot_score = pcauchy(cauchy_stat_score, lower.tail = FALSE)
  return(list(p_cwot = p_cwot, p_cwot_lrt = p_cwot_lrt, p_cwot_score = p_cwot_score,
              PVAL_lrt_grids = PVAL_lrt_grids, PVAL_score_grids = PVAL_score_grids))
}

