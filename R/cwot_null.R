#' Prepare null model for cwot.
#' @param y - Binary response variable.
#' @param tr - Binary treatment variable, 1 for treated, 0 for placebo.
#' @param x - Covariates in addition to treatment variable.
#' @return A list of objects needed for cwot.
#' @references Hong Zhang, Devan Mehrotra and Judong Shen, "AWOT and CWOT for Genotype and Genotype by Treatment Interaction Joint Analysis in Pharmacogenetics GWAS".
#' @examples
#' n = 100
#' y = rbinom(n, 1, 0.3)
#' x = data.frame(x1=rnorm(n))
#' tr = rbinom(n, 1, 0.5)
#' nullmod = cwot_null(y, tr, x)
#' @export
#' @import stats SPAtest

cwot_null = function(y, tr, x=NULL){
  if(is.null(x)){
    nullcovar = data.frame(tr=tr)
  }else{
    nullcovar = data.frame(tr=tr,x=x)
  }
  nullmod_spa = ScoreTest_wSaddleApprox_NULL_Model(y~., data=nullcovar)
  nullmod_lrt = glm(y~., data=nullcovar, family="binomial")
  return(list(nullmod_spa=nullmod_spa,
              nullmod_lrt=nullmod_lrt,
              nullcovar = nullcovar,
              y=y))
}
