getCOV = function(Tr, G, P0, rho_grids, M1=NULL, M2=NULL, M3=NULL, M4=NULL){

  G = as.matrix(G);
  Tr = as.matrix(Tr);
  GT = G*Tr
  tmp1 = P0%*%G
  tmp2 = P0%*%GT
  v1 = c(t(G)%*%tmp1)
  v2 = c(t(GT)%*%tmp2)
  v3 = c(t(G)%*%tmp2)

  rho_grids = as.matrix(rho_grids)
  one = as.matrix(rep(1,length(rho_grids)),ncol=1)
  if(is.null(M1)){
    M1 = rho_grids%*%t(rho_grids)
  }
  if(is.null(M2)){
    M2 = (1-abs(rho_grids))%*%t(1-abs(rho_grids))
  }
  if(is.null(M3)){
    M3 = rho_grids%*%t(one)
  }
  if(is.null(M4)){
    M4 = rho_grids%*%t(abs(rho_grids))
  }
  return(M1*v1 + M2*v2 + (M3+t(M3)-M4-t(M4))*v3)
}

# getWeights <- function(n, a){
#   id = 1:n
#   (abs(2*id-n-1)/(n-1))^a*sign(2*id-n-1)
# }
#
# getWeights_v1 <- function(n, a){
#   n = n - 1
#   f2 = function(x)2^(2*a-1)*(abs(x/n-0.5))^a*sign(x/n-0.5)
#   f3 = function(x)1-f2(1.5*n-x)
#   f1 = function(x)-1+f2(0.5*n+x)
#   f = function(x, n) f1(x)*(x<n/4) + f2(x)*(x>=0.25*n)*(x<=0.75*n) + f3(x)*(x>0.75*n)
#   id = 0:n
#   f(id, n)
# }
