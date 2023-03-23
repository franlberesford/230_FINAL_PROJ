
###function to use for parallel loop in theta proposal part of MCMC #### 

par_functions_theta_prop <- function(ind, y, Sigma){
  if (ind ==1 ){
    out <- t(y)%*%chol2inv(chol(Sigma))%*%y
  } else if (ind ==2){
    out <-  chol2inv(chol(Sigma))
  } else if (ind ==3){
    out <-  prod(diag(chol(Sigma))^2)
  }
  return(out)
}

