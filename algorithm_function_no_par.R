### algorithm function ### 

mcmc_alg_meth_1 <- function(num_locs, num_its, dist_mat, X, t_cov, beta_prop_sigma, theta_prop_sigma){
  
    #length of covariance parameter vector 
    length_theta <- length(theta_prop_sigma)
    
    ## set up matrices to store chains 
    y_sample <- matrix(NA, ncol = num_locs, nrow = num_its)
    theta_sample <- matrix(NA, ncol = length_theta, nrow = num_its)
    beta_sample <- numeric(length = num_its) #i think this is just one number so we only need it to be a vector but i could be wrong ? 
    
    #initialise sample 
    y_sample[1,] <- y_init <- rnorm(nlocs, 0, 25)
    theta_sample[1,] <- theta_init <-  rgamma(length_theta, shape=9 ,scale = 0.5)
    beta_sample[1] <- beta_init <-  rnorm(1, 0, 25)
    
    #theta proposal memory allocation 
    theta_prop <- numeric(length = length_theta)
  
  for (t in 2:num_its){
    #### y
    
    ### prep matrices 
    temp_Sigma = theta_sample[(t-1),2]* exp(-1 * theta_sample[(t-1),1]* dist_mat) #Covariance matrix using theta values from previous iteration. 
    #temp_Sigma 
    
    ## make prior need to update Sigma matrix 
    if (t==2){
      ABC <- find_abc(y_sample[t-1,], delta =0.2)
      C <- exp(beta[t-1])* diag(t_cov*ABC[,3])
      B <- exp(beta[t-1]) * t(t_cov*ABC[,2])
    }
    
    ### matrices to parallise
    Sigma_R <- chol(temp_Sigma) #this gives upper triangular but paper uses lower triangular 
    inv_temp_Sigma <- chol2inv(Sigma_R)
    #yt_temp_Sigma_y <- t(y_sample[t-1,])%*%solve(temp_Sigma)%*%y_sample[t-1,] #not needed at this step i don't think 
    det_Sigma <- prod(diag(Sigma_R)^2)
    det_sigma_C <-  prod(diag(chol(inv_temp_Sigma+ C))^2) #sigma using theta values from time t-1 
    
    
    #prior_Sigma <-  solve((solve(temp_Sigma) + C))
    #prior_mu <- prop_sigma %*%t(X-B)
    
    
    ### setting up parameters for y proposal 
    Sigma_C_R <- chol(inv_temp_Sigma + C ) 
    inv_Sigma_C <- chol2inv(Sigma_C_R)
    prop_sigma <- inv_Sigma_C  #sigma uses thetas from previous iteration 
    prop_mu = prop_sigma %*%t(X-B) 
    
    ### Z from current y values   
    Z = (X-B) %*% inv_Sigma_C %*% t(X-B)
    
    ### sample y 
    prop_y <- rmvnorm(1, mean = prop_mu, sigma = prop_sigma )
    
    ### calculate new ABC values as well as Z and determinant 
    ABC_prime <- find_abc(prop_y, delta = 0.2)
    C_prime <- exp(beta[t-1])* diag(t_cov*ABC_prime[,3])
    B_prime <- exp(beta[t-1] ) * t(t_cov*ABC_prime[,2])
    inv_temp_Sigma_C_prime <- chol2inv(chol(inv_temp_Sigma+ C_prime))
    Z_prime = t(X-B_prime) %*% inv_temp_Sigma_C_prime %*% (x-B_prime)  
    det_sigma_C_prime <-  prod(diag(chol(inv_temp_Sigma + C_prime))^2) #Determinant of Sigma+C' where C' uses proposed y values 
    
    ##(paper computes the two determinants Sigma+C and Sigma+C' parallelly but i dont see how that makes sense given the order of the code so maybe something is wrong here )
    
    
    ### acceptance probability for proposal of y 
    acc_prob_y <- 0.5*t(prop_y)%*%C%*%prop_y + t(B)%*%prop_y -0.5*t(prior_y)%*%C%*%prior_y - t(B_prime)%*%prior_y - t(t)%*%exp(beta[t-1] + prop_y) 
    + t(t_cov)%*%exp(beta[t-1] + prior_y) -0.5* log(det_sigma_C) + 0.5 * log(det_sigma_C_prime) + 0.5*Z - 0.5*Z_prime 
    
    
    u <- runif(1,0,1)
    if (acc_prob_y >= log(u) ){
      y_sample[t,] <- prop_y
      A <- ABC_prime[,1]
      B <- B_prime
      C <- C_prime 
    } else {
      y_sample[t,] <- prior_y
    }
    
    #### beta 
    #(note that we do not need any of the letters for this proposal )
    if (t==2) {
      prior_beta <- beta_sample[t-1]
    }
    prop_beta <- rnorm(1, prior_beta, beta_prop_sigma)
    
    acc_prob_beta <- exp(prior_beta - prop_beta) * t(t_cov)%*% exp(y[t,]) + (prior_beta - prop_beta)* t(log(t_cov))%*% X
    
    u <- runif(1,0,1) 
    if (acc_prob_beta >= log(u)){
      beta_sample[t] <- prop_beta
      prior_beta <- prop_beta
    } else {
      beta[t] <- prior_beta 
    }
    
    
    #### theta 
    #(theta based on 2 independent random walk proposals and I've forgotten what they are )
    
    #propose new theta value
    for (i in 1:length(theta_init)){
      theta_prop[i] <- rnorm(1, mean = theta_sample[t-1, i], sd = theta_prop_sd[i])
    }
    
    #prior and proposed Sigma matrices 
    temp_Sigma_prime = theta_prop[2]* exp(-1 * theta_prop[1]* dist_mat)
    temp_Sigma_prior = theta_sample[t-1,2]* exp(-1 * theta_sample[t-1,1]* dist_mat)
    
    #calculate determinants 
    det_temp_Sigma_prime <-  prod(trace(chol(temp_Sigma_prime))^2)
    det_temp_Sigma_prior <-  prod(trace(chol(temp_Sigma_prior))^2)
    
    
    #calculate joint acceptance probability 
    acc_prob_theta <- 0.5*(-t(y_sample[t,])%*%solve(temp_Sigma_prime)%*%y_sample[t,] +  t(y_sample[t,])%*%solve(temp_Sigma_prior)%*%y_sample[t,] - log( det_temp_Sigma_prime) + log(det_temp_Sigma_prop ) ) 
    
    u <- runif(1,0,1) 
    if (acc_prob_theta >= log(u)){
      theta_sample[t,] <- prop_theta
      prior_theta <- prop_theta
    } else {
      theta_sample[t-1,] <- prior_theta 
    }
  }
  
  return(list("ys" = y_sample, "beta" = beta_sample, "theta" = theta_sample))
}
