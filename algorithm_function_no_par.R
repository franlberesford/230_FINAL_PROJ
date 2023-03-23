### algorithm function ### 

mcmc_alg_meth_0 <- function(num_locs, num_its, dist_mat, X, t_cov, abc_delta,  beta_prop_sigma, theta_prop_sigma){
  
    acc_counts <- rep(0,3)
  
    #length of covariance parameter vector 
    length_theta <- length(theta_prop_sigma)
    
    ## set up matrices to store chains 
    y_sample <- matrix(NA, ncol = num_locs, nrow = num_its)
    theta_sample <- matrix(NA, ncol = length_theta, nrow = num_its)
    beta_sample <- numeric(length = num_its) #i think this is just one number so we only need it to be a vector but i could be wrong ? 
    
    #initialise sample 
    theta_sample[1,] <- theta_init <- c(rgamma(1, 9, 2) ,rinvgamma(1, 0.5, 0.5))
    y_sample[1,] <- y_init <- rmvnorm(1, mean = rep(0,num_locs),  sigma = theta_init[2]* exp(-1 * theta_init[1]* dist_mat) )
    beta_sample[1] <- beta_init <-  rnorm(1, mean(X/t_cov), sd(X/t_cov))
    
    #theta proposal memory allocation 
    theta_prop <- numeric(length = length_theta)
  
  for (t in 2:num_its){
    #print(t)
    #### y
    
    ### prep matrices 
    Sigma = theta_sample[(t-1),2]* exp(-1 * theta_sample[(t-1),1]* dist_mat) #Covariance matrix using theta values from previous iteration. 
    #temp_Sigma 
    
    ## make prior need to update Sigma matrix 
    if (t==2){
      prior_y <- t(matrix(y_init))
      ABC <- find_abc(y_sample[t-1,], delta =abc_delta )
    }
    
    ### recalculate B an C based on new beta values 
    C <- exp(beta_sample[t-1])* diag(t_cov*ABC[,3])
    B <- exp(beta_sample[t-1]) * t(t_cov*ABC[,2])
    
    #print(ABC[,3])
    
    
    ### matrices to parallise
    Sigma_R <- chol(Sigma) #this gives upper triangular but paper uses lower triangular 
    inv_Sigma <- chol2inv(Sigma_R)
    #yt_temp_Sigma_y <- t(y_sample[t-1,])%*%solve(temp_Sigma)%*%y_sample[t-1,] #not needed at this step i don't think 
    det_Sigma <- prod(diag(Sigma_R)^2)
    
    #print(inv_temp_Sigma)
    #print(inv_temp_Sigma+C)
    
    det_sigma_C <-  prod(diag(chol(inv_Sigma+ C))^2) #sigma using theta values from time t-1 
    
    
    #prior_Sigma <-  solve((solve(temp_Sigma) + C))
    #prior_mu <- prop_sigma %*%t(X-B)
    
    
    ### setting up parameters for y proposal 
    Sigma_C_R <- chol(inv_Sigma + C ) 
    inv_Sigma_C <- chol2inv(Sigma_C_R)
    prop_sigma <- inv_Sigma_C  #sigma uses thetas from previous iteration 
    prop_mu = prop_sigma %*%t(X-B) 
    
    
    ### Z from current y values   
    Z = (X-B) %*% inv_Sigma_C %*% t(X-B)
    
    ### sample y 
    prop_y <- rmvnorm(1, mean = prop_mu, sigma = prop_sigma )
    
    #print(prop_y)
    
    ### calculate new ABC values as well as Z and determinant 
    ABC_prime <- find_abc(prop_y, delta = abc_delta )
    C_prime <- exp(beta_sample[t-1])* diag(t_cov*ABC_prime[,3])
    B_prime <- exp(beta_sample[t-1] ) * t(t_cov*ABC_prime[,2])
    #print(ABC_prime[,3])
    #print(inv_temp_Sigma + C_prime)
    
    inv_Sigma_C_prime <- chol2inv(chol(inv_Sigma+ C_prime))
    Z_prime = (X-B_prime) %*% inv_Sigma_C_prime %*% t(X-B_prime)  
    det_sigma_C_prime <-  prod(diag(chol(inv_Sigma + C_prime))^2) #Determinant of Sigma+C' where C' uses proposed y values 
    
    ##(paper computes the two determinants Sigma+C and Sigma+C' parallelly but i dont see how that makes sense given the order of the code so maybe something is wrong here )
    
    
    
    ### acceptance probability for proposal of y 
    acc_prob_y <- 0.5*prop_y%*%C%*%t(prop_y) + B%*%t(prop_y) -0.5*prior_y%*%C_prime%*%t(prior_y) - B_prime%*%t(prior_y) - t(t_cov)%*%t(exp(beta_sample[t-1] + prop_y)) + t(t_cov)%*%t(exp(beta_sample[t-1] + prior_y)) -0.5* log(det_sigma_C) + 0.5 * log(det_sigma_C_prime) + 0.5*Z - 0.5*Z_prime 
    
    
    u <- runif(1,0,1)
    if (is.na(acc_prob_y >= log(u) ) ){
      y_sample[t,] <- prior_y
    } else if(acc_prob_y >= log(u)) {
      y_sample[t,] <- prop_y
      prior_y <- prop_y 
      print(t)
      acc_counts[1] <- acc_counts[1]+1 
      ABC <- ABC_prime
    } else {
      y_sample[t,] <- prior_y
    }
    
    #### beta 
    #(note that we do not need any of the letters for this proposal )
    if (t==2) {
      prior_beta <- beta_sample[t-1]
    }
    prop_beta <- rnorm(1, prior_beta, beta_prop_sigma)
    #print(prop_beta)
    
    acc_prob_beta <- (exp(prior_beta) - exp(prop_beta) )* t(t_cov)%*% exp(y_sample[t,]) + (prop_beta - prior_beta)* t(log(t_cov))%*% X
    
    u <- runif(1,0,1) 
    if (acc_prob_beta >= log(u)){
      acc_counts[2] <- acc_counts[2] + 1 
      beta_sample[t] <- prop_beta
      prior_beta <- prop_beta
    } else {
      beta_sample[t] <- prior_beta 
    }
    
    
    #### theta 
    #(theta based on 2 independent random walk proposals and I've forgotten what they are )
    if (t==2){
      prior_theta <- theta_init 
    }
    
    
    #propose new theta value
    for (i in 1:length(theta_init)){
      theta_prop[i] <- rtruncnorm(1, a=0, b=Inf , mean = theta_sample[t-1, i], sd = theta_prop_sigma[i])
    }
    
    #prior and proposed Sigma matrices 
    Sigma_prime = theta_prop[2]* exp(-1 * theta_prop[1]* dist_mat)
    Sigma_prior = theta_sample[t-1,2]* exp(-1 * theta_sample[t-1,1]* dist_mat)
    
    #calculate determinants 
    det_Sigma_prime <-  prod(diag(chol(Sigma_prime))^2)
    det_Sigma_prior <-  prod(diag(chol(Sigma_prior))^2)
    
    #inverse of the Sigma matrices 
    inv_Sigma_prime <- chol2inv(chol(Sigma_prime))
    inv_Sigma_prior <- chol2inv(chol(Sigma_prior))
    
    #calculate joint acceptance probability 
    acc_prob_theta <- 0.5*(-t(y_sample[t,])%*%inv_Sigma_prime%*%y_sample[t,] +  t(y_sample[t,])%*%inv_Sigma_prior%*%y_sample[t,] - log( det_Sigma_prime) + log(det_Sigma_prior ) ) 
    
    u <- runif(1,0,1) 
    if (is.na(acc_prob_theta >= log(u))){
      print(acc_prob_theta)
      print(paste("y^TSigma^-1 prime y = ",t(y_sample[t,])%*%inv_Sigma_prime%*%y_sample[t,]  ))
      print(paste("y^TSigma^-1 prior y = ",t(y_sample[t,])%*%inv_Sigma_prior%*%y_sample[t,]  ))
      print(paste("log det Sigma prime = ",log( det_Sigma_prime)  ))
      print(paste("log det Sigma prior = ",log( det_Sigma_prior)  ))
      theta_sample[t,] <- prior_theta 
      
      #re define the acceptance probability ignoring the -Infs 
      acc_prob_theta <-  0.5*(-t(y_sample[t,])%*%inv_Sigma_prime%*%y_sample[t,] +  t(y_sample[t,])%*%inv_Sigma_prior%*%y_sample[t,]) 
    } else if (acc_prob_theta >= log(u)) {
      acc_counts[3] <- acc_counts[3] + 1 
      theta_sample[t,] <- theta_prop
      prior_theta <- theta_prop
    } else {
      theta_sample[t,] <- prior_theta  
    }
   if (floor(t/1000) == t/1000 ){print(paste0(t, " out of ", num_its, " iterations done."))}
  }
  acc_rates <- acc_counts / num_its
  
  return(list("ys" = y_sample, "beta" = beta_sample, "theta" = theta_sample, "acceptance_rate" = acc_rates))
}
