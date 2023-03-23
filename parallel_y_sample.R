### method 2 - splitting up our spatial area into subsets #### 


parallel_y_sample <- function(M, K, label_vec, abc_mat, abc_delta, beta, theta, full_prior_y){
  
  for (m in 1:M){
    ##filter data into only points in the mth group 
    prior_y_group <- full_prior_y[which(label_vec== m )]
    
    ## for each subgroup k in group m we want a list of the neighbouring subgroups from other groups so we know to include these points in the respective Sigma, so we will be taking a subset of the distance matrix ?? then using this to re construct our Sigma 
    
    abc_prop <- ### want a way to save this and then put it all back together so can keep it as the prior for the next round
        
    sample_list  <- foreach(k=1:K) %dopar% {  
      ### propose new y value for subset k of group m (in M)
      ### or work all of them out and subset ?? 
      prior_y <- prior_y_group(## which subgroup = k )
      
      ### prep matrices 
      Sigma = theta[2]* exp(-1 * theta[1]* dist_mat[### subset by columns we need ])
  
      abc_prior <- abc_mat[#subset columns that are in our subgroup ]
      b_prior <- exp(beta) * t(t_cov*abc_prior[,2])
      c_prior <- exp(beta)* diag(t_cov*abc_prior[,3])
        
      inv_Sigma_C <- chol2inv(chol( chol2inv(chol(Sigma)) + c_prior ))
      
      prop_sigma <-  inv_Sigma_C
      prop_mu = prop_sigma %*%t(x_sub-b_prior) 
      
      ### Z from current y values   
      Z_prior = (x_sub- b_prior) %*% inv_Sigma_C %*% t(x_sub -b_prior)
      
      ### sample y 
      prop_y <- rmvnorm(1, mean = prop_mu, sigma = prop_sigma )  
      
      
      ### re calc abc, z, det etc with proposal 
      
      ### Z from current y values  
      abc_prop <-   find_abc(prop_y, delta = abc_delta ) #y at this point already correct dimensions 
      b_prop <- exp(beta_sample[t-1]) * t(t_cov*abc_prop[,2])
      c_prop <- exp(beta_sample[t-1])* diag(t_cov*abc_prop[,3])      
      
      inv_temp_Sigma_C_prop <- chol2inv(chol(inv_temp_Sigma+ c_prop))
      Z_prop = (x_sub - b_prop)%*% inv_temp_Sigma_C_prop %*%  t(x_sub - b_prop)
      
      C_list <- list(c_prior, c_prime)
      
      list_det_Sigma_C <-foreach(i=1:2) %dopar% {  
        prod(diag(chol(inv_temp_Sigma + C_list[[i]]))^2) #Determinant of Sigma+{C,C'} where C' uses proposed y values
      }
      
      ### acceptance probability for proposal of y 
      acc_prob_y <- 0.5*prop_y%*%c_prior%*%t(prop_y) + b_prior%*%t(prop_y) -0.5*prior_y%*%c_prop%*%t(prior_y) - b_prop%*%t(prior_y) - t(t_cov)%*%t(exp(beta + prop_y)) + t(t_cov)%*%t(exp(beta + prior_y)) -0.5* log(list_det_Sigma_C[[1]]) + 0.5 * log(list_det_Sigma_C[[2]]) + 0.5*Z_prior - 0.5*Z_prop
      
      
      u <- runif(1,0,1)
      if (acc_prob_y >= log(u) ){
        y_sample[t,] <- prop_y
        prior_y <- prop_y 
        print(t)
        acc_counts[1] <- acc_counts[1]+1 
        ABC <- ABC_prime
      } else {
        y_sample[t,] <- prior_y
      }
      
    }
      
    

  }
  
  
}

M = 
K = detectCores(logical=FALSE)




