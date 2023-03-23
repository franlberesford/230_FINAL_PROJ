### function that runs the test and then does the plots ####


alg1_and_plots <- function(num_locs=nlocs, num_its= 1000, dist_mat = distmat, X=Xs, t_cov= ts, abc_delta= 1e-30,beta_prop_sigma = 0.01, theta_prop_sigma = c(1,5)){
  test <- mcmc_alg_meth_1(num_locs, num_its= num_its, dist_mat = distmat, X=Xs, t_cov= ts, abc_delta,beta_prop_sigma , theta_prop_sigma)
  
  
  #currently this does converge to values of a similar value to what we want but once we get there it doesn't move at all 
  plot_y_trace(test)
  
  #should be 1
  plot(test$beta, type = "l" )
  
  ## looking for 10 in left and 0.5 in right 
  plot_theta_trace(test)

  print(test$acceptance_rate)
}