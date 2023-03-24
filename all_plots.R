### function that runs the test and then does the plots ####


all_plots <- function(test){
  #currently this does converge to values of a similar value to what we want but once we get there it doesn't move at all 
  plot_y_trace(test)
  
  #should be 1
  plot(test$beta, type = "l" )
  
  ## looking for 10 in left and 0.5 in right 
  plot_theta_trace(test)

  print(test$acceptance_rate)
}