### function to simulate the data for different sizes under the same process #### 



sim_data <- function(len, range = 10 , sigma2= 0.5, beta = 1, min_d = 0 , max_d = 1/sqrt(2)){
  # Simulate the data
  
  xycoords <- seq(from=min_d,to= max_d ,length=len)
  allcoords <- expand.grid(xycoords,xycoords)
  nlocs <- nrow(allcoords)
  
  # make the covariance function for the GPs
  true.range <- range ## this is what they call alpha 
  true.sigma2 <- sigma2
  distmat <- as.matrix(dist(allcoords,diag=TRUE,upper=TRUE))
  Sigma <- true.sigma2 * exp(-true.range*distmat)
  
  # simulate the GP at each location
  ys <- chol(Sigma)%*%rnorm(nlocs) # leave it at mean=0
  
  # simulate poisson process based on that
  
  #ts <- runif(nlocs,min=1,max=10) # observation durations
  ts <- rep(100, nlocs)
  beta <- beta # will not be differentiated from the mean of the GP,so use mean=0
  lambdas <- ts*exp(ys+beta)
  
  Xs <- rpois(nlocs,lambdas)
  
  # # plot heat map of the data
  # plotdata <- data.frame(allcoords,Xs,ys)
  # ggplot(plotdata, aes(Var1, Var2, fill= Xs)) + 
  #   geom_tile()
  # 
  # # we want to model y, etc from above data to compare to true surface
  # ggplot(plotdata, aes(Var1, Var2, fill= ys)) + 
  #   geom_tile()
  # 
  # allcoords$id <- 1:length(allcoords[,1])
  # plot(allcoords[1:2])
  # text(allcoords[1:2], labels = allcoords[,3])
  # 
  # 
  
  return(list("ys" = ys, "Xs" = Xs , "ts" = ts, "distance" = distmat))
}