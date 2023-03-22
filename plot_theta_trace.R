### function to plot theta ssample side by side as lines ### 

plot_theta_trace <- function(object){
  num_th <- dim(object$theta)[2]
  theta <- object$theta
  par(mfrow=c(1,num_th))
  title = c("alpha", "sigma2")
  for (i in 1:num_th){
    plot(theta[,i], type = "l", main = title[i] )
  }
  par(mfrow=c(1,1))
}
