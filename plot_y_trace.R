### function to plot y sample as lines ### 

plot_y_trace <- function(object){
  num_y <- dim(object$ys)[2]
  ys <- object$ys 
  range_y <- range(ys)
  plot(ys[,1], col = rand_color(1), type = "l", ylim = range_y)
  for (i in 2:num_y){
    lines(ys[,i], col = rand_color(1))
  }
}
