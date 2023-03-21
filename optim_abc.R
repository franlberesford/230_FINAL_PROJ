### find ABC function ### 


## returns vector format and then will transform in the code 

find_abc <- function(y, delta){
  a_vec <- b_vec <- c_vec <- numeric(length = length(y))
  
  for (i in 1:length(y)){
    y_p <- y[i] + delta 
    y_m <- y[i] - delta 
    
    y_p_2 <- y_p^2
    y_p_3 <- y_p^3
    y_p_4 <- y_p^4
    y_p_5 <- y_p^5
    
    y_m_2 <- y_m^2
    y_m_3 <- y_m^3
    y_m_4 <- y_m^4
    y_m_5 <- y_m^5
    
    y_d <- y_p - y_m
    y_d_2 <- y_p_2 - y_m_2 
    y_d_3 <- y_p_3 - y_m_3 
    y_d_4 <- y_p_4 - y_m_4 
    y_d_5 <- y_p_5 - y_m_5
    
    
    t1 <- 8/3* y_d_3/y_d_4 *(exp(y_p)*y_p - exp(y_m)*y_m )
    t2 <- 2*y_d_4
    t3 <- 2*(exp(y_p)*(y_p-1)^2 - exp(y_m)*(y_m -1) )
    t4 <- 4/y_d_4 *(exp(y_p)*y_p - exp(y_m) *y_m)
    
    b <- (1/(6*delta)* (t1+t2) -t3 - 2/5*y_d_5*t4) / (-1/(6*delta) *(8/9*y_d_2 + 4/3 *y_d + y_d_2  ) - y_d_4/2 + 2/5 *y_d_5 *(4/3* 1/y_d + 2*y_d_2/ y_d_4))
    
    a <- -1/(4*delta) *(8/9 * b * y_d_3^2 / y_d_4 + 8/3* y_d_3/y_d_4 *(exp(y_p)*y_p - exp(y_m)*y_m) + 4/3 *b * y_d_2*y_d_3 / y_d_4 + 2*y_d_4+ b * y_d_2     ) 
    
    c <-  4/3 * b * y_d_3/y_d_4 + 4 / y_d_4 *(exp(y_p)*y_p - exp(y_m) * y_m) + 2 * b *y_d_2 / y_d_4  
    
    a_vec[i] <- ifelse(is.na(a), runif(1, -1, 1), a)
    b_vec[i] <- ifelse(is.na(b), runif(1, -1, 1), b)
    c_vec[i] <- ifelse(is.na(c), rgamma(1, 9,2), ifelse(c<0 , runif(1, 1e-15, 1e-10), c))
  }
  
  return(matrix(c(a_vec,b_vec,c_vec), ncol = 3))
}
