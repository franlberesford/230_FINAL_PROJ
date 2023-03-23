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
    e_y_d <- exp(y_p)-exp(y_m)
    
    t1 <- 2*exp(y_p)*(y_p-1) - 2*exp(y_m)*(y_m-1)
    t2 <- -exp(y_p)*(y_p^2-2*y_p+2) + exp(y_m)*(y_m^2-2*y_m+2)
    t3 <- y_d_3/(2*y_d)
    t4 <- y_d_2/(2*y_d)
    t5 <- (2/3)*y_d_3 - (y_d_2^2)/(2*y_d)
    t6 <- ((1/3)*y_d_3*t4/t5 - (1/(4*t5)) * y_d_4)
    t7 <- 2*e_y_d*(t4/t5) + (t1/t5)
    t8 <- (-t3/9)*y_d_3 - t6*y_d_2*t3/3 + (1/4)*t6*y_d_4 + (1/10)*y_d_5
    t9 <- -t7*y_d_2*t3/3 - (2/3)*e_y_d*t3 + (t7/4)*y_d_4 - t2
    
    c_vec[i] <- max(t9/t8,0)
    b_vec[i] <- c_vec[i]*t6 - t7
    a_vec[i] <- 1
    
    #opts <- optim(par=c(1,1,1),fn=fntomin,method="SANN")
    
    
    #a_vec[i] <- opts$par[1]#ifelse(is.na(a), runif(1, -1, 1), a)
    #b_vec[i] <- opts$par[2]#ifelse(is.na(b), runif(1, -1, 1), b)
    #c_vec[i] <- opts$par[3]#ifelse(is.na(c), rgamma(1, 9,2), ifelse(c<0 , runif(1, 1e-15, 1e-10), c))
  }
  
  return(matrix(c(a_vec,b_vec,c_vec), ncol = 3))
  
  
  
  
  
  }
