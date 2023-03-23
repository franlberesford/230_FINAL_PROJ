#### funnction that calculates determinant of a matrix + inverse sigma matrix ####

det_func <- function(x, inv_Sigma){ prod(diag(chol(inv_Sigma + x))^2) }   

