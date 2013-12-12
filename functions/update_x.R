
# ----------------------------------------------------------------------------
#update_x function adds the update point into the abscissas vector, and also checks if the derivative at the update point is between its neighbors(i.e., checking the concaveness)

update_x <- function(f, x, hx, hpx, update, ...) {
  if (is.na(update)){
    return (list(x=x, hx=hx, hpx=hpx))
  }
  else {
    hx_update <- log(f(update, ...))
    xx <- update
    hpx_update <- attributes(numericDeriv(quote(log(f(xx, ...))), "xx"))$gradient[1, 1]
    rank_x <- rank(c(x, update))
    new_x <- rep(NA, length(x)+1)
    new_hx <- rep(NA, length(x)+1)
    new_hpx <- rep(NA, length(x)+1)
    
    new_x[rank_x] <- c(x, update)
    new_hx[rank_x] <- c(hx, hx_update)
    new_hpx[rank_x] <- c(hpx, hpx_update)
    
    #Check concaveness by testing if h'(x) is decreasing
    if (is.unsorted(-new_hpx)) stop("The sample function is not log-concaved!")
    
    return(list(x=new_x, hx=new_hx, hpx=new_hpx))
  }
}

