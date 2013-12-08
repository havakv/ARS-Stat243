#update_x function adds the update point into the abscissaes vector, and also checks if the derivative at the update point is between its neighbours(i.e., checking the cancaveness)

update_x <- function(f, x, hx, hpx, update) {
  if (is.na(update)){
    return (list(x=x, hx=hx, hpx=hpx))
  }
  else {
    hx_update <- log(f(update))
    xx <- update
    hpx_update <- attributes(numericDeriv(quote(log(f(xx))), "xx"))$gradient[1, 1]
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

##Test case for update_x###
#f <- function(x) return(x^2)
#x <- seq(1.5, 10.5, by=1)
#hx <- log(f(x))
#hpx <- 2/x

#update_x(f, x, hx, hpx, 5.1)
#update_x(f, x, hx, hpx, 12)
#update_x(f, x, hx, hpx, 0.5)

#ff <- function(x) return(exp(x^2))
#update_x(ff, x, hx, hpx, 5.1)
