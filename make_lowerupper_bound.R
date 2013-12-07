#make_lower_bound function returns the lower_bound function
#lower_bound function is set to be -Inf when argument falls outside the abscissaes range.
make_lower_bound <- function(x, hx, left_bound, right_bound) {
  x <- c(left_bound, x, right_bound)
  hx <- c(-Inf, hx, -Inf)
  lower_bound <- function(xx) {
    k=length(x)
    index <- rowSums(outer(xx, x, function(x1,x2) x1>x2))
    index[index==0] <- 1
    index[index==k] <- k - 1
    lx <- ((x[index+1]-xx)*hx[index] + (xx-x[index])*hx[index+1])/(x[index+1]-x[index])
    lx[index %in% c(1,k-1)] <- -Inf
    return(lx)
  }
  return(lower_bound)
}

#make_upper_bound function returns the upper_bound function
make_upper_bound <- function(x, hx, hpx, z) {
  upper_bound <- function(xx) {
    index <- rowSums(outer(xx, z, function(x1,x2) x1>x2))
    index[index==0] <- 1
    index[index==length(z)] <- length(z) - 1
    ux <- hx[index] + (xx-x[index])*hpx[index]
    return (ux)
  }
  return(upper_bound)
}

#################################################################
make_lower_bound2 <- function(x, hx, left_bound, right_bound) {
  x <- c(left_bound, x, right_bound)
  hx <- c(-Inf, hx, -Inf)
  lower_bound <- approxfun(x,hx)
  return(lower_bound)
}

#------------------------------------------------
