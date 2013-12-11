
# --------------------------------------------------------------
#make_lower_bound function returns the lower_bound function
#lower_bound function is set to be -Inf when argument falls outside the abscissaes range.

make_lower_bound <- function(x, hx) {
  min <- min(x)
  max <- max(x)
  left <- min - abs(min*.Machine$double.eps)
  right <- max + abs(max*.Machine$double.eps)
  x <- c(left, x, right)
  hx <- c(-Inf, hx, -Inf)
  lower_bound <- approxfun(x, hx, rule = c(2,2))
  return(lower_bound)
}

make_upper_bound <- function(x, hx, hpx, z) {
  upper_bound_eval_z <- function(z) {
    index <- rowSums(outer(z, z, function(x1,x2) x1>x2))
    index[index==0] <- 1
    index[index==length(z)] <- length(z) - 1
    ux <- hx[index] + (z-x[index])*hpx[index]
    return (ux)
  }
  # Because of numerical issues
  # Need to find a better way to do to this. Is x sorted? 
  m <- length(z)
  #if (z[1] == -Inf) z[1] <- -1000*max(abs(x))
  if (z[1] == -Inf) z[1] <- -1e300
  #if (z[m] == Inf) z[m] <- 1000*max(abs(x))
  if (z[m] == Inf) z[m] <- 1e300
  y <- upper_bound_eval_z(z)
  upper_bound <- approxfun(z,y)
  return(upper_bound)
}

#------------------------------------------------
