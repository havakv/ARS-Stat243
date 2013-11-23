make_z <- function(x, fx, fpx, left_bound, right_bound) {
  k <- length(x)
  xm1 <- c(0, x)
  fxm1 <- c(0, fx)
  fpxm1 <- c(0, fpx)
  x <- c(x, 0)
  fx <- c(fx, 0)
  fpx <- c(fpx, 0)
  
  z <- (fx - fxm1 - x*fpx + xm1*fpxm1) / (fpxm1 - fpx)
  z[1] <- left_bound
  z[k+1] <- right_bound
  
  return (z)
}

#Test case for make_z
# x <- seq(-10, 10, by = 2)
# fx <- x^2
# fpx <- 2*x
# plot(x, fx)
# points(make_z(x, fx, fpx, -12, 12), make_z(x, fx, fpx, -12, 12)^2, col="red")