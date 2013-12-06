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

#Test case for make_lower_bound function
x <- seq(-10.5, 10.5, by=4)
hx <- -x^2
hpx <- -2*x
z <- make_z(x,hx,hpx,-12.0, 12.0)
lower <- make_lower_bound(x, hx, -12.0, 12.0)
upper <- make_upper_bound(x, hx, hpx, z)
xx <- seq(-12, 12, by=0.2)
plot(xx, -xx^2, type='l')
lines(xx, lower(xx), col="red")
lines(xx, upper(xx), col="blue")

x <- seq(-10.5, 10.5, by=4)
hx <- -x^2
hpx <- -2*x
z <- make_z(x,hx,hpx,-Inf, Inf)
lower <- make_lower_bound(x, hx, -Inf, Inf)
upper <- make_upper_bound(x, hx, hpx, z)
xx <- seq(-20, 20, by=0.2)
plot(xx, -xx^2, type='l')
lines(xx, lower(xx), col="red")
lines(xx, upper(xx), col="blue")



