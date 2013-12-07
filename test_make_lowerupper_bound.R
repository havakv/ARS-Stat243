# Testing make_lower_bound and make_upper_bound

source("ARS.R")

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


#################################################################

x <- seq(-10.5, 10.5, by=4)
hx <- -x^2
hpx <- -2*x
z <- make_z(x,hx,hpx,-12.0, 12.0)
lower <- make_lower_bound(x, hx, -12.0, 12.0)
lower2 <- make_lower_bound2(x, hx, -12.0, 12.0)

xx <- runif(n = 1e5, min = -12, max = 12)

library(rbenchmark)
benchmark(lower(xx), lower2(xx),
	  replications = 10, columns=c('test', 'elapsed', 'replications'))


