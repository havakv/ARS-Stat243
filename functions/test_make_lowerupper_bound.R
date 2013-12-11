# Testing make_lower_bound and make_upper_bound

source("../ARS.R")

x <- seq(-10.5, 10.5, by=4)
hx <- -x^2
hpx <- -2*x
z <- make_z(x,hx,hpx,-12.0, 12.0)
lower <- make_lower_bound(x, hx)
upper <- make_upper_bound(x, hx, hpx, z)
lower_old <- make_lower_bound_old(x, hx, -12.0, 12.0)
upper_old <- make_upper_bound_old(x, hx, hpx, z)
xx <- seq(-12, 12, by=0.2)

plot(xx, -xx^2, type='l')
lines(xx, lower(xx), col="red")
lines(xx, upper(xx), col="blue")


x <- seq(-10.5, 10.5, by=4)
hx <- -x^2
hpx <- -2*x
z <- make_z(x,hx,hpx,-Inf, Inf)
lower <- make_lower_bound(x, hx)
upper <- make_upper_bound(x, hx, hpx, z)
lower_old <- make_lower_bound_old(x, hx, -Inf, Inf)
upper_old <- make_upper_bound_old(x, hx, hpx, z)
xx <- seq(-20, 20, by=0.2)
plot(xx, -xx^2, type='l')
lines(xx, lower(xx), col="red")
lines(xx, upper(xx), col="blue")


#################################################################

# Benchmark lower
source("make_lowerupper_bound.R")
x <- seq(-10.5, 10.5, by=4)
hx <- -x^2
hpx <- -2*x
z <- make_z(x,hx,hpx,-Inf, Inf)
lower_old <- make_lower_bound_old(x, hx, -Inf, Inf)
lower <- make_lower_bound(x, hx)

xx <- runif(n = 1e5, min = -20, max = 20)

library(rbenchmark)
benchmark(lower_old(xx), lower(xx),
	  replications = 10, columns=c('test', 'elapsed', 'replications'))


# Benchmark upper
upper_old <- make_upper_bound_old(x, hx, hpx, z)
upper <- make_upper_bound(x, hx, hpx, z)

benchmark(upper_old(xx), upper(xx),
	  replications = 10, columns=c('test', 'elapsed', 'replications'))

