#Test case for make_z
source("make_z.R")

x <- seq(-10, 10, by = 2)
fx <- x^2
fpx <- 2*x
plot(x, fx)
points(make_z(x, fx, fpx, -12, 12), make_z(x, fx, fpx, -12, 12)^2, col="red")
