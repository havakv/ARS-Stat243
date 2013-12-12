#############################################################################
# ARS function test
set.seed(0)
hist(ars(dnorm, 10000, x_init = c(-1,1)), main = "ARS sampling from a truncated normal distribution", freq = 0)

hist(ars(dnorm, 10000, x_init = c(-1,1)), main = "ARS sampling from a truncated normal distribution", freq = 0)

cat("The following does not satisfy our input requirements. \n")
hist(ars(dnorm, 2000, x_init = runif(2)), main = "ARS sampling from a normal distribution")

# qqplot
ars_sample <- ars(dnorm, 10000, x_init=c(-1,1))
qqnorm(ars_sample)
qqline(ars_sample)

#Test for extra arguments
ars_sample <- ars(dnorm, 10000, x_init=c(1,3), mean = 2, sd = 5)
qqnorm(ars_sample)
qqline(ars_sample)

#Test for bounded distribution
ars_sample <- ars(dnorm, 10000, left_bound=-2, right_bound=2,x_init=c(-1,1))
hist(ars_sample)

#To check if it can catch bad initial points
ars_sample <- ars(dnorm, 10000, x_init=c(1,2))
ars_sample <- ars(dnorm, 10000, left_bound=0.1, x_init=c(1,2))
hist(ars_sample)

#############################################################################
# Make_z function test

set.seed(0)
x <- sort(rnorm(100))
fx <- dnorm(x, log = 1)
fpx <- diag(attributes(numericDeriv(quote(dnorm(x, log = 1)), "x"))$gradient)
plot(x, fx, type = "p", pch = 16)
points(make_z(x, fx, fpx, -12, 12), dnorm(make_z(x, fx, fpx, -12, 12),log = 1), pch = 17, col = "red")
legend("topright", c("log-normal plot", "make_z generation plot"), pch = c(16, 17))


cat("Since z is just the point from the density function(in this case, lognormal), so we can just plot two set of points together to see if they fall into the same distribution trend. \n")

#############################################################################
# make_upper/lower bound tests
# Test 1
set.seed(0)
x <- runif(100)
hx <- dnorm(x, log = 1)
hpx <- diag(attributes(numericDeriv(quote(dnorm(x, log = 1)), "x"))$gradient)
z <- make_z(x,hx,hpx,-12.0, 12.0)
lower <- make_lower_bound(x, hx)
upper <- make_upper_bound(x, hx, hpx, z)
xx <- seq(-5, 5, by=0.1)
cat("This shows that all lower bound should less than upper bound by our generation. \n")
all((lower(xx) <= upper(xx)) == TRUE)

#Test 2
cat("Now we use -x^2 as a function to plot upper and lower bound. \n")
x <- seq(-10, 10, by=4)
cat("Notice here we need to define a function and input the name of this function \n")
f<-function(x){-x^2}
fpx <- diag(attributes(numericDeriv(quote(f(x)), "x"))$gradient)
z <- make_z(x,f(x),fpx,-Inf, Inf)
lower <- make_lower_bound(x, f(x))
upper <- make_upper_bound(x, f(x), fpx, z)
xx <- seq(-10, 10, by=0.2)
plot(xx, -xx^2, type='l')
lines(xx, lower(xx), lty = 2)
lines(xx, upper(xx), lty = 6)
legend("topright", c("upper_bound", "lower_bound"), lty = c(6,2))


cat("In test 1, we basically just check if the corresponding lower bound is less than upper bound, and the function passes the test as long as the output is TRUE. Then we tried to plot our upper/lower bound using previous function we defined: $make_z$.We use a logconcave function $f(x)=x^2$, and it turns out the plot is what we had expected. \n")

#############################################################################
# sample_upper_bound test

set.seed(0)
x <- sort(rnorm(10000))
cat("Notice here we need to define a function and input the name of this function \n")
f<-function(x){-x^2}
fpx <- diag(attributes(numericDeriv(quote(f(x)), "x"))$gradient)
z <- make_z(x,f(x),fpx,-Inf, Inf)
upper <- make_upper_bound(x, f(x), fpx, z)
xx <- seq(-10, 10, by=0.2)
hist(sample_upper_bound(2000, x, f(x), fpx, z), freq = 0)

cat("Since the histogram of sample_upper_bound function matchs the trend of $f(x) = -x^2$, so it is obvious that our sample_upper_bound function works as we expected. \n")

#############################################################################
# filter function test
set.seed(0)
f <- dnorm
n <- 10
x <- sort(rnorm(7))
left_bound   <- -10
right_bound  <- 10
hx <- log(f(x))
hpx <- diag(attributes(numericDeriv(quote(log(f(x))), "x"))$gradient)
sample <- rep(NA, n)
count <- 0
z <- make_z(x, hx, hpx, left_bound, right_bound) 
lower_bound <- make_lower_bound(x, hx)
upper_bound <- make_upper_bound(x, hx, hpx, z)
cand <- sample_upper_bound(n - count, x, hx, hpx, z)
cand_filtered <- filter(cand, lower_bound, upper_bound, n - count)

cat("Basically, we have our candidates and cand_filtered, and the function does what it is suppose to do. \n")
cand
cand_filtered

#############################################################################
# update_sample test

set.seed(0)
f<-dnorm
x <- sort(rnorm(7))
hx <- log(f(x))
hpx <- diag(attributes(numericDeriv(quote(log(f(x))), "x"))$gradient)
sample <- c(1,2,3, rep(NA, 10))
cat("When accepted point is not 'NA' \n")
accepted <- cand_filtered$accepted
update <- cand_filtered$update
count <- length(na.omit(sample))
upper_bound <- make_upper_bound(x, hx, hpx, z)
update_sample(cand = cand, accepted=accepted, update = update, f = f, upper_bound = upper_bound)
cat("When accepted point is 'NA' \n")
accepted <- NA
update_sample(cand, accepted, update,  f, upper_bound)

cat("In this case, the update_sample_function returns the index of the last accepted point, and it returns 1 if the point satisfies the test "flag" and we accepted no point. It returns 0 if the point does not satisfy "flag" and we accepted no point. According to the test, it works correctly with specific output. \n")

#############################################################################
# update_x test

set.seed(0)
x<-seq(0.0001,1, by = 0.2)
f<-dnorm
hx <- log(f(x))
hpx <- diag(attributes(numericDeriv(quote(log(f(x))), "x"))$gradient)
update <- cand_filtered$update
upper_bound <- make_upper_bound(x, hx, hpx, z)
update_x(f = f, x = x, hx = hx, hpx = hpx, update = update$cand)

cat("The update_x outputs the updated points, their evaluations under function h(x) and their corresponding derivatives \n")


#############################################################################
