# Test of filter() and update_sample

# Test functions
set.seed(0)
source("../ars_all.R")
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

# filter()
set.seed(0)
cand_filtered <- filter(cand, lower_bound, upper_bound, n-count)
cand_filtered

accepted <- cand_filtered$accepted
update <- cand_filtered$update

# update_sample()
update_sample(cand, accepted, update, count, f, upper_bound)
sample
count

