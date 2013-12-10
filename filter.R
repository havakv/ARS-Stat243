#!/usr/bin/env Rscript
# the function filter() used in adaptive_rejection_sampling_main.R
# It performs a squeezing test between the lower and upper bound.

# Returns a list with
# $accepted: vector of accepted candidates
# $update:   first element to not be accepted


filter <- function(cand, lower_bound, upper_bound){
  w <- runif(n = length(cand)) # Should use value from script
  squeeze <- lower_bound(cand) - upper_bound(cand)

  # Gives NA if all are accepted
  i <- which(log(w) > squeeze)[1]

  if (is.na(i)) {
    return(list(accepted=cand, update=NA))
  } 
  else{
    if (i==1) {
      return(list(accepted=NA, update=list(x=cand[i], w=w[i])))
    }
    else {
      return(list(accepted=cand[1:(i-1)], update=list(x=cand[i], w=w[i])))
    }
  }
}

#update_sample function adds the accepted sample into the final sample, and it also check whether to add update into the final sample
update_sample <- function(accepted, update, count, f, upper_bound) {
  # Update should contain w
  flag <- (log(update$w) < log(f(update$x))/upper_bound(update$x))
  if (is.na(accepted[1])) {
    if (flag){
      sample[(count+1)] <<- update$x
      count <<- count + 1
    }
  }
  else {
    k <- length(accepted)
    if (flag) {
      sample[(count+1):(count+k)] <<- accepted
      sample[count+k+1] <<- update$x
      count <<- count + k + 1
    }
    else {
      sample[(count+1):(count+k)] <<- accepted
      count <<- count + k
    }
  }
  return(NULL)
}

#---------------------------------------------------
# FASTER!!!
# Both methods require some differences in the main funcrtion

filter2 <- function(cand, lower_bound, upper_bound, m){
  #m <- length(cand) # Should use value from script
  #m should be (n - cand)
  w <- runif(m) 
  squeeze <- lower_bound(cand) - upper_bound(cand)

  # Gives NA if all are accepted
  i <- which(log(w) > squeeze)[1]

  if (is.na(i)) {
    return(list(accepted=m, update=NA))
  } 
  else{
    if (i==1) {
      return(list(accepted=NA, update=list(x=cand[i], w=w[i])))
    }
    else {
      return(list(accepted=(i-1), update=list(x=cand[i], w=w[i])))
    }
  }
}
# Have to rewrite 
update_sample2 <- function(cand, accepted, update, count, f, upper_bound) {
  update$w <- runif(1) # REMOVE THIS!!!!!!!!!!!!!!!
  # Update should contain w !!!!!!!!!!!!!!!!
  flag <- (log(update$w) < log(f(update$x))/upper_bound(update$x))
  if (is.na(accepted)) {
    if (flag){
      sample[(count+1)] <<- update$x
      count <<- count + 1
    }
  }
  else {
    k <- accepted
    if (flag) {
      sample[(count+1):(count+k)] <<- cand[1:accepted]
      sample[count+k+1] <<- update$x
      count <<- count + k + 1
    }
    else {
      sample[(count+1):(count+k)] <<- cand[1:accepted]
      count <<- count + k 
    }
  }
  return(NULL)
}

# ------------------------------------------
# Test functions
set.seed(0)
source("ars_all functions.R")
# Now run filter2 and update_sample2!!!!!!!!!!!!!!!!!!!!!
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
lower_bound <- make_lower_bound(x, hx, left_bound, right_bound)
upper_bound <- make_upper_bound(x, hx, hpx, z)
cand <- sample_upper_bound(n - count, x, hx, hpx, z)

set.seed(0)
cand_filtered <- filter(cand, lower_bound, upper_bound)
set.seed(0)
cand_filtered2 <- filter2(cand, lower_bound, upper_bound, n-count)
cand_filtered
cand_filtered2

accepted <- cand_filtered$accepted
update <- cand_filtered$update
set.seed(0)
sample1 <- update_sample(sample, accepted, update, count, f, upper_bound)
accepted <- cand_filtered2$accepted
update <- cand_filtered2$update
set.seed(0)
update_sample2(cand, accepted, update, count, f, upper_bound)
sample1
sample


############################################
x <- rep(NA, 1e7)
x[2:1001] <- rnorm(1000)
noLength <- function(x) {length(x)}
Length <- function(x) {length(na.omit(x))}
ifLength <- function(x) {
  if (is.na(x[1])) {
    x = x[-1]
    length(x)
  }
  else length(x)
}
library(rbenchmark)

benchmark({noLength(x)}, {Length(x)}, {ifLength(x)},
	  replications = 10, columns=c('test', 'elapsed', 'replications'))

benchmark({noLength(x)}, 
	  replications = 1000, columns=c('test', 'elapsed', 'replications'))

# ------------------------------------------
byVal <- function(x) { 
  x[1] <- 5
  return(x)
}
byRef <- function() { 
  x[1] <<- 5
  return(NULL)
}
benchmark({x <- byVal(x)}, byRef(),
	  replications = 10, columns=c('test', 'elapsed', 'replications'))
