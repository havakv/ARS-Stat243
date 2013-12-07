#!/usr/bin/env Rscript
# the function filter() used in adaptive_rejection_sampling_main.R
# It performs a squeezing test between the lower and upper bound.

# Returns a list with
# $accepted: vector of accepted candidates
# $update:   first element to not be accepted


filter <- function(cand, lower_bound, upper_bound){
  w <- runif(n = length(cand)) # Should use value from script
  squeeze <- lower_bound(cand) - upper_bound(cand)

  # Gives NA if the all is accepted
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
  if (log(update$w) < log(f(update$x))/upper_bound(update)) {
    accepted <- c(accepted, update)
  }
  accepted <- na.exclude(accepted)
  k <- length(accepted)
  if (k != 0) {
    sample[(count+1):(count+k)] <<- accepted
  }
  count <<- count + k
  return(NULL)
}

# -------------------------------------------
x <- rep(NA, 1e7)
x[1:1000] <- rnorm(1000)
noLength <- function(x) {length(x)}
Length <- function(x) {length(na.omit(x))}
library(rbenchmark)

benchmark({noLength(x)},
	  {Length(x)},
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
