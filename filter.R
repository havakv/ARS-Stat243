#!/usr/bin/env Rscript
# the function filter() used in adaptive_rejection_sampling_main.R
# It performs a squeezing test between the lower and upper bound.

# Should we make lower and upper bond functions or what?
# What looks nice, and is most effective?

# Returns a list with
# $accepted: vector of accepted candidates
# $update:   first element to not be accepted

filter <- function(cand, lower_bound, upper_bound){
  w <- runif(n = length(cand)) # Should use value from script...
  i <- 0
  accepted <- TRUE
  while (accepted && i < length(cand)){ # Should use value from script...
    i <- 1
    if (w[i]>exp(lower_bound[i] - upper.bound[i])) {
    	accepted <- FALSE
    }
  }
  if (i == 1) {
    return(list(accepted=NULL, update=cand(1)))
  } 
  else{
    return(list(accepted=cand[1:(i-1)], update=cand(i)))
  }
  # consider returning index for update as well...
}


# -----------------------
# FASTER ???
# -----------------------

# Maybe it is faster to evaluate the bounds as functions. The we do not need to evaluate more than up to update.
filter1 <- function(candd){
  w <- runif(n = length(cand)) # Should use value from script...
  i <- 0
  accepted <- TRUE
  while (accepted && i < length(cand)){ # Should use value from script...
    i <- 1
    # Now the bounds are functions to evaluate one point.
    # Can we use make_lower_bound() directly? In that case, change name.
    if (w[i]>exp(lower_bound(i) - upper.bound(i))) {
    	accepted <- FALSE
    }
  }
  if (i == 1) {
    return(list(accepted=NULL, update=cand(1)))
  } 
  else{
    return(list(accepted=cand[1:(i-1)], update=cand(i)))
  }
  # consider returning index for update as well...
}

# ---------------------------------------------------------

# Using vectorized calculations 
# Have to evaluate all, but probably faster. 
# Do not think which should take much work, but we should probably test everything
filter2 <- function(cand, lower_bound, upper_bound){
  w <- runif(n = length(cand)) # Should use value from script...
  squeeze <- exp(lower_bound - upper.bound)

  # Gives NA if the all is accepted
  i <- which(w > squeeze)[1]

  if (is.na(i)) {
    return(list(accepted=NULL, update=cand(1)))
  } 
  else {
    return(list(accepted=cand[1:(i-1)], update=cand(i)))
  }
  # consider returning index for update as well...
}
