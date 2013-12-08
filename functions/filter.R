# the function filter() used in adaptive_rejection_sampling_main.R
# It performs a squeezing test between the lower and upper bound.

# Returns a list with
# $accepted: integer: position of last accepted candidate
# $update:   list: $cand is fist value not accepted
# 		   $w is corresponding draw form uniform dist
#                  Or if all candidates are accepted returns NA

filter <- function(cand, lower_bound, upper_bound, m){
  w <- runif(m) 
  squeeze <- lower_bound(cand) - upper_bound(cand)

  # Gives NA if all are accepted
  i <- which(log(w) > squeeze)[1]

  if (is.na(i)) {
    return(list(accepted=m, update=NA))
  } 
  else{
    if (i==1) {
      return(list(accepted=NA, update=list(cand=cand[i], w=w[i])))
    }
    else {
      return(list(accepted=(i-1), update=list(cand=cand[i], w=w[i])))
    }
  }
}
##########################

# update_sample()
# Check if first point not excepted in filter is accepted, and updates sample
# output NULL. Assign to sample and count in parent scope.
# Have to rewrite 
update_sample <- function(cand, accepted, update, count, f, upper_bound) {
  flag <- (log(update$w) < log(f(update$cand))/upper_bound(update$cand))
  if (is.na(accepted)) {
    if (flag){
      sample[(count+1)] <<- update$cand
      count <<- count + 1
    }
  }
  else {
    k <- accepted
    if (flag) {
      sample[(count+1):(count+k)] <<- cand[1:accepted]
      sample[count+k+1] <<- update$cand
      count <<- count + k + 1
    }
    else {
      sample[(count+1):(count+k)] <<- cand[1:accepted]
      count <<- count + k 
    }
  }
  invisible(NULL)
}

