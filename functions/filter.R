# -----------------------------------------------------------------------
# The function filter() used in adaptive_rejection_sampling_main.R
# It performs a squeezing test between the lower and upper bound, on the candidates in cand.

# Returns a list with
# $accepted: integer: position of last accepted candidate
# $update:   list: $cand is fist value not accepted
# 		   $w is corresponding draw form uniform dist
#                  Or if all candidates are accepted returns NA
# If all the cand are accepted, the update value is NA
# If no cand is accepted, the accepted value is NA

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

#----------------------------------------------------------------------------
# update_sample 
# Tests if the first candidate rejected by filer should be accepted after all
# Returns the index of the last accepted candidate.
update_sample <- function(cand, accepted, update, f, upper_bound) {
  flag <- (log(update$w) < log(f(update$cand))/upper_bound(update$cand))
  if (is.na(accepted)) {
    if (flag){
      return(1)
    }
  }
  else {
    if (flag) {
      return(accepted + 1)
    }
    else {
      return(accepted)
    }
  }
  return(0)
}

