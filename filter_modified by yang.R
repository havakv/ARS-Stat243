#Function filer takes in the cand generated from sample_upper_bound and return the accepted sample and update point. 
#If all the cand are accepted, the update value is NA
#If no cand is accepted, the accepted value is NA

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
      return(list(accepted=NA, update=cand[i]))
    }
    else {
      return(list(accepted=cand[1:(i-1)], update=cand[i]))
    }
  }
}


#update_sample function adds the accepted sample into the final sample, and it also check whether to add update into the final sample
update_sample <- function(sample, accepted, update, count, f, upper_bound) {
  w <- runif(1)
  flag <- (log(w) < log(f(update))/upper_bound(update))
  if (is.na(accepted[1])) {
    if (flag) sample[(count+1)] <- update
  }
  else {
    k <- length(accepted)
    if (flag) {
      sample[(count+1):(count+1+k)] <- c(accepted, update)
    }
    else {
      sample[(count+1):(count+k)] <- c(accepted)
    }
  }
  return(sample)
}

#Test case for update_sample
sample <- c(1,2,3, rep(NA, 10))
accepted <- c(4,5,6)
update <- 7
count <- length(na.omit(sample))
update_sample(sample, accepted, update, count)
accepted <- NA
update_sample(sample, accepted, update, count)