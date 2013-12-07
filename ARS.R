## STAT 243 Final Project: adaptive rejection sampling ##
## Jinze Gu, Havard Kvamme, Reid Stevens and Yang Wang ##

## This is the primary function that performs the sampling from a log-concave probability density function using the adaptive rejection sampling algorithm developed by Gilks and Wild (Appl. Statist. 41, 337, 1992) ##

#f is the target density function, n is the number of required sample, test indicates whether to perform the test, k is the initial number of abscissaes

ars <- function(f, n, left_bound = -Inf, right_bound = Inf, x_init) {
  #Generating the initial abscissaes x
  x <- x_init
  hx <- log(f(x))
  hpx <- diag(attributes(numericDeriv(quote(log(f(x))), "x"))$gradient)
  sample <- rep(NA, n)
  count <- 0
  z <- make_z(x, hx, hpx, left_bound, right_bound) 

  while (count < n) {
    #Make the lower_bound, upper bound functions
    lower_bound <- make_lower_bound(x, hx, left_bound, right_bound)
    upper_bound <- make_upper_bound(x, hx, hpx, z)
    
    #Draw samples from the upper bound function
    cand <- sample_upper_bound(n - count, x, hx, hpx, z)
    
    #"update" records the first point in the sample that can not be accepted based on lowerbound
    #"accepted" records the candidates drawn from upper_bound function that are accepted by only using the lowerbound until the "update" point
    cand_filtered <- filter(cand, lower_bound, upper_bound, n-count)
    accepted <- cand_filtered$accepted
    update <- cand_filtered$update
    
    if (is.na(update[1])) {
      sample[(count+1):n] <- cand
      return(sample)
    }
    
    #Check the log-concaveness by checking if the log(f(update)) falls in between lower and upper bound function 
    if ((log(f(update)) < lower_bound(update)) | (log(f(update)) > upper_bound(update))) stop("The sample function is not log-concaved!")
    
    #Update the sample using cand_filtered
    update_sample(cand, accepted, update, count, f, upper_bound)
    
    #Update the abscissas x
    update_absci <- update_x(f, x, hx, hpx, update$cand)
    x <- update_absci$x
    hx <- update_absci$hx
    hpx <- update_absci$hpx
    z <- make_z(x, hx, hpx, left_bound, right_bound)
  }
  return (sample)
}
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

#make_lower_bound function returns the lower_bound function
#lower_bound function is set to be -Inf when argument falls outside the abscissaes range.
make_lower_bound_old <- function(x, hx, left_bound, right_bound) {
  x <- c(left_bound, x, right_bound)
  hx <- c(-Inf, hx, -Inf)
  lower_bound <- function(xx) {
    k=length(x)
    index <- rowSums(outer(xx, x, function(x1,x2) x1>x2))
    index[index==0] <- 1
    index[index==k] <- k - 1
    lx <- ((x[index+1]-xx)*hx[index] + (xx-x[index])*hx[index+1])/(x[index+1]-x[index])
    lx[index %in% c(1,k-1)] <- -Inf
    return(lx)
  }
  return(lower_bound)
}

#make_upper_bound function returns the upper_bound function
make_upper_bound_old <- function(x, hx, hpx, z) {
  upper_bound <- function(xx) {
    index <- rowSums(outer(xx, z, function(x1,x2) x1>x2))
    index[index==0] <- 1
    index[index==length(z)] <- length(z) - 1
    ux <- hx[index] + (xx-x[index])*hpx[index]
    return (ux)
  }
  return(upper_bound)
}

#################################################################
make_lower_bound <- function(x, hx, left_bound, right_bound) {
  x <- c(left_bound, x, right_bound)
  hx <- c(-Inf, hx, -Inf)
  lower_bound <- approxfun(x,hx)
  return(lower_bound)
}

make_upper_bound <- function(x, hx, hpx, z) {
  upper_bound_eval_z <- function(xx) {
    index <- rowSums(outer(xx, z, function(x1,x2) x1>x2))
    index[index==0] <- 1
    index[index==length(z)] <- length(z) - 1
    ux <- hx[index] + (xx-x[index])*hpx[index]
    return (ux)
  }
  y <- upper_bound_eval_z(z)
  upper_bound <- approxfun(z,y)
  return(upper_bound)
}

#------------------------------------------------
#Calculating the z based on abscissae x. z[1]=left_bound, z[k+1]=right_bound
make_z <- function(x, fx, fpx, left_bound, right_bound) {
  k <- length(x)
  xm1 <- c(0, x)
  fxm1 <- c(0, fx)
  fpxm1 <- c(0, fpx)
  x <- c(x, 0)
  fx <- c(fx, 0)
  fpx <- c(fpx, 0)
  
  z <- (fx - fxm1 - x*fpx + xm1*fpxm1) / (fpxm1 - fpx)
  z[1] <- left_bound
  z[k+1] <- right_bound
  
  return (z)
}

#sample_upper_bound function samples m points from the upper bound funtion using inverse cdf method. The inverse cdf is calculated analytically and implemented in the auxillary function inversecdf.
sample_upper_bound <- function(m, x, hx, hpx, z) {
  k <- length(x)
  
  #Calcualte the normalized integration of upper_bound funtion at each linear interval
  factor1 <- exp(hx - x*hpx) 
  factor2 <- (exp(c(0, hpx) * z)[2:(k+1)] - exp(hpx*z[1:k])) / hpx
  I <- factor1 * factor2
  Inormalize <- sum(I)
  
  intervalcdf <- cumsum(I)
  sample <- runif(m, min=0, max=Inormalize)
  #calculate which interval the sample falls into
  sample_interval <- sapply(sample, function(t) sum(intervalcdf<t)) + 1
  t <- sample - c(0, intervalcdf)[sample_interval]
  
  sample_x <- inversecdf(t, sample_interval, factor1, hpx, z)
  return (sample_x)
}


inversecdf <- function(t, j, factor1, hpx, z) {
  return ( log(t/factor1[j]*hpx[j] + exp(hpx[j]*z[j])) / hpx[j] )
}#update_x function adds the update point into the abscissaes vector, and also checks if the derivative at the update point is between its neighbours(i.e., checking the cancaveness)

update_x <- function(f, x, hx, hpx, update) {
  if (is.na(update)){
    return (list(x=x, hx=hx, hpx=hpx))
  }
  else {
    hx_update <- log(f(update))
    xx <- update
    hpx_update <- attributes(numericDeriv(quote(log(f(xx))), "xx"))$gradient[1, 1]
    rank_x <- rank(c(x, update))
    new_x <- rep(NA, length(x)+1)
    new_hx <- rep(NA, length(x)+1)
    new_hpx <- rep(NA, length(x)+1)
    
    new_x[rank_x] <- c(x, update)
    new_hx[rank_x] <- c(hx, hx_update)
    new_hpx[rank_x] <- c(hpx, hpx_update)
    
    #Check concaveness by testing if h'(x) is decreasing
    if (is.unsorted(-new_hpx)) stop("The sample function is not log-concaved!")
    
    return(list(x=new_x, hx=new_hx, hpx=new_hpx))
  }
}

##Test case for update_x###
#f <- function(x) return(x^2)
#x <- seq(1.5, 10.5, by=1)
#hx <- log(f(x))
#hpx <- 2/x

#update_x(f, x, hx, hpx, 5.1)
#update_x(f, x, hx, hpx, 12)
#update_x(f, x, hx, hpx, 0.5)

#ff <- function(x) return(exp(x^2))
#update_x(ff, x, hx, hpx, 5.1)
