## STAT 243 Final Project: adaptive rejection sampling ##
## Jinze Gu, Haavard Kvamme, Reid Stevens and Yang Wang ##

## This is the primary function that performs the sampling from a log-concave probability density function using the adaptive rejection sampling algorithm developed by Gilks and Wild (Appl. Statist. 41, 337, 1992) ##

#f is the target density function, n is the number of required sample, test indicates whether to perform the test, k is the initial number of abscissas

ars <- function(f, n, left_bound = -Inf, right_bound = Inf, x_init, ...) {
  #Generating the initial abscissaes x
  x <- x_init
  x <- sort(x)
  hx <- log(f(x, ...))
  hpx <- diag(attributes(numericDeriv(quote(log(f(x, ...))), "x"))$gradient)
  if (((hpx[1] < 0) & (left_bound==-Inf)) | ((hpx[length(hpx)]>0) & (right_bound==Inf))) stop("The derivatie at the first/last initial point must be positive/negative")
  sample <- rep(NA, n)
  count <- 0
  z <- make_z(x, hx, hpx, left_bound, right_bound) 

  while (count < n) {
    #Make the lower_bound, upper bound functions
    lower_bound <- make_lower_bound(x, hx)
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
    if ((log(f(update$cand, ...)) < lower_bound(update$cand)) | (log(f(update$cand, ...)) > upper_bound(update$cand))) stop("The sample function is not log-concaved!")
    
    #Update the sample using cand_filtered
    index <- update_sample(cand, accepted, update, f, upper_bound)
    if (index != 0) {
      sample[(count+1):(count+index)] <- cand[1:index]
      count <- count + index
    }
    #Update the abscissas x
    update_absci <- update_x(f, x, hx, hpx, update$cand, ...)
    x <- update_absci$x
    hx <- update_absci$hx
    hpx <- update_absci$hpx
    z <- make_z(x, hx, hpx, left_bound, right_bound)
  }
  return (sample)
}
