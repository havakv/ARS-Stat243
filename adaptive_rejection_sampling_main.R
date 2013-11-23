## STAT 243 Final Project: adaptive rejection sampling ##
## Jinze Gu, Havard Kvamme, Reid Stevens and Yang Wang ##

## This is the primary function that performs the sampling from a log-concave probability density function using the adaptive rejection sampling algorithm developed by Gilks and Wild (Appl. Statist. 41, 337, 1992) ##

#f is the target density function, n is the number of required sample, test indicates whether to perform the test, k is the initial number of abscissaes

ars <- function(f, n, test=TRUE, left_bound = -Inf, right_bound = Inf, k=3) {
  #Generating the initial abscissaes x
  x <- initial(left_bound, right_bound, k)
  hx <- log(f(x))
  hpx <- log(fp(f, x))
  sample <- rep(NA, n)
  count <- 0
  z <- make_z(x, hx, hpx, left_bound, right_bound) 

  while (count < n) {
    #Make the lower_bound, upper bound functions
    lower_bound <- make_lower_bound(x, hx, z, left_bound, right_bound)
    upper_bound <- make_upper_bound(x, hx, hpx, z, left_bound, right_bound)
    
    #Draw samples from the upper bound function
    cand <- sample_upper_bound(n - count, x, hx, hpx, z, left_bound, right_bound)
    
    #"update" records the first point in the sample that can not be accepted based on lowerbound
    #"accepted" records the candidates drawn from upper_bound function that are accepted by only using the lowerbound until the "update" point
    cand_filtered <- filter(cand, lower_bound, upper_bound)
    accepted <- cand_filtered[1]
    update <- cand_filtered[2]
    
    #Update the sample using cand_filtered
    sample <- update_sample(sample, cand_filtered)
    count <- length(na.omit(sample))
    
    #Update the abscissaes x
    update_absci <- update_x(x, hx, hpx, update)
    x <- update_absci[1]
    hx <- update_absci[2]
    hpx <- update_absci[3]
    z <- make_z(x, hx, hpx, left_bound, right_bound)
  }
  return (sample)
}