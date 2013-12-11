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
  #sample_interval <- sapply(sample, function(t) sum(intervalcdf<t)) + 1
  sample_interval <- rowSums(outer(sample, intervalcdf, function(x,y) x>y)) + 1
  t <- sample - c(0, intervalcdf)[sample_interval]
  
  sample_x <- inversecdf(t, sample_interval, factor1, hpx, z)
  return (sample_x)
}


inversecdf <- function(t, j, factor1, hpx, z) {
  return ( log(t/factor1[j]*hpx[j] + exp(hpx[j]*z[j])) / hpx[j] )
}
