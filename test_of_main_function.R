#Test for the main function ars()

source("ars_all_functions.R")
#should return a straight line
set.seed(0)
system.time(ars_sample <- ars(dnorm, 1e4, x_init=c(-1,1)))
qqnorm(ars_sample)
qqline(ars_sample)

#My code becomes really slow for n=1e5
ars_sample <- ars(dnorm, 1e5, x_init=c(-1,1))
qqnorm(ars_sample)
qqline(ars_sample)

#Test for extra arguments
ars_sample <- ars(dnorm, 1e4, x_init=c(1,3), mean=2, sd=5)
qqnorm(ars_sample)
qqline(ars_sample)

#Test for bounded distribution
ars_sample <- ars(dnorm, 1e4, left_bound=-2, right_bound=2,x_init=c(-1,1))
hist(ars_sample)

#To check if it can catch bad initial points
ars_sample <- ars(dnorm, 1e4, x_init=c(1,2))
ars_sample <- ars(dnorm, 1e4, left_bound=0.1, x_init=c(1,2))
hist(ars_sample)

#To check if it can catch non-log-concaveness
testFun <- function(x) exp(x^2)
ars_sample <- ars(testFun, 1e4, -2, 2, x_init=c(-1,1))


#----------------------------------------------
# New
source("ARS.R")
set.seed(0)
system.time(ars_sample2 <- ars(dnorm, 1e4, x_init=c(-1,1)))
par(mfrow = c(2,1))
qqnorm(ars_sample2)
qqline(ars_sample2)

#My code becomes really slow for n=1e5
ars_sample2 <- ars(dnorm, 1e5, x_init=c(-1,1))
qqnorm(ars_sample2)
qqline(ars_sample2)

#Test for extra arguments
ars_sample2 <- ars(dnorm, 1e4, x_init=c(1,3), mean=2, sd=5)
qqnorm(ars_sample2)
qqline(ars_sample2)

#Test for bounded distribution
ars_sample2 <- ars(dnorm, 1e4, left_bound=-2, right_bound=2,x_init=c(-1,1))
hist(ars_sample2)

#To check if it can catch bad initial points
ars_sample2 <- ars(dnorm, 1e4, x_init=c(1,2))
ars_sample2 <- ars(dnorm, 1e4, left_bound=0.1, x_init=c(1,2))
hist(ars_sample2)

#To check if it can catch non-log-concaveness
testFun <- function(x) exp(x^2)
ars_sample2 <- ars(testFun, 1e4, -2, 2, x_init=c(-1,1))
