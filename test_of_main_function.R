#Test for the main function ars()

#should return a straight line
ars_sample <- ars(dnorm, 1e4, x_init=c(-1,1))
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

