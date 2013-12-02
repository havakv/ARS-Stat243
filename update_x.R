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
    
    return(list(x=new_x, hx=new_hx, hpx=new_hpx))
  }
}

##Test case for update_x###
f <- function(x) return(x^2)
x <- seq(-10.5, 10.5, by=1)
hx <- log(f(x))
hpx <- 2/x

update_x(f, x, hx, hpx, 5.1)
update_x(f, x, hx, hpx, 12)
update_x(f, x, hx, hpx, -100)


