initial<-seq(1,4,by = 0.03)
x_star<-seq(0.97,3.57, by = 0.3)


hpx<-function(hx, x, parameter = NA){ 
  # Problem here is that I used a stupid way to cover function that have less than or equal to
  # 2 parameters, and do you guys have any idea about how to revise this function?
  x = x
  if(sum(is.na(parameter)) == 1){
    h_prime<-diag(attr(numericDeriv(quote(hx(x)),"x"),"gradient"))
  }
  if(length(parameter) == 1){
    h_prime<-diag(attr(numericDeriv(quote(hx(x,parameter[1])),"x"),"gradient"))
  }
  if(length(parameter) == 2){
    h_prime<-diag(attr(numericDeriv(quote(hx(x,parameter[1],parameter[2])),"x"),"gradient"))
  }
  if(all(diff(h_prime) <= 0)) {return(h_prime)}
  else {warning("The log-concativity of this function is not good")}
}

# Here I used log-normal to test 
f<-function(x, mean = 0, sd = 1){dnorm(x, mean, sd, log = 1)}
hpx(f,initial,c(0,1))
# Here I used log-gamma to test
g<-function(x, shape, rate){dgamma(x, shape, rate, log = 1)}
hpx(g,initial,c(1,2))
hpx(g,initial,c(2,3))


make_z<-function(x, hx, hpx, left_bound = -Inf, right_bound = Inf){
  x = sort(x);k<-length(x)
  z_mid<-(hx[2:k]-hx[1:(k-1)]-x[2:k]*hpx[2:k]+x[1:(k-1)]*hpx[1:(k-1)]) / (hpx[1:(k-1)]-hpx[2:k])
  z<-c(left_bound, z_mid, right_bound)
  return(z)
}
hx<-dnorm(initial,log = 1)
hprime<-hpx(f,initial,c(0,1))
make_z(initial,hx,hprime)
z<-make_z(initial,hx,hprime)

make_upper_bound<-function(x, initial, z, hx, hpx = hprime, left_bound = -Inf, right_bound = Inf){
  x_star<-sort(x)
  x_star_interval<-sapply(x_star,function(x) return(sum(z<=x)))
  u_k<-vector("numeric", length(x_star))
  u_k <- hx[x_star_interval] + (x_star - initial[x_star_interval])*hprime[x_star_interval]
  #Possible issues in the value that may exceed max(initial)
  return(u_k)
}

make_upper_bound(x_star, initial, z, hx, hprime, left_bound = 0.001, right_bound = 0.5)

make_lower_bound<-function(x, initial, hx, hpx = hprime, left_bound = -Inf, right_bound = Inf){
  x_star<-sort(x)
  x_star_interval<-sapply(x_star,function(x) return(sum(initial<=x)))
  l_k<-vector("numeric", length(x_star))
  l_k<-((initial[x_star_interval+1] - x_star)*hx[x_star_interval] + (x_star - initial[x_star_interval])*hx[x_star_interval])/(initial[x_star_interval+1] - initial[x_star_interval])
  #Possible issues in the value that may exceed max(initial)
  return(l_k)
}

make_lower_bound(x_star,initial, hx, hprime)