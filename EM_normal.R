#Data Simulation
n=1e4
mu=2
v=1
y=rnorm(n,mean = mu, sd = sqrt(v))
d=y<2.2 
m <- sum(d)
yCen=y 
yCen[!d]=2.2 
lambda=rep(NA,50)
lambda[1]=0
completeData=y
#suppose the censored part sigma is known,which is 1
for(i in 2:length(lambda)){
  # E-step
  t <- 2.2-lambda[i-1]
  completeData[!d] <- lambda[i-1]+dnorm(t,0,1)/(1-pnorm(t,0,1))
  # M-step
  lambda[i] <- (sum(completeData[d])+sum(completeData[!d]))/n
}

tail(lambda)

#suppose sigma is unknown, sigma needs to be initialed as var(y!=2.2)
T_1 <- sum(completeData[d])
T_2 <- sum(completeData[d]^2)
lambda=rep(NA,100)
lambda[1]=0
sigma=rep(NA,100)
sigma[1]=var(y[d])
for(i in 2:length(lambda)){
  # E-step
  t <- (2.2-lambda[i-1])/sqrt(sigma[i-1])
  S_1 <- (lambda[i-1]+sqrt(sigma[i-1])*dnorm(t,0,1)/(1-pnorm(t,0,1)))*(n-m)
  S_2 <- (n-m)*(lambda[i-1]^2 + sigma[i-1]) + (n-m)*(lambda[i-1]+2.2)*(sqrt(sigma[i-1])*dnorm(t,0,1))/(1-pnorm(t,0,1))
  # M-step
  lambda[i] <- (T_1 + S_1)/n
  sigma[i] <- (T_2+S_2)/n - ((T_1+S_1)^2)/(n^2)
}
tail(lambda)
tail(sigma)

#the real mean and variance of data is 2 and 1, when censored part sigma is known, this can
#perfectly recover the original estimator. When sigma is unknown, the algorithm works not as
#good as before, the result has less accuracy.
