setwd("C:/Users/User/Documents/Studium_MA/3. Semester/Statistical Consulting/R")
seal_clean <- read.csv("seal_data_cleaned.csv")
seal1 <- seal_clean[which(seal_clean$sealID=="1"),]


# functions
# for (i in 1:21) {
#   assign(paste0("seal",i),seal_clean[which(seal_clean$sealID==i), ])
# }

mle <- function(obs,mu0,sigma0,gamma0,N){
  parvect <- pn2pw(mu0,sigma0,gamma0,N)
  mod <- nlm(mllk,parvect,obs,N,print.level=2,iterlim=1000,stepmax=5)
  pn <- pw2pn(mod$estimate,N)
  return(list(mu=pn$mu,sigma=pn$sigma,gamma=pn$gamma,delta=pn$delta,mllk=mod$minimum))
}

N <- 2
parvect <- c( c(50, 150), c(20, 40), matrix(rep(0.5, 4), nrow = 2) )
#for(i in 1:21) {
#  obs <- seal_clean[which(seal_clean$sealID==i), ]
#  obs <- obs[1:100, 17]
#  assign(paste0("lk",i), mllk(parvect, obs, N))
#  }

mllk <- function(parvect,obs,N){
  vec <- rep(NA, 21)
  for(i in 1:21) {
    obs <- seal_clean[which(seal_clean$sealID==i), ]
    obs <- obs[1:100, 17]
  lpn <- pw2pn(parvect,N)
  gamma<-lpn$gamma
  n <- length(obs)
  allprobs <- matrix(rep(1,N*n),nrow=n)
  ind.step <- which(!is.na(obs))
  for (j in 1:N){
    step.prob <- rep(1,n)
    step.prob[ind.step] <- dgamma(obs[ind.step],shape=lpn$mu[j]^2/lpn$sigma[j]^2,scale=lpn$sigma[j]^2/lpn$mu[j])
    allprobs[,j] <- step.prob
  }
  foo <- lpn$delta  
  lscale <- 0
  lscale <- forwardalgo(foo, gamma, allprobs, lscale, n)
  vec[i] <- lscale
  }
  lk <- prod(vec)
  return(-lk)
}

mod <- nlm(mllk, parvect=parvect,obs=seal_clean, N=2,print.level=2,iterlim=1000,stepmax=5)
