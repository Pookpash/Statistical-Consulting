setwd("C:/Users/User/Documents/Studium_MA/3. Semester/Statistical Consulting/R")
seal_clean <- read.csv("seal_data_cleaned.csv")
seal1 <- seal_clean[which(seal_clean$sealID=="1"),]

obs <- cbind(seal1$surf.dur[1:100], seal1$dive.dur[1:100], seal1$max.dep[1:100], seal1$steplen[1:100])
colnames(obs) <- c("surf.dur", "dive.dur", "max.dep", "steplen")

## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
pn2pw <- function(mu1,mu2,mu3,mu4,sigma1,sigma2,sigma3,sigma4,gamma,N){
  for(i in 1:4){
    mu <- cbind(mu1, mu2, mu3, mu4)
    assign(paste0("tmu", i), log(mu[,i]))    
  }
  for(i in 1:4){
    sigma <- cbind(sigma1, sigma2, sigma3, sigma4)
    assign(paste0("tsigma", i), log(sigma[,i]))    
  }  
  tgamma <- NULL
  if(N>1){
    foo <- log(gamma/diag(gamma))           ### gamma ist eine Matrix!
    tgamma <- as.vector(foo[!diag(N)])
  }
  parvect <- c(tmu1,tmu2,tmu3,tmu4,tsigma1,tsigma2,tsigma3,tsigma4,tgamma)
  return(parvect)
}
## function that performs the inverse transformation
pw2pn <- function(parvect,N){
  mu <- exp(parvect[1:(4*N)])
  for(i in 1:4) {
    assign(paste0("mu", i), mu[(2*i-1):(2*i)])    
  }   
  sigma <- exp(parvect[(4*N+1):(4*2*N)])
  for(i in 1:4) {
    assign(paste0("sigma", i), sigma[(2*i-1):(2*i)])    
  }   
  gamma <- diag(N)
  if(N>1){
    gamma[!gamma] <- exp(parvect[(4*2*N+1):(length(parvect))])
    gamma <- gamma/apply(gamma,1,sum)                 ### gamma = Matrix
  }
  delta <- solve(t(diag(N)-gamma+1),rep(1,N))
  return(list(mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,sigma4=sigma4,gamma=gamma,delta=delta))
}

mllk <- function(parvect,obs,N){
  lpn <- pw2pn(parvect,N)
  gamma <- lpn$gamma
  n <- length(obs[,1])
  allprobs <- matrix(rep(1,N*n),nrow=n)
  ind.surf <- which(!is.na(obs[,1]))
  ind.dive <- which(!is.na(obs[,2]))
  ind.dep <- which(!is.na(obs[,3]))
  ind.step <- which(!is.na(obs[,4]))
    for (j in 1:N){
    surf.prob <- dive.prob <- dep.prob <- step.prob <- rep(1,n)
    surf.prob[ind.surf] <- dgamma(obs[ind.surf, 1],shape=lpn$mu1[j]^2/lpn$sigma1[j]^2,scale=lpn$sigma1[j]^2/lpn$mu1[j])
    dive.prob[ind.dive] <- dgamma(obs[ind.dive, 2],shape=lpn$mu2[j]^2/lpn$sigma2[j]^2,scale=lpn$sigma2[j]^2/lpn$mu2[j])
    dep.prob[ind.dep] <- dgamma(obs[ind.dep, 3],shape=lpn$mu3[j]^2/lpn$sigma3[j]^2,scale=lpn$sigma3[j]^2/lpn$mu3[j])
    step.prob[ind.step] <- dgamma(obs[ind.step, 4],shape=lpn$mu4[j]^2/lpn$sigma4[j]^2,scale=lpn$sigma4[j]^2/lpn$mu4[j])
    allprobs[,j] <- surf.prob*dive.prob*dep.prob*step.prob
  }
  foo <- lpn$delta  
  lscale <- 0
  lscale <- forwardalgo(foo, gamma, allprobs, lscale, n)
  return(-lscale)
}

mle <- function(obs,mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,gamma0,N){
  parvect <- pn2pw(mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,gamma0,N)
  mod <- nlm(mllk,parvect,obs,N,print.level=2,iterlim=1000,stepmax=5)
  pn <- pw2pn(mod$estimate,N)
  return(list(mu1=pn$mu1,mu2=pn$mu2,mu3=pn$mu3,mu4=pn$mu4,sigma1=pn$sigma1,sigma2=pn$sigma2,sigma3=pn$sigma3,sigma4=pn$sigma4,gamma=pn$gamma,delta=pn$delta,mllk=mod$minimum))
}

# Beispiel
mod <- mle(obs, c(80, 160), c(210, 190), c(13, 17),c(30,40), c(10, 20),c(15,30),c(5,3),c(10,20), matrix(rep(c(0.5, 0.5), 2), nrow=2, byrow = T), 2)
