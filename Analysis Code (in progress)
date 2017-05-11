library(CircStats) # for von Mises distribution
library(boot) # for logit

## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
# mu & kappa: von Mises distr.; theta(mean) & phi(concentration) of Beta distr.
pn2pw <- function(mu1,mu2,mu3,sigma1,sigma2,sigma3,mu,kappa,theta,phi,gamma,N){   
  for(i in 1:3){
    mu <- cbind(mu1,mu2,mu3)
    assign(paste0("tmu", i), log(mu[,i]))    
  }
  for(i in 1:3){
    sigma <- cbind(sigma1,sigma2,sigma3)
    assign(paste0("tsigma", i), log(sigma[,i]))    
  }
  tmu <- NULL   # Transformation fehlt noch!
  tkappa <- NULL   # Transformation fehlt noch!
  ttheta <- logit(theta)
  tphi <- log(phi)
  tgamma <- NULL
  if(N>1){
    foo <- log(gamma/diag(gamma))
    tgamma <- as.vector(foo[!diag(N)])
  }
  parvect <- c(tmu1,tmu2,tmu3,tsigma1,tsigma2,tsigma3,tmu,tkappa,ttheta,tphi,tgamma)
  return(parvect)
}

## function that performs the inverse transformation
pw2pn <- function(parvect,N){
  mu <- exp(parvect[1:(3*N)])
  for(i in 1:3) {
    assign(paste0("mu", i), mu[(i*N-N+1):(i*N)])
  }   
  sigma <- exp(parvect[(3*N+1):(3*2*N)])
  for(i in 1:4) {
    assign(paste0("sigma", i), sigma[(i*N-N+1):(i*N)])
  } 
  mu <- parvect[(3*2*N+1):(7*N)]   # Transformation fehlt noch!
  kappa <- parvect[(7*N+1):(8*N)]   # Transformation fehlt noch!
  theta <- inv.logit(parvect[(8*N+1):(9*N)])
  phi <- exp(parvect[(9*N+1):(10*N)])
  gamma <- diag(N)
  if(N>1){
    gamma[!gamma] <- exp(parvect[(10*N+1):(length(parvect))])
    gamma <- gamma/apply(gamma,1,sum)           
  }
  delta <- solve(t(diag(N)-gamma+1),rep(1,N))
  return(list(mu1=mu1,mu2=mu2,mu3=mu3,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,
              mu=mu,kappa=kappa,theta=theta,phi=phi,gamma=gamma,delta=delta))
}

mllk <- function(parvect,obsl,N){ #obsl has to be a list!
  vec <- rep(NA,length(obsl))
  for(i in 1:length(obsl)){
    lpn <- pw2pn(parvect,N)
    gamma <- lpn$gamma
    n <- length(obsl[[i]][,1])
    allprobs <- matrix(rep(1,N*n),nrow=n)
    ind.surf <- which(!is.na(obsl[[i]][,2]))
    ind.dive <- which(!is.na(obsl[[i]][,3]))
    ind.dep <- which(!is.na(obsl[[i]][,4]))
    ind.step <- which(!is.na(obsl[[i]][,5]))
    for (j in 1:N){
      surf.prob <- dive.prob <- dep.prob <- step.prob <- rep(1,n)
      surf.prob[ind.surf] <- dgamma(obsl[[i]][ind.surf, 2],shape=lpn$mu1[j]^2/lpn$sigma1[j]^2,scale=lpn$sigma1[j]^2/lpn$mu1[j])
      dive.prob[ind.dive] <- dgamma(obsl[[i]][ind.dive, 3],shape=lpn$mu2[j]^2/lpn$sigma2[j]^2,scale=lpn$sigma2[j]^2/lpn$mu2[j])
      dep.prob[ind.dep] <- dgamma(obsl[[i]][ind.dep, 4],shape=lpn$mu3[j]^2/lpn$sigma3[j]^2,scale=lpn$sigma3[j]^2/lpn$mu3[j])
      step.prob[ind.step] <- dgamma(obsl[[i]][ind.step, 5],shape=lpn$mu4[j]^2/lpn$sigma4[j]^2,scale=lpn$sigma4[j]^2/lpn$mu4[j])
      allprobs[,j] <- surf.prob*dive.prob*dep.prob*step.prob
    }
    foo <- lpn$delta  
    lscale <- 0
    lscale <- forwardalgo(foo, gamma, allprobs, lscale, n)
    vec[i] <- lscale
  }
  lk <- sum(vec)
  return(-lk)
}

mle <- function(obs,mu01,mu02,mu03,sigma01,sigma02,sigma03,mu0,kappa0,theta0,phi0,gamma0,N){
  parvect <- pn2pw(mu01,mu02,mu03,sigma01,sigma02,sigma03,mu0,kappa0,theta0,phi0,gamma0,N)
  obsl <- create_obslist(obs)
  mod <- nlm(mllk,parvect,obsl,N,print.level=2,iterlim=1000,stepmax=5)
  pn <- pw2pn(mod$estimate,N)
  return(list(mu1=pn$mu1,mu2=pn$mu2,mu3=pn$mu3,sigma1=pn$sigma1,sigma2=pn$sigma2,sigma3=pn$sigma3,
              mu=pn$mu,kappa=pn$kappa,theta=pn$theta,phi=pn$phi,gamma=pn$gamma,delta=pn$delta,
              mllk=mod$minimum))
}

create_obslist <- function(obs){
  obslist <- list()
  for (i in 1:length(unique(obs[,"ID_burst"]))) {
    obslist[[i]]<-obs[which(obs[,"ID_burst"]==unique(obs[,"ID_burst"])[i]),]   # so gelöst, da ID_burst nicht 1,2,3... benannt
  }
  obslist <- Filter(length,obslist)   # Was macht der Befehl?
  return(obslist)
}

#Code to run a specific N-state model multiple times with different starting values (with all 4 variables, adjust accordingly)
fitmult <- function(obs,n_fits,N){
  modl <- list()
  for (i in 1:n_fits){
    mat <- matrix(runif(N^2,0,1), nrow = N)
    mat <- mat/apply(mat,1,sum) 
    modl[[i]] <- mle(obs,c(runif(N,20,200)),c(runif(N,40,300)),c(runif(N,1,50)),c(runif(N,10,200)),
                     c(runif(N,2,200)),c(runif(N,20,70)),c(runif(N,1,30)),c(runif(N,10,150)),
                     mat, N)
  }
  return(modl)
}
