setwd("C:/Users/Pook/Documents/Robben")
seal_clean <- read.csv("seal_data_cleaned.csv")

obs <- cbind(seal_clean$sealID, seal_clean$surf.dur, seal_clean$dive.dur,
             seal_clean$max.dep, seal_clean$steplen)
colnames(obs) <- c("sealID","surf.dur", "dive.dur", "max.dep", "steplen")
#split dataset for harbour and grey seals (1=grey seal)
obs_1 <- obs[which(obs[,1]<=11),]
obs_0 <- obs[which(obs[,1]>=12),]

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
    assign(paste0("mu", i), mu[(i*N-N+1):(i*N)])    ####### geänderte Indizierung!
  }   
  sigma <- exp(parvect[(4*N+1):(4*2*N)])
  for(i in 1:4) {
    assign(paste0("sigma", i), sigma[(i*N-N+1):(i*N)])   ####### geänderte Indizierung!
  }   
  gamma <- diag(N)
  if(N>1){
    gamma[!gamma] <- exp(parvect[(4*2*N+1):(length(parvect))])
    gamma <- gamma/apply(gamma,1,sum)                 ### gamma = Matrix
  }
  delta <- solve(t(diag(N)-gamma+1),rep(1,N))
  return(list(mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,sigma4=sigma4,gamma=gamma,delta=delta))
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

mle <- function(obs,mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,gamma0,N){
  parvect <- pn2pw(mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,gamma0,N)
  obsl <- create_obslist(obs)
  mod <- nlm(mllk,parvect,obsl,N,print.level=2,iterlim=1000,stepmax=5)
  pn <- pw2pn(mod$estimate,N)
  return(list(mu1=pn$mu1,mu2=pn$mu2,mu3=pn$mu3,mu4=pn$mu4,sigma1=pn$sigma1,sigma2=pn$sigma2,sigma3=pn$sigma3,sigma4=pn$sigma4,gamma=pn$gamma,delta=pn$delta,mllk=mod$minimum))
}

#aic and bic only work for pdfs with two parameters! good enough for us atm.
aic.mod <- function(mod,n_variables){
  llk <- mod$mllk
  no_states <- length(mod$gamma[1,])
  params <- 2*n_variables*no_states+no_states*(no_states-1)
  aic <- 2*llk+2*params
  return(aic)
}

#aic and bic only work for pdfs with two parameters! good enough for us atm.
bic.mod <- function(mod,n_variables,len){ #len should be the length of the data
  llk <- mod$mllk
  no_states <- length(mod$gamma[1,])
  params <- 2*n_variables*no_states+no_states*(no_states-1)
  bic <- 2*llk+log(len)*params
  return(bic)
}

create_obslist <- function(obs){
  lmin <- min(obs[,1])
  lmax <- max(obs[,1])
  obslist <- list()
  for (i in lmin:lmax) {
    obslist[[i]]<-obs[which(obs[,1]==i),]
  }
  fixval <- length(obslist)#fix for the problem with the harbour seals
  if(fixval>=15){ 
    obslist <- obslist[-c(1:11)]
  }
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
results <- fitmult(obs_0,1,2)

#WIP
plotresults <- function(res,N,densvec=c(rep(0.1,4))){
  par(mfrow=c(2,2))
  x1 <- seq(1,150,length=1000)
  x2 <- seq(1,250,length=1000)
  x3 <- seq(0.1,40,length=1000)
  x4 <- seq(1,180,length=1000)
  plot(x1,dgamma(x1,res$mu1[1]^2/res$sigma1[1]^2,scale=res$sigma1[1]^2/res$mu1[1]), type="l", col=1, ylim=c(0,densvec[1]), main="Surface Duration", ylab="density", xlab= "seconds", lwd=2)
  for(i in 2:N){
    points(x1,dgamma(x1,res$mu1[i]^2/res$sigma1[i]^2,scale=res$sigma1[i]^2/res$mu1[i]), type="l", col=i, lwd=2)
  }
  
  plot(x2,dgamma(x2,res$mu2[1]^2/res$sigma2[1]^2,scale=res$sigma2[1]^2/res$mu2[1]), type="l", col=1, ylim=c(0,densvec[2]), main="Dive Duration", ylab="density", xlab= "seconds", lwd=2)
  for(i in 2:N){
    points(x2,dgamma(x2,res$mu2[i]^2/res$sigma2[i]^2,scale=res$sigma2[i]^2/res$mu2[i]), type="l", col=i, lwd=2)
  }
  
  plot(x3,dgamma(x3,res$mu3[1]^2/res$sigma3[1]^2,scale=res$sigma3[1]^2/res$mu3[1]), type="l", col=1, ylim=c(0,densvec[3]), main="Maximum Depth", ylab="density", xlab= "meter", lwd=2)
  for(i in 2:N){
    points(x3,dgamma(x3,res$mu3[i]^2/res$sigma3[i]^2,scale=res$sigma3[i]^2/res$mu3[i]), type="l", col=i, lwd=2)
  }
  
  plot(x4,dgamma(x4,res$mu4[1]^2/res$sigma4[1]^2,scale=res$sigma4[1]^2/res$mu4[1]), type="l", col=1, ylim=c(0,densvec[4]), main="Steplength", ylab="density", xlab= "meter", lwd=2)
  for(i in 2:N){
    points(x4,dgamma(x4,res$mu4[i]^2/res$sigma4[i]^2,scale=res$sigma4[i]^2/res$mu4[i]), type="l", col=i, lwd=2)
  }
  par(mfrow=c(1,1))
}
