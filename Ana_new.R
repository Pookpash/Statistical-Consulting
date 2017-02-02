setwd("C:/Users/Pook/Documents/Robben")
seal_clean <- read.csv("seal_data_cleaned.csv")

seal_clean$X <- 0

obs <- cbind(seal_clean$sealID, seal_clean$surf.dur, seal_clean$dive.dur,
             seal_clean$max.dep, seal_clean$steplen,seal_clean$X)
colnames(obs) <- c("sealID","surf.dur", "dive.dur", "max.dep", "steplen","div.ID")

#split dataset for harbour and grey seals (1=grey seal)
obs_1 <- obs[which(obs[,1]<=11),]
obs_0 <- obs[which(obs[,1]>=12),]
for (i in 12:21){
  obs_0[,1][obs_0[,1]==i]<-(i-11)
}

set_divID<- function(observs){
  observs[,5][observs[,2]>=600] <- 0
  observs[,5][observs[,3]>=600] <- 0
  observs[,5][observs[,5]>=500] <- 0
  curdiv <- 1 #init diveID
  
  for(i in 1:dim(observs)[1]){
    if(observs[i,5]>0){
      observs[i,6]<-curdiv
    }else{
      curdiv <- curdiv + 1
      observs[i,6]<-curdiv
    }
  }
  return(observs)
}

obs_1 <- set_divID(observs = obs_1)
obs_0 <- set_divID(observs = obs_0)

obs_0 <- obs_0[obs_0[,5] > 0,]
obs_0[,6]<- obs_0[,6]-2 # fix for the harbour seals since they start at dive 3
obs_1 <- obs_1[obs_1[,5] > 0,]

fix_it2 <- function(observs){ #remove one data point time series
  for (i in unique(observs[,6])){
    if(sum(observs[,6]==i)<2){
      val<-which(observs[,6]==i)
      observs[val,6] <- 0
    }
  }
  observs <- observs[observs[,6]>0,]
  return(observs)
}

obs_1<- fix_it2(obs_1)
obs_0<- fix_it2(obs_0)

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
  lmin <- min(obs[,6])
  lmax <- max(obs[,6])
  obslist <- list()
  for (i in lmin:lmax) {
    obslist[[i]]<-obs[which(obs[,6]==i),]
  }
  return(obslist)
}

create_obslist <- function(obs){
  obslist <- list()
  for (i in unique(obs[,6])) {
    obslist[[i]]<-obs[which(obs[,6]==i),]
  }
  obslist <- Filter(length,obslist)
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

plotresults <- function(res,N,sealtype,densvec=c(rep(0.1,4))){
  par(mfrow=c(2,2),oma=c(0,0,2,0))
  x1 <- seq(1,150,length=1000)
  x2 <- seq(1,250,length=1000)
  x3 <- seq(0.1,45,length=1000)
  x4 <- seq(1,200,length=1000)
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
  if(sealtype==1){
    title1 <- paste0("Grey Seal state dependent distributions for ", N , " states")
  } else {
    title1 <- paste0("Harbour Seal state dependent distributions for ", N , " states")
  }
  mtext(title1, outer = TRUE, cex = 1.5)
  par(mfrow=c(1,1))
}

viterbi<-function(obs,mod,N){
  Gamma <- mod$gamma
  delta <- mod$delta
  allStates <- NULL
  obsl <- create_obslist(obs)
  for(i in 1:length(obsl)){
    T <- dim(obsl[[i]])[1]
    allprobs <- matrix(rep(1,N*T),nrow=T)
    ind.surf <- which(!is.na(obsl[[i]][,2]))
    ind.dive <- which(!is.na(obsl[[i]][,3]))
    ind.dep <- which(!is.na(obsl[[i]][,4]))
    ind.step <- which(!is.na(obsl[[i]][,5]))
    for (j in 1:N){
      surf.prob <- dive.prob <- dep.prob <- step.prob <- rep(1,T)
      surf.prob[ind.surf] <- dgamma(obsl[[i]][ind.surf, 2],shape=mod$mu1[j]^2/mod$sigma1[j]^2,scale=mod$sigma1[j]^2/mod$mu1[j])
      dive.prob[ind.dive] <- dgamma(obsl[[i]][ind.dive, 3],shape=mod$mu2[j]^2/mod$sigma2[j]^2,scale=mod$sigma2[j]^2/mod$mu2[j])
      dep.prob[ind.dep] <- dgamma(obsl[[i]][ind.dep, 4],shape=mod$mu3[j]^2/mod$sigma3[j]^2,scale=mod$sigma3[j]^2/mod$mu3[j])
      step.prob[ind.step] <- dgamma(obsl[[i]][ind.step, 5],shape=mod$mu4[j]^2/mod$sigma4[j]^2,scale=mod$sigma4[j]^2/mod$mu4[j])
      allprobs[,j] <- surf.prob*dive.prob*dep.prob*step.prob
    }
    xi <- matrix(0,as.integer(T),N)
    u <- delta%*%diag(N)
    xi[1,] <- u/sum(u)
    for (t in 2:T){
      u<-apply(xi[t-1,]*Gamma,2,max)%*%diag(allprobs[t,])
      xi[t,] <- u/sum(u)
    }
    iv<-numeric(T)
    iv[T] <-which.max(xi[T,])
    for (t in (T-1):1){ 
      iv[t] <- which.max(Gamma[,iv[t+1]]*xi[t,])
    }
    allStates <- c(allStates,iv)
  }
  return(allStates)
}

plot_viterbi <- function(stateobj,sealtype,obs,mod,N){
  states_surf <- as.data.frame(stateobj)
  states_surf$col <- states_surf[,1]+7
  states_div <- as.data.frame(stateobj)
  states_div$col <- states_div[,1]+7
  states_dep <- as.data.frame(stateobj)
  states_dep$col <- states_dep[,1]+7
  states_step <- as.data.frame(stateobj)
  states_step$col <- states_step[,1]+7
  for(i in 1:N){
    states_surf[,1][states_surf[,1]==i]<- mod$mu1[i]
  }
  for(i in 1:N){
    states_div[,1][states_div[,1]==i]<- mod$mu2[i]
  }
  for(i in 1:N){
    states_dep[,1][states_dep[,1]==i]<- mod$mu3[i]
  }
  for(i in 1:N){
    states_step[,1][states_step[,1]==i]<- mod$mu4[i]
  }
  par(mfrow=c(2,2),oma=c(0,0,2,0))
  plot(obs[,2],type="l",col="blue",lwd=2, ylim= c(0,200), xlab = "Observation no.", ylab= "seconds",main="Surface Duration")
  points(states_surf[,1],pch=15,col=states_surf$col)
  plot(obs[,3],type="l",col="blue",lwd=2,xlab = "Observation no.", ylab= "seconds",main="Dive Duration")
  points(states_div[,1],pch=15,col=states_div$col)
  plot(obs[,4],type="l",col="blue",lwd=2,xlab = "Observation no.", ylab= "meters",main="Dive Depth")
  points(states_dep[,1],pch=15,col=states_dep$col)
  plot(obs[,5],type="l",col="blue",lwd=2,xlab = "Observation no.", ylab= "meters",main="Steplength")
  points(states_step[,1],pch=15, col=states_step$col)
  if(sealtype==1){
    title1 <- paste0("Grey Seal 1 decoded states (N = ", N,", first 250 obs.)")
  } else {
    title1 <- paste0("Harbour Seal 1 decoded states (N = ", N,", first 250 obs.)")
  }
  mtext(title1, outer = TRUE, cex = 1.5)
  par(mfrow=c(1,1))
}

results <- fitmult(obs_0,1,2)
states <- viterbi(obs = obs_0,mod = results[[1]],N=2)
plot_viterbi(states[1:250],sealtype=0,obs=obs_0[1:250,],mod=results[[1]],N=2)
