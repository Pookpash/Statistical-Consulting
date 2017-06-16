
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
