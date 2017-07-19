setwd("C:/Users/User/Documents/Studium_MA/4. Semester/Statistical Consulting/R")

library(CircStats) # for von Mises distribution

lforward <- function(data,mod,N,covs){
  obsl <-  create_obslist(data)
  lalpha.all <- rep(NA,sum(sapply(obsl,nrow)))
  lalpha.list <- list()
  mumat <- t(conv2mat(mod,N,var=F))
  sigmat <- t(conv2mat(mod,N,var=T))
  for (i in 1:length(obsl)){
    n <- dim(obsl[[i]])[1]
    lalpha <- matrix(NA,n,N)
    allprobs <- allprobs_rcpp(N,n,as.matrix(obsl[[i]][,c(13,14,7,6,5)]),mumat,sigmat)
    covsvec <- covsfix(obsl[[i]],covs)
    covsvec <- c(rep(1,n),covsvec)
    covs.mat <- matrix(covsvec,ncol=length(covs)+1,byrow = F)
    gamma <- trMatrix_rcpp(N, mod$beta, covs.mat)
    u <- mod$delta*allprobs[1,]
    l <- log(sum(u))
    phi <- u/sum(u)   
    lalpha[1,] <- l+log(phi)
    for(t in 2:n) {
      u <- phi%*%gamma[,,t]*allprobs[t,]
      l <- l+log(sum(u))
      phi <- u/sum(u)
      lalpha[t,] <- l+log(phi)
    }
    lalpha.list[[i]] <- lalpha
  }
  return(lalpha.list)
}

psres <- function(data,mod,N,covs){
  obsl <- create_obslist(data)
  surf.res.list <- c()
  dive.res.list <- c()
  max.res.list <- c()
  step.res.list <- c()
  angle.res.list <- c()
  for (j in 1:length(obsl)){
    obs <- obsl[[j]]
    n <- dim(obs)[1]
    covsvec <- covsfix(obsl[[j]],covs)
    covsvec <- c(rep(1,n),covsvec)
    covs.mat <- matrix(covsvec,ncol=length(covs)+1,byrow = F)
    gamma <- trMatrix_rcpp(N, mod$beta, covs.mat)
    surf.res <- rep(NA,n)
    dive.res <- rep(NA,n)
    max.res <- rep(NA,n)
    step.res <- rep(NA,n)
    angle.res <- rep(NA,n)
    psurf <- matrix(NA,n,N)
    pdive <- matrix(NA,n,N)
    pmax <- matrix(NA,n,N)
    pstep <- matrix(NA,n,N)
    pangle <- matrix(NA,n,N)
    la <- lforward(data,mod,N,covs)
    for(i in 1:N) {
      for(t in 1:n) {
        psurf[t,i] <- pgamma(obs[t,"surf.dur"],shape=mod$mu1[i]^2/mod$sigma1[i]^2,scale=mod$sigma1[i]^2/mod$mu1[i])
        pdive[t,i] <- pgamma(obs[t,"dive.dur"],shape=mod$mu2[i]^2/mod$sigma2[i]^2,scale=mod$sigma2[i]^2/mod$mu2[i])
        pmax[t,i] <- pgamma(obs[t,"max.dep"],shape=mod$mu3[i]^2/mod$sigma3[i]^2,scale=mod$sigma3[i]^2/mod$mu3[i])
        pstep[t,i] <- pgamma(obs[t,"steplen"],shape=mod$mu4[i]^2/mod$sigma4[i]^2,scale=mod$sigma4[i]^2/mod$mu4[i])
        pangle[t,i] <- pvm(obs[t,"turnA"],0,kappa=mod$kappa[i])
      }
    }
    surf.res[1] <- qnorm(mod$delta%*%psurf[1,]) 
    dive.res[1] <- qnorm(mod$delta%*%pdive[1,])
    max.res[1] <- qnorm(mod$delta%*%pmax[1,])
    step.res[1] <- qnorm(mod$delta%*%pstep[1,])
    angle.res[1] <- qnorm(mod$delta%*%pangle[1,])
    for(t in 2:n) {
      lat <- la[[j]][t-1,]  
      c <- max(lat) 
      elat <- exp(lat-c)
      surf.res[t] <- qnorm(t(elat)%*%(gamma[,,t]/sum(elat))%*%psurf[t,])
      dive.res[t] <- qnorm(t(elat)%*%(gamma[,,t]/sum(elat))%*%pdive[t,])
      max.res[t] <- qnorm(t(elat)%*%(gamma[,,t]/sum(elat))%*%pmax[t,])
      step.res[t] <- qnorm(t(elat)%*%(gamma[,,t]/sum(elat))%*%pstep[t,])
      angle.res[t] <- qnorm(t(elat)%*%(gamma[,,t]/sum(elat))%*%pangle[t,])
    }
    surf.res.list <- c(surf.res.list,surf.res)
    dive.res.list <- c(dive.res.list,dive.res)
    max.res.list <- c(max.res.list,max.res)
    step.res.list <- c(step.res.list,step.res)
    angle.res.list <- c(angle.res.list,angle.res)
  }
  
  return(list(surf.residuals=surf.res.list,dive.residuals=dive.res.list,
              max.residuals=max.res.list,step.residuals=step.res.list,
              angle.residuals=angle.res.list))
}

psres(seal10[1:10000,],mod1,3,c(8))
