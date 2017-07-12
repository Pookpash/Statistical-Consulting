library(CircStats) # for von Mises distribution
library(boot) # for logit

setwd("C:/Users/Admin/Documents/Robben")
seal10 <- read.csv("subsample_10scaled.csv",header=T)
seal10 <- seal10[,-1] #aus irgendeinem Grund eine Spalte zuviel...

## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
# mu & kappa: von Mises distr.
# beta: Parameter in folgender Reihenfolge angeben: erst alle beta0, dann alle beta1, ...
pn2pw <- function(mu1,mu2,mu3,mu4,sigma1,sigma2,sigma3,sigma4,kappa,delta,beta,N){   
    for(i in 1:4){
        mu <- cbind(mu1, mu2, mu3, mu4)
        assign(paste0("tmu", i), log(mu[,i]))    
    }
    for(i in 1:4){
        sigma <- cbind(sigma1, sigma2, sigma3, sigma4)
        assign(paste0("tsigma", i), log(sigma[,i]))    
    } 
    tkappa <- log(kappa)
    tbeta <- matrix(beta,byrow=T,ncol=N*(N-1))
    tdelta <- log(delta)
    parvect <- c(tmu1,tmu2,tmu3,tmu4,tsigma1,tsigma2,tsigma3,tsigma4,tkappa,tdelta,tbeta)
    return(parvect)
}

## function that performs the inverse transformation
pw2pn <- function(parvect,N){
    mu <- exp(parvect[1:(4*N)])
    for(i in 1:4) {
        assign(paste0("mu", i), mu[(i*N-N+1):(i*N)])
    }   
    sigma <- exp(parvect[(4*N+1):(4*2*N)])
    for(i in 1:4) {
        assign(paste0("sigma", i), sigma[(i*N-N+1):(i*N)])
    } 
    kappa <- exp(parvect[(4*2*N+1):(9*N)])
    delta <- exp(parvect[(9*N+1):(10*N)])
    delta <- delta/sum(delta)
    beta <- matrix(parvect[(10*N+1):length(parvect)],ncol=N*(N-1),byrow=T)
    return(list(mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,sigma4=sigma4,
                kappa=kappa,delta=delta,beta=beta))
}

conv2mat <- function(plist, N, var = F){
    if (var == T){
        vecv <- c(plist$sigma1,plist$sigma2,plist$sigma3,plist$kappa,plist$sigma4)
        mat <- matrix(vecv, ncol = N, nrow = 5,byrow=T) # falls surf.dur entfernt wird, Spaltenanzahl -1
    } else {
        vecm <- c(plist$mu1,plist$mu2,plist$mu3,plist$mu4)
        mat <- matrix(vecm, ncol = N, nrow = 4,byrow=T) # falls surf.dur entfernt wird, Spaltenanzahl -1
        return(mat)
    }
}


mllk <- function(parvect,obsl,N,covs){ #covs as a a vector of columns, e.g. c(1,2,5)]
    vec <- rep(NA,length(obsl))
    lpn <- pw2pn(parvect,N)
    foo <- lpn$delta
    mumat <- t(conv2mat(lpn, N, var = F))
    sigmat <- t(conv2mat(lpn, N, var = T))
    for(i in 1:length(obsl)){
        n <- length(obsl[[i]][,1])
        covs.mat <- matrix(c(rep(1,dim(obsl[[i]])[1]),obsl[[i]][,covs]),ncol=length(covs)+1,byrow = F)
        gamma <- trMatrix_rcpp(N, lpn$beta, covs.mat)
        allprobs <- allprobs_rcpp(N,n,as.matrix(obsl[[i]][,c(13,14,7,6,5)]),mumat,sigmat) 
        lscale <- 0
        lscale <- forwardalgo_w_cov(foo,gamma, allprobs, lscale, n)
        vec[i] <- lscale
    }
    lk <- sum(vec)
    return(-lk)
}

create_obslist <- function(obs){
    obslist <- list()
    for (i in 1:length(unique(obs[,"ID_burst"]))) {
        obslist[[i]]<-obs[which(obs[,"ID_burst"]==unique(obs[,"ID_burst"])[i]),]
    }
    obslist <- Filter(length,obslist)
    return(obslist)
}

mle <- function(obs,mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,kappa0,delta0,beta0,covs,N){
    parvect <- pn2pw(mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,kappa0,delta0,beta0,N)
    obsl <- create_obslist(obs)
    mod <- nlm(mllk,parvect,obsl,N,covs,hessian=T,print.level=2,iterlim=1000,stepmax=5)
    pn <- pw2pn(mod$estimate,N)
    return(list(mu1=pn$mu1,mu2=pn$mu2,mu3=pn$mu3,mu4=pn$mu4,sigma1=pn$sigma1,sigma2=pn$sigma2,sigma3=pn$sigma3,sigma4=pn$sigma4,
                kappa=pn$kappa,delta=pn$delta,beta=pn$beta,mllk=mod$minimum,hessmat=mod$hessian))
}

viterbi<-function(obs,mod,covs,N){ #mod h.t.b. defined as mod[[x]] otherwise the $-operator does not work
    beta <- mod$beta
    delta <- mod$delta
    mumat <- t(conv2mat(mod,N,F))
    sigmat <- t(conv2mat(mod,N,T))
    allStates <- NULL
    obsl <- create_obslist(obs)
    for(i in 1:length(obsl)){
        n <- length(obsl[[i]][,1])
        covs.mat <- matrix(c(rep(1,dim(obsl[[i]])[1]),obsl[[i]][,covs]),ncol=length(covs)+1,byrow = F)
        gamma <- trMatrix_rcpp(N, beta, covs.mat)
        allprobs <- allprobs_rcpp(N,n,as.matrix(obsl[[i]][,c(13,14,7,6,5)]),mumat,sigmat)
        
        xi <- matrix(0,as.integer(n),N)
        u <- delta%*%diag(N)
        xi[1,] <- u/sum(u)
        for (t in 2:n){
            u<-apply(xi[t-1,]*gamma[,,t],2,max)%*%diag(allprobs[t,])
            xi[t,] <- u/sum(u)
        }
        iv<-numeric(n)
        iv[n] <-which.max(xi[n,])
        for (t in (n-1):1){ 
            iv[t] <- which.max(gamma[,iv[t+1],t+1]*xi[t,]) #slice "t+1" correct?
        }
        allStates <- c(allStates,iv)
    }
    return(allStates)
}

confints <- function(mod){ #first just mu1 to test
    muvec <- mod$mu1
    fisher_info <- solve(-mod$hessmat)
    prop_sigma<-sqrt(diag(fisher_info))
    prop_sigma<-diag(prop_sigma)
    upper<-muvec+1.96*prop_sigma[1:3]
    lower<-muvec-1.96*prop_sigma[1:3]
    interval<-data.frame(upper=upper, lower=lower)
}

#Code to run a specific N-state model multiple times with different starting values (with all 4 variables, adjust accordingly)
# n_cov covariates and intercept
fitmult <- function(obs,n_fits,covs,N){ #covs as a a vector of columns, e.g. c(1,2,5)]
    modl <- list()
    for (i in 1:n_fits){
        modl[[i]] <- mle(obs,mu01=c(runif(N,3,35)),mu02=c(runif(N,80,175)),mu03=c(runif(N,30,55)),mu04=c(runif(N,10,150)),
                         sigma01=c(runif(N,2,10)),sigma02=c(runif(N,25,100)),sigma03=c(runif(N,10,30)),sigma04=c(runif(N,10,80)),
                         kappa0=c(runif(N,2,7)),delta0 = c(rep(1,N)),beta0=c(runif(N*(N-1),-3,-1),rnorm(length(covs)*N*(N-1))),covs,N)
        assign(paste0("mod",i),modl[[i]],env = .GlobalEnv) # save each model
    }
    return(modl)
}

test<-fitmult(seal10[1:1000,],1,c(8),3)

testvit <- viterbi(seal10[1:500,],test[[1]],c(8),3)

testconf <- confints(test[[1]]) #hessian almost never invertible and if so NaN bc of sqrt(-...)
                                #so i think hessian has also to be transformed back maybe?
