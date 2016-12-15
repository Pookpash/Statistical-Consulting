setwd("C:/Users/Svenja/Uni/Master StatWi/WS 16-17/Statistical Consulting/Robben")
seal_clean <- read.csv("seal_data_cleaned.csv")
seal1 <- seal_clean[which(seal_clean$sealID=="1"),]

mle <- function(obs,mu0,sigma0,gamma0,N){
    parvect <- pn2pw(mu0,sigma0,gamma0,N)
    mod <- nlm(mllk,parvect,obs,N,print.level=2,iterlim=1000,stepmax=5)
    pn <- pw2pn(mod$estimate,N)
    return(list(mu=pn$mu,sigma=pn$sigma,gamma=pn$gamma,delta=pn$delta,mllk=mod$minimum))
}

## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
pn2pw <- function(mu,sigma,gamma,N){
    tmu <- log(mu)
    tsigma <- log(sigma)
    tgamma <- NULL
    if(N>1){
        foo <- log(gamma/diag(gamma))           ### gamma ist eine Matrix!
        tgamma <- as.vector(foo[!diag(N)])
    }
    parvect <- c(tmu,tsigma,tgamma)
    return(parvect)
}
## function that performs the inverse transformation
pw2pn <- function(parvect,N){
    mu <- exp(parvect[1:N])
    sigma <- exp(parvect[(N+1):(2*N)])
    gamma <- diag(N)
    if(N>1){
        gamma[!gamma] <- exp(parvect[(2*N+1):(N+N^2)])
        gamma <- gamma/apply(gamma,1,sum)                 ### gamma = Matrix
    }
    delta <- solve(t(diag(N)-gamma+1),rep(1,N))
    return(list(mu=mu,sigma=sigma,gamma=gamma,delta=delta))
}


obs <- list()
for (i in 1:21) {
    obs[[i]]<-seal_clean[which(seal_clean$sealID==i),17]
    }

mllk <- function(parvect,obs,N){
    vec <- rep(NA, length(obs))
    for(i in 1:length(obs)) {
        lpn <- pw2pn(parvect,N)
        gamma<-lpn$gamma
        n <- length(obs[[i]])
        allprobs <- matrix(rep(1,N*n),nrow=n)
        ind.step <- which(!is.na(obs[[i]]) & obs[[i]]!=0)
        for (j in 1:N){
            step.prob <- rep(1,n)
            step.prob[ind.step] <- dgamma(obs[[i]][ind.step],shape=lpn$mu[j]^2/lpn$sigma[j]^2,scale=lpn$sigma[j]^2/lpn$mu[j])
            allprobs[,j] <- step.prob
        }
        foo <- lpn$delta  
        lscale <- 0
        lscale <- forwardalgo(foo, gamma, allprobs, lscale, n)
        vec[i] <- lscale
    }
    lk <- sum(vec)
    return(-lk)
}

N <- 2
parvect <- pn2pw( c(50, 150), c(20, 200), matrix(rep(0.5, 4), nrow = 2),N)
mod <- nlm(mllk,parvect,obs=obs,N,print.level=2,iterlim=1000,stepmax=5)
# zurück convertieren!!!
