library(CircStats) # for von Mises distribution
# besser circular statt CircStats?
library(boot) # for logit

setwd("C:/Users/Admin/Documents/Robben")
data <- read.csv("subsample_10.csv")


# theta(mean) & phi(concentration) of Beta distr.
# ttheta <- logit(theta)
# tphi <- log(phi)
# theta <- inv.logit(parvect[(8*N+1):(9*N)])
# phi <- exp(parvect[(9*N+1):(10*N)])


## function that converts 'natural' parameters (possibly constrained) to 'working' parameters (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
# mu & kappa: von Mises distr.
pn2pw <- function(mu1,mu2,mu3,mu4,sigma1,sigma2,sigma3,sigma4,mu,kappa,gamma,N){   
        mu.vonmises <- mu
        for(i in 1:4){
                mu <- cbind(mu1, mu2, mu3, mu4)
                assign(paste0("tmu", i), log(mu[,i]))    
        }
        for(i in 1:4){
                sigma <- cbind(sigma1, sigma2, sigma3, sigma4)
                assign(paste0("tsigma", i), log(sigma[,i]))    
        } 
        tmu <-  kappa*cos(mu.vonmises) # https://github.com/TheoMichelot/moveHMM/blob/master/R/n2w.R
        tkappa <- kappa*sin(mu.vonmises)   
        tgamma <- NULL
        if(N>1){
                foo <- log(gamma/diag(gamma))
                tgamma <- as.vector(foo[!diag(N)])
        }
        parvect <- c(tmu1,tmu2,tmu3,tmu4,tsigma1,tsigma2,tsigma3,tsigma4,tmu,tkappa,tgamma)
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
        x <- parvect[(4*2*N+1):(9*N)]   # https://github.com/TheoMichelot/moveHMM/blob/master/R/w2n.R
        y <- parvect[(9*N+1):(10*N)]   
        mu <- Arg(x+1i*y)
        kappa <- sqrt(x^2+y^2)
        gamma <- diag(N)
        if(N>1){
                gamma[!gamma] <- exp(parvect[(10*N+1):(length(parvect))])
                gamma <- gamma/apply(gamma,1,sum)           
        }
        delta <- solve(t(diag(N)-gamma+1),rep(1,N))
        return(list(mu1=mu1,mu2=mu2,mu3=mu3,mu4=mu4,sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,sigma4=sigma4,
                    mu=mu,kappa=kappa,gamma=gamma,delta=delta))
}

conv2mat <- function(plist, N, var = F){
        plist <- plist[1:(length(plist)-2)]
        if (var == T){
                vecv <- c(plist$sigma1,plist$sigma2,plist$sigma3,plist$kappa,plist$sigma4)
                mat <- matrix(vecv, ncol = N, nrow = (length(plist)/2),byrow=T)
        } else {
                vecm <- c(plist$mu1,plist$mu2,plist$mu3,plist$mu,plist$mu4)
                mat <- matrix(vecm, ncol = N, nrow = (length(plist)/2),byrow=T)
        return(mat)
        }
}

mllk <- function(parvect,obsl,N){ #obsl has to be a list!
        vec <- rep(NA,length(obsl))
        lpn <- pw2pn(parvect,N)
        mumat <- t(conv2mat(lpn, N, var = F))
        sigmat <- t(conv2mat(lpn, N, var = T))
        gamma <- lpn$gamma
        for(i in 1:length(obsl)){
                n <- length(obsl[[i]][,1])
                allprobs <- allprobs_rcpp(N,n,as.matrix(obsl[[i]][,c(13,14,7,6,5)]),mumat,sigmat)
                foo <- lpn$delta  
                lscale <- 0
                lscale <- forwardalgo(foo, gamma, allprobs, lscale, n)
                vec[i] <- lscale
        }
        lk <- sum(vec)
        return(-lk)
}

mle <- function(obs,mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,mu0,kappa0,gamma0,N){
        parvect <- pn2pw(mu01,mu02,mu03,mu04,sigma01,sigma02,sigma03,sigma04,mu0,kappa0,gamma0,N)
        obsl <- create_obslist(obs)
        mod <- nlm(mllk,parvect,obsl,N,print.level=2,iterlim=1000,stepmax=5)
        pn <- pw2pn(mod$estimate,N)
        return(list(mu1=pn$mu1,mu2=pn$mu2,mu3=pn$mu3,mu4=pn$mu4,sigma1=pn$sigma1,sigma2=pn$sigma2,sigma3=pn$sigma3,sigma4=pn$sigma4,
                    mu=pn$mu,kappa=pn$kappa,gamma=pn$gamma,delta=pn$delta,mllk=mod$minimum))
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
                modl[[i]] <- mle(obs,c(runif(N,3,45)),c(runif(N,80,200)),c(runif(N,25,60)),c(runif(N,5,200)),c(runif(N,1,10)),
                                 c(runif(N,25,100)),c(runif(N,5,50)),c(runif(N,10,100)),c(runif(N,-3,3)),c(runif(N,1,10)),
                                 mat, N)
        }
        return(modl)
}

#test

fitmult(data,1,3)

#troubleshooting

N <- 2 
obsl <- create_obslist(data)
mat <- matrix(runif(N^2,0,1), nrow = N)
mat <- mat/apply(mat,1,sum) 
parvect <- pn2pw(c(runif(N,3,35)),c(runif(N,80,175)),c(runif(N,30,55)),c(runif(N,10,100)),c(runif(N,5,10)),
                 c(runif(N,25,100)),c(runif(N,10,50)),c(runif(N,10,100)),c(runif(N,0.1,1)),c(runif(N,4,7)),
                 mat,N)

vec <- rep(NA,length(obsl))
lpn <- pw2pn(parvect,N)
mumat <- t(conv2mat(lpn, N, var = F))
sigmat <- t(conv2mat(lpn, N, var = T))
gamma <- lpn$gamma

i <- 1

n <- length(obsl[[i]][,1])
allprobs <- allprobs_rcpp(N,n,as.matrix(obsl[[i]][,c(13,14,7,6,5)]),mumat,sigmat)#gives NaNs, and therefore loglike -inf later!
allprobs
foo <- lpn$delta
lscale <- 0
lscale <- forwardalgo(foo, gamma, allprobs, lscale, n)
vec[i] <- lscale
