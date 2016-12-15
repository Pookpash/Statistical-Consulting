setwd("C:/Users/Pook/Documents/Robben")
seal_clean <- read.csv("seal_data_cleaned.csv")

obs <- cbind(seal_clean$sealID, seal_clean$surf.dur, seal_clean$dive.dur,
             seal_clean$max.dep, seal_clean$steplen)
colnames(obs) <- c("sealID","surf.dur", "dive.dur", "max.dep", "steplen")

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
        obslist <- list()
        for (i in 1:21) {
                obslist[[i]]<-obs[which(obs[,1]==i),]
        }
        return(obslist)
}

#example

# Beispiel
mod <- mle(obs, c(93.93751, 43.75735), c(95.60691, 170.46261), c(7.363427, 22.342912),c(72.04225, 107.97769), c(113.34030,  10.81537),c(76.22124, 41.65749),c(4.926802, 2.656713),c(78.01318, 68.72758), matrix(rep(c(0.5, 0.5), 2), nrow=2, byrow = T), 2)
mod3 <- mle(obs,c(runif(2,40,200)),c(runif(2,40,500)),c(runif(2,5,30)),c(runif(2,10,200)),
            c(runif(2,2,40)),c(runif(2,5,50)),c(runif(2,1,10)),c(runif(2,3,200)),
            matrix(rep(c(0.5, 0.5), 2), nrow=2, byrow = T), 2)
aic.mod(mod,n_variables=4)
bic.mod(mod,n_variables=4,len=100)

#Code to run a specific N-state model multiple times with different starting values (with all 4 variables, adjust accordingly)
fitmult <- function(obs,n_fits,N){
        modl <- list()
        for (i in 1:n_fits){
                mat <- matrix(runif(N^2,0,1), nrow = N)
                mat <- mat/apply(mat,1,sum) 
                modl[[i]] <- mle(obs,c(runif(N,20,200)),c(runif(N,40,300)),c(runif(N,4,40)),c(runif(N,10,200)),
                             c(runif(N,2,200)),c(runif(N,20,70)),c(runif(N,1,20)),c(runif(N,10,150)),
                             mat, N)
        }
        return(modl)
}
results <- fitmult(obs,2,2)

