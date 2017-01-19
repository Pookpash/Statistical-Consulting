setwd("C:/Users/Admin/Documents/Robben")
seal_clean <- read.csv("seal_data_cleaned.csv")

obs <- cbind(seal_clean$sealID, seal_clean$surf.dur, seal_clean$dive.dur,
             seal_clean$max.dep, seal_clean$steplen)
colnames(obs) <- c("sealID","surf.dur", "dive.dur", "max.dep", "steplen")
#split dataset for harbour and grey seals (1=grey seal)
obs_1 <- obs[which(obs[,1]<=11),]
obs_0 <- obs[which(obs[,1]>=12),]
for (i in 12:21){
        obs_0[,1][obs_0[,1]==i]<-(i-11)
}

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

res2sl <- list()
res2sl$mu1<- c(38.63129, 111.38412)
res2sl$mu2<- c(156.91947, 71.18653)
res2sl$mu3<- c(22.29625, 6.10472)
res2sl$mu4<- c(120.10310, 26.85120)
res2sl$sigma1<- c(11.83331, 144.80094)
res2sl$sigma2<- c(57.62893, 62.17002)
res2sl$sigma3<- c(13.64301, 4.78521)
res2sl$sigma4<- c(82.68534, 30.19057)
res2sl$gamma <- matrix(c(0.931261, 0.068739, 0.121455, 0.878546),2,byrow=T)
res2sl$delta <- c(0.6385836, 0.3614164)

res3sl <- list()
res3sl$mu1<- c(141.40020, 46.54668, 36.31873)
res3sl$mu2<- c(39.63630, 148.27180, 156.74800)
res3sl$mu3<- c(4.04023, 17.42329, 22.92530)
res3sl$mu4<- c(19.53895, 52.51476, 157.32574)
res3sl$sigma1<- c(201.91241, 25.40830, 8.86741)
res3sl$sigma2<- c(27.81127, 82.01774, 43.82091)
res3sl$sigma3<- c(2.54377, 12.74487, 13.89699)
res3sl$sigma4<- c(21.75779, 43.74214, 66.28669)
res3sl$gamma <- matrix(c(0.824601, 0.129573, 0.045826, 0.089276, 0.845705,
                        0.065019, 0.022303, 0.064740, 0.912958),nrow=3,byrow=T)
res3sl$delta <- c(0.2366674, 0.3656206, 0.3977120)

res4sl <- list()
res4sl$mu1<- c(31.88682,  43.06180,  52.36523, 172.35877)
res4sl$mu2<- c(141.15560, 170.96454, 124.07681,  30.61768)
res4sl$mu3<- c(13.453759, 32.833299, 10.621648,  3.267052)
res4sl$mu4<- c(153.40901, 114.75732,  44.92601,  13.27642)
res4sl$sigma1<- c(6.640114,  10.799599,  39.948326, 257.437716)
res4sl$sigma2<- c(38.30029, 58.10965, 85.06068, 19.16952)
res4sl$sigma3<- c(6.580927, 11.820816,  7.563079,  1.796753)
res4sl$sigma4<- c(65.40738, 79.61849, 41.91209, 14.24946)
res4sl$gamma <- matrix(c(0.884416338, 0.008210614, 0.10084565, 0.006527402,
                         0.006025303, 0.937861930, 0.04595484, 0.010157922,
                         0.080060102, 0.040454540, 0.77944006, 0.100045301,
                         0.018009185, 0.022190366, 0.16236337, 0.797437079),4,byrow=T)
res4sl$delta <- c(0.2471444, 0.2865394, 0.2971978, 0.1691183)

hres2sl <- list()
hres2sl$mu1<- c(40.81682, 113.60000)
hres2sl$mu2<- c(162.9102, 105.0562)
hres2sl$mu3<- c(21.836744,  4.603524)
hres2sl$mu4<- c(48.74350, 21.43112)
hres2sl$sigma1<- c(12.53681, 156.83692)
hres2sl$sigma2<- c(62.6588, 101.1413)
hres2sl$sigma3<- c(10.11698,  3.13330)
hres2sl$sigma4<- c(56.09231, 30.72411)
hres2sl$gamma <- matrix(c(0.96609502, 0.03390498, 0.06668338, 0.93331662),nrow=2,byrow=T)
hres2sl$delta <- c(0.6629333, 0.3370667)

hres3sl <- list()
hres3sl$mu1<- c(43.17531,  38.56498, 133.67557)
hres3sl$mu2<- c(155.98686, 166.67972,  91.54371)
hres3sl$mu3<- c(30.452695, 13.377403,  4.162385)
hres3sl$mu4<- c(45.79058, 50.33951, 14.99569)
hres3sl$sigma1<- c(12.50534,  12.72083, 194.91207)
hres3sl$sigma2<- c(43.73357, 80.95269, 91.70760)
hres3sl$sigma3<- c(4.794423, 6.432709, 2.920787)
hres3sl$sigma4<- c(49.66673, 60.30456, 21.24480)
hres3sl$gamma <- matrix(c(0.96173300, 0.01427472, 0.02399228, 0.01049263, 0.94707429,
                           0.04243308, 0.02583921, 0.06992465, 0.90423614),nrow=3,byrow=T)
hres3sl$delta <- c(0.2992625, 0.4336228, 0.2671147)

hres4sl <- list()
hres4sl$mu1<- c(38.24525,  37.90367, 402.57188,  42.40279)
hres4sl$mu2<- c(165.61557,  92.11664, 150.37865, 156.21455)
hres4sl$mu3<- c(14.262110,  3.215758, 11.306048, 30.649192)
hres4sl$mu4<- c(49.23663, 20.37289, 26.62026, 46.26177)
hres4sl$sigma1<- c(11.69191,  29.32017, 623.58415,  10.96548)
hres4sl$sigma2<- c(75.94096,  81.42040, 176.53191,  43.55601)
hres4sl$sigma3<- c(5.612117,  1.347204, 10.025119,  4.631337)
hres4sl$sigma4<- c(59.91370, 29.02331, 39.35957, 49.84760)

plotresults(N=2,res2sl,sealtype = 1, c(0.036,0.011,0.12,0.052))
plotresults(N=3,res3sl,sealtype = 1, c(0.047,0.0182,0.186,0.04))
plotresults(N=4,res4sl,sealtype = 1, c(0.062,0.025,0.25,0.04))
plotresults(N=2,hres2sl,sealtype = 0, c(0.036,0.009,0.17,0.052))
plotresults(N=3,hres3sl,sealtype = 0, c(0.033,0.0093,0.172,0.042))
plotresults(N=4,hres4sl,sealtype = 0, c(0.037,0.0095,0.315,0.042))

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

states <- viterbi(obs = obs_1,mod = res4sl,N=4)

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
        plot(obs[,2],type="l",col="blue",lwd=2, ylim= c(0,500), xlab = "Observation no.", ylab= "seconds",main="Surface Duration")
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

plot_viterbi(states[1:250],sealtype=1,obs=obs_1[1:250,],mod=res4sl,N=4)

#find delta with eigenvalues
P<- res3sl$gamma

eigen(t(P)) # here it is!
eigen_vect = eigen(t(P))$vectors[,1]
stat_dist = eigen_vect/sum(eigen_vect) # as there is subspace of them, 
# but we need the one with sum = 1
stat_dist
