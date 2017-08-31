seal_ns<- read.csv("subsample_10.csv")

df <- data.frame(rep(NA,5),rep(NA,5))
colnames(df) <- c("mean","sd")
rownames(df) <- colnames(seal_ns)[8:12]
for(i in 1:5){
  df[i,1] <- mean(seal_ns[,i+7])
}
for(i in 1:5){
  df[i,2] <- sd(seal_ns[,i+7])
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



plotstationarymitcov <- function(mod, k, df= df, col=cbPalette[2:4]){  
  beta0 <- mod$beta[1,]
  betasex <- mod$beta[4,]
  betak <- mod$beta[k,]
  if(k == 3){
    covs <- "dist2coast"
    xunscale <- seq(1800,117000, length = 2000)
    x <- (xunscale - df[covs,"mean"])/(df[covs, "sd"])
    achse <- seq(0,120000, by=20000)
    achse <- (achse - df[covs,"mean"])/df[covs,"sd"]
  }
  if(k == 4){
    covs <- "dist2haulout"
    xunscale <- seq(1800,290000, length = 2000)
    x <- (xunscale - df[covs,"mean"])/(df[covs, "sd"])
    achse <- seq(1800,290000, by=20000)
    achse <- (achse - df[covs,"mean"])/df[covs,"sd"]
  }
  if(k == 5){
    covs <- "dist2fishingnet"
    xunscale <- seq(300,510000, length = 2000)
    x <- (xunscale - df[covs,"mean"])/(df[covs, "sd"])
    achse <- seq(300,510000, by=50000)
    achse <- (achse - df[covs,"mean"])/df[covs,"sd"]
  }
  if(k == 6){
    covs <- "temperature"
    xunscale <- seq(-5,20, length = 2000)
    x <- (xunscale - df[covs,"mean"])/(df[covs, "sd"])
    achse <- seq(-5,20, by=5)
    achse <- (achse - df[covs,"mean"])/df[covs,"sd"]
  }
  if(k == 7){
    covs <- "salinty"
    xunscale <- seq(0,27, length = 2000)
    x <- (xunscale - df[covs,"mean"])/(df[covs, "sd"])
    achse <- seq(0,27, by=5)
    achse <- (achse - df[covs,"mean"])/df[covs,"sd"]
  }
  male <- matrix(0,3,length(x))
  for(i in 1:length(x)){
    gammamale <- gammafemale <- matrix(0,3,3)
    gammamale[1,2] <- exp(beta0[1]+betak[1]*x[i])/(1+exp(beta0[1]+betak[1]*x[i])+exp(beta0[2]+betak[2]*x[i]))
    gammamale[1,3] <- exp(beta0[2]+betak[2]*x[i])/(1+exp(beta0[1]+betak[1]*x[i])+exp(beta0[2]+betak[2]*x[i]))
    gammamale[1,1] <-1/(1+exp(beta0[1]+betak[1]*x[i])+exp(beta0[2]+betak[2]*x[i]))
    
    gammamale[2,1] <- exp(beta0[3]+betak[3]*x[i])/(1+exp(beta0[3]+betak[3]*x[i])+exp(beta0[4]+betak[4]*x[i]))
    gammamale[2,3] <- exp(beta0[4]+betak[4]*x[i])/(1+exp(beta0[3]+betak[3]*x[i])+exp(beta0[4]+betak[4]*x[i]))
    gammamale[2,2] <- 1/(1+exp(beta0[3]+betak[3]*x[i])+exp(beta0[4]+betak[4]*x[i]))
    
    gammamale[3,1] <- exp(beta0[5]+betak[5]*x[i])/(1+exp(beta0[5]+betak[5]*x[i])+exp(beta0[6]+betak[6]*x[i]))
    gammamale[3,2] <- exp(beta0[6]+betak[6]*x[i])/(1+exp(beta0[5]+betak[5]*x[i])+exp(beta0[6]+betak[6]*x[i]))
    gammamale[3,3] <- 1/(1+exp(beta0[5]+betak[5]*x[i])+exp(beta0[6]+betak[6]*x[i]))
    
    deltamale <- solve( t(diag(3)-gammamale+matrix(1,3,3)), c(1,1,1) )
    male[,i] <- deltamale
  }
  
  female <- matrix(0,3,length(x))
  for(i in 1:length(x)){
    gammafemale <- gammafemale <- matrix(0,3,3)
    gammafemale[1,2] <- exp(beta0[1]+betasex[1]+betak[1]*x[i])/(1+exp(beta0[1]+betasex[1]+betak[1]*x[i])+exp(beta0[2]+betasex[2]+betak[2]*x[i]))
    gammafemale[1,3] <- exp(beta0[2]+betasex[2]+betak[2]*x[i])/(1+exp(beta0[1]+betasex[1]+betak[1]*x[i])+exp(beta0[2]+betasex[2]+betak[2]*x[i]))
    gammafemale[1,1] <-1/(1+exp(beta0[1]+betasex[1]+betak[1]*x[i])+exp(beta0[2]+betasex[2]+betak[2]*x[i]))
    
    gammafemale[2,1] <- exp(beta0[3]+betasex[3]+betak[3]*x[i])/(1+exp(beta0[3]+betasex[3]+betak[3]*x[i])+exp(beta0[4]+betasex[4]+betak[4]*x[i]))
    gammafemale[2,3] <- exp(beta0[4]+betasex[4]+betak[4]*x[i])/(1+exp(beta0[3]+betasex[3]+betak[3]*x[i])+exp(beta0[4]+betasex[4]+betak[4]*x[i]))
    gammafemale[2,2] <- 1/(1+exp(beta0[3]+betasex[3]+betak[3]*x[i])+exp(beta0[4]+betasex[4]+betak[4]*x[i]))
    
    gammafemale[3,1] <- exp(beta0[5]+betasex[5]+betak[5]*x[i])/(1+exp(beta0[5]+betasex[5]+betak[5]*x[i])+exp(beta0[6]+betasex[6]+betak[6]*x[i]))
    gammafemale[3,2] <- exp(beta0[6]+betasex[6]+betak[6]*x[i])/(1+exp(beta0[5]+betasex[5]+betak[5]*x[i])+exp(beta0[6]+betasex[6]+betak[6]*x[i]))
    gammafemale[3,3] <- 1/(1+exp(beta0[5]+betasex[5]+betak[5]*x[i])+exp(beta0[6]+betasex[6]+betak[6]*x[i]))
    
    deltafemale <- solve( t(diag(3)-gammafemale+matrix(1,3,3)), c(1,1,1) )
    female[,i] <- deltafemale
  }
  ma <- list(male,female)
  pdf("blubber.pdf",width = 6 , height= 7)
  par(mfrow=c(2,1))
  plot(x,ma[[1]][1,],typ="l",ylim = c(0,1),ylab = "probability",xaxt="n", col=col[1],xlab= covs, main = "male")
  axis(1,achse,round(achse*df[covs,"sd"]+df[covs,"mean"])/1000)
  lines(x,ma[[1]][2,],col=col[2])
  lines(x,ma[[1]][3,],col=col[3])
  
  plot(x,ma[[2]][1,],typ="l",ylim = c(0,1),ylab = "probability",xaxt="n", col=col[1],xlab= covs, main = "female")
  axis(1,achse,round(achse*df[covs,"sd"]+df[covs,"mean"]))
  lines(x,ma[[2]][2,],col=col[2])
  lines(x,ma[[2]][3,],col=col[3])
  graphics.off()
  return(list(male,female))
}  

plotstationarymitcov(modwithall[[1]],3, df= df, col=cbPalette[2:4])
