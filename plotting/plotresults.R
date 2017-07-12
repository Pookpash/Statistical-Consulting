plotresults <- function(res,N, colvec,bildvec){ #bildvec hat 4 Einträge, jeweils die rechten grenzen im Plot, linke Grenzen und Turning angle sind klar
  par(mfrow=c(3,2),oma=c(0,0,2,0))
  x1 <- seq(1,bildvec[1],length=1000)
  x2 <- seq(1,bildvec[2],length=1000)
  x3 <- seq(0.1,bildvec[3],length=1000)
  x4 <- seq(0.1,bildvec[4],length=1000)
  x5 <- seq(-pi,pi,length=1000)
  
  m1 <- rep(0,N)
  for(i in 1:N){
    m1[i]<-max(dgamma(x1,res$mu1[i]^2/res$sigma1[i]^2,scale=res$sigma1[i]^2/res$mu1[i]))
  }
  plot(x1,dgamma(x1,res$mu1[1]^2/res$sigma1[1]^2,scale=res$sigma1[1]^2/res$mu1[1]), type="l", col= colvec[1], ylim=c(0,max(m1)), main="Maximum Depth", ylab="density", xlab= "meter", lwd=2)
  for(i in 2:N){
    points(x1,dgamma(x1,res$mu1[i]^2/res$sigma1[i]^2,scale=res$sigma1[i]^2/res$mu1[i]), type="l", col=colvec[i], lwd=2)
  }
  
  m2 <- rep(0,N)
  for(i in 1:N){
    m2[i]<-max(dgamma(x2,res$mu2[i]^2/res$sigma2[i]^2,scale=res$sigma2[i]^2/res$mu2[i]))
  }
  plot(x2,dgamma(x2,res$mu2[1]^2/res$sigma2[1]^2,scale=res$sigma2[1]^2/res$mu2[1]), type="l", col=colvec[1], ylim=c(0,max(m2)), main="Dive Duration", ylab="density", xlab= "seconds", lwd=2)
  for(i in 2:N){
    points(x2,dgamma(x2,res$mu2[i]^2/res$sigma2[i]^2,scale=res$sigma2[i]^2/res$mu2[i]), type="l", col=colvec[i], lwd=2)
  }
  
  m3 <- rep(0,N)
  for(i in 1:N){
    m3[i]<-max(dgamma(x3,res$mu3[i]^2/res$sigma3[i]^2,scale=res$sigma3[i]^2/res$mu3[i]))
  }
  plot(x3,dgamma(x3,res$mu3[1]^2/res$sigma3[1]^2,scale=res$sigma3[1]^2/res$mu3[1]), type="l", col=colvec[1], ylim=c(0,max(m3)), main="Surface Duration", ylab="density", xlab= "seconds", lwd=2)
  for(i in 2:N){
    points(x3,dgamma(x3,res$mu3[i]^2/res$sigma3[i]^2,scale=res$sigma3[i]^2/res$mu3[i]), type="l", col=colvec[i], lwd=2)
  }
  
  m4 <- rep(0,N)
  for(i in 1:N){
    m4[i]<-max(dgamma(x4,res$mu4[i]^2/res$sigma4[i]^2,scale=res$sigma4[i]^2/res$mu4[i]))
  }
  plot(x4,dgamma(x4,res$mu4[1]^2/res$sigma4[1]^2,scale=res$sigma4[1]^2/res$mu4[1]), type="l", col=colvec[1], ylim=c(0,max(m4)), main="Steplength", ylab="density", xlab= "meter", lwd=2)
  for(i in 2:N){
    points(x4,dgamma(x4,res$mu4[i]^2/res$sigma4[i]^2,scale=res$sigma4[i]^2/res$mu4[i]), type="l", col=colvec[i], lwd=2)
  }
  
  m5 <- rep(0,N)
  for(i in 1:N){
    m5[i]<-max(dvm(x5,res$mu[i],res$kappa[i]))
  }
  plot(x5,dvm(x5,res$mu[1],res$kappa[1]), type="l", col=colvec[1], ylim=c(0,max(m5)), main="turning angle", ylab="density", xlab= "radient", lwd=2)
  for(i in 2:N){
    points(x5,dvm(x5,res$mu[i],res$kappa[i]), type="l", col=colvec[i], lwd=2)
  }
  mtext(title1, outer = TRUE, cex = 1.5)
  par(mfrow=c(1,1))
}
