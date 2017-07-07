## Funktionen Plots

resids_hist <- function(resids,var=5){
  par(mfrow=c(2,3))
  title <- c("Surface Duration","Dive Duration","Maximum Depth","Steplength","Turning Angle")
  for(i in 1:var){
    hist(resids[[i]],prob=T,xlab=NULL,main=title[i],col="grey")
    xfit <- seq(min(na.omit(resids[[i]])),max(na.omit(resids[[i]])),length=100)
    yfit <- dnorm(xfit)
    lines(xfit, yfit, lwd=2, col=2)
  }
}

resids_qq <- function(resids,var=5){
  par(mfrow=c(2,3))
  title <- c("Surface Duration","Dive Duration","Maximum Depth","Steplength","Turning Angle")
  for(i in 1:var){
    res <- sample(unique(resids[[i]]),size=50000)
    qqnorm(res,main=title[i])
    qqline(res,col=2,lwd=3)
  }
}

resids_acf <- function(resids,var=5){
  par(mfrow=c(2,3))
  title <- c("Surface Duration","Dive Duration","Maximum Depth","Steplength","Turning Angle")
  for(i in 1:var){
    acf(na.omit(resids[[i]]),main=title[i])
  }
}
