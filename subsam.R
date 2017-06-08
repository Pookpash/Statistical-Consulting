setwd("E:\\Master\\Robbenkloppen")
data <- read.csv("seal_data_cleaned.csv")

colnames(data)[5] <- "x" #some renaming could be moved to data cleaning
colnames(data)[6] <- "y"


subsam <- 20 #subsample size
numvar <- 5 #number of variables that will be included in the data set produced




#turning angle function (from Roland)

turnAngle <- function(x,y,z){
  v <- c(y[1]-x[1],y[2]-x[2])
  w <- c(z[1]-y[1],z[2]-y[2])
  angle <- atan2(w[2],w[1])-atan2(v[2],v[1])
  if (angle <= -pi) angle <- angle + 2*pi
  if (angle > pi) angle <- angle -2*pi
  return(angle)
}



#function that has one ID_burst as input and handles the calculations for it

handleburst <- function(oneburstdata,subsam){
  numleftover <- nrow(oneburstdata) %% subsam
  leftover <- oneburstdata[(nrow(oneburstdata)-numleftover+1):nrow(oneburstdata),]
  eaten <- oneburstdata[1:(nrow(oneburstdata)-numleftover),]
  nss <- nrow(eaten)/subsam
  out <- matrix(rep(NA,(nss+1)*numvar), nrow = nss+1, ncol = numvar) #adjust according to no of variables (possible max is ncol(eaten))
  for (i in 1:nss){
    out[i,1] <-mean(eaten$dive_dep[((i-1)*subsam+1):(i*subsam)])
    out[i,2] <-mean(eaten$dive.dur[((i-1)*subsam+1):(i*subsam)])
    out[i,3] <-mean(eaten$surf.dur[((i-1)*subsam+1):(i*subsam)])
  }
 
  out[i+1,1] <- mean(leftover$dive_dep)
  out[i+1,2] <- mean(leftover$dive.dur)
  out[i+1,3] <- mean(leftover$surf.dur)

eatenss <- eaten[seq(1,nrow(eaten),by=subsam),]
steps <- turns <- rep(NA,nss)
  for (k in 2:(nss-1)){
    x<-c(eatenss$x[k-1],eatenss$y[k-1])
    y<-c(eatenss$x[k],eatenss$y[k])
    z<-c(eatenss$x[k+1],eatenss$y[k+1])
    if (!is.na(x[1]) & !is.na(y[1]) & !is.na(z[1])) turns[k]<-turnAngle(x,y,z)
    steps[k]<-sqrt(sum((z-y)^2))
  }  
 
out[1:nrow(eatenss),4] <- turns
out[1:nrow(eatenss),5] <- steps
 
    return(out)
}

split2list <- split(data, f = data$ID_burst) #split ID_bursts into a list of matricies
subsamplelist <- lapply(split2list, function(x) handleburst(x,subsam)) #apply function to each ID_burst

subsample <- do.call(rbind,subsamplelist) #reunite into final dataset


write.csv(subsample, file = "subsample_20.csv")




