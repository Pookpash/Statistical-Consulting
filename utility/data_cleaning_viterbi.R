setwd("E:\\Master\\Robbenkloppen")
data <- read.delim("seal_data_20170420.txt")
library(dplyr)

#remove unneeded columns
data <- data[,c(1,2,4:18)]

### datacleaning
#first column
data[,1] <- gsub("S01","12",data[,1])
data[,1] <- gsub("S02","13",data[,1])
data[,1] <- gsub("S03","14",data[,1])
data[,1] <- gsub("S04","15",data[,1])
data[,1] <- gsub("S05","16",data[,1])
data[,1] <- gsub("S06","17",data[,1])
data[,1] <- gsub("S07","18",data[,1])
data[,1] <- gsub("S08","19",data[,1])
data[,1] <- gsub("S09","20",data[,1])
data[,1] <- gsub("S10","21",data[,1])
data[,1] <- gsub("G", "",data[,1])   #remove "G"s
data[,1] <- gsub("01","1", data[,1]) #remove unnecessary "0"s
data[,1] <- gsub("02","2", data[,1])
data[,1] <- gsub("03","3", data[,1])
data[,1] <- gsub("04","4", data[,1])
data[,1] <- gsub("05","5", data[,1])
data[,1] <- gsub("06","6", data[,1])
data[,1] <- gsub("07","7", data[,1])
data[,1] <- gsub("08","8", data[,1])
data[,1] <- gsub("09","9", data[,1])
data[,1] <- as.factor(data[,1])
#Seal 1 to 11 are the grey seals

#sex column
data[,11] <- gsub("Male","0",data[,11])
data[,11] <- gsub("Female","1",data[,11])
data[,11] <- as.integer(data[,11])
#male = 0

#steplen function (in km)
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}
steplen <-earth.dist(long1=data[,4],lat1=data[,5],long2=data[,6],lat2=data[,7])
steplen <- steplen *1000 # convert to meters
data$steplen <- steplen
colnames(data)[colnames(data) == "V16"] <- "steplen"

#calculate % of bathymetry
dive_dep2 <- data$max.dep/data$bathymetry
dive_dep2[dive_dep2>=1.0] <- 0.999999 #to prevent values > 1
data$dive_dep <- dive_dep2

#dive id
levels(data$ID_burst)<-c(1:1020) 
#this may lead to problems when analyzing them by species
#in that case, just rerun the code above for the Harbour Seals dataset

data <- data[,-c(3,6:7,14,19)]

grey_end <- match(unique(data$ID), data$ID)[12]
data <- data[1:(grey_end-1),]

data$ID_burst <- as.numeric(data$ID_burst)
data$ID <- as.numeric(data$ID)

#subsample

colnames(data)[3] <- "x" #some renaming could be moved to data cleaning
colnames(data)[4] <- "y"


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
  if (numleftover > 0) {
    leftover <- oneburstdata[(nrow(oneburstdata)-numleftover+1):nrow(oneburstdata),]
  }
  eaten <- oneburstdata[1:(nrow(oneburstdata)-numleftover),]
  nss <- nrow(eaten)/subsam
  if (numleftover > 0) {
    out <- matrix(rep(NA,(nss+1)*(numvar+2)), nrow = nss+1, ncol = numvar+2) #adjust according to no of variables (possible max is ncol(eaten))
  } else {
    out <- matrix(rep(NA,(nss)*(numvar+2)), nrow = nss, ncol = numvar+2)
  }
  for (i in 1:nss){
    out[i,4] <- mean(eaten$max.dep[((i-1)*subsam+1):(i*subsam)])
    out[i,5] <- mean(eaten$dive.dur[((i-1)*subsam+1):(i*subsam)])
    out[i,6] <- mean(eaten$surf.dur[((i-1)*subsam+1):(i*subsam)])
    out[i,7] <- mean(eaten$dist2coast[((i-1)*subsam+1):(i*subsam)])
    out[i,8] <- mean(eaten$dist2haulout[((i-1)*subsam+1):(i*subsam)])
    out[i,9] <- mean(eaten$dist2fishingnet[((i-1)*subsam+1):(i*subsam)])
    out[i,10] <- mean(eaten$sea_surface_temperature[((i-1)*subsam+1):(i*subsam)])
    out[i,11] <- mean(eaten$sea_surface_salinity[((i-1)*subsam+1):(i*subsam)])
    out[i,12] <- mean(eaten$steplen[((i-1)*subsam+1):(i*subsam)])
    out[i,2] <- eaten$ID_burst[1]
    out[i,3] <- eaten$sex[1]
    out[i,1] <- eaten$ID[1]
    out[i,13] <- eaten$x[i*subsam]
    out[i,14] <- eaten$y[i*subsam]
    
  }
  if (numleftover > 0){
    
    out[i+1,4] <- mean(leftover$max.dep)
    out[i+1,5] <- mean(leftover$dive.dur)
    out[i+1,6] <- mean(leftover$surf.dur)
    out[i+1,7] <- mean(leftover$dist2coast)
    out[i+1,8] <- mean(leftover$dist2haulout)
    out[i+1,9] <- mean(leftover$dist2fishingnet)
    out[i+1,10] <- mean(leftover$sea_surface_temperature)
    out[i+1,11] <- mean(leftover$sea_surface_salinity)
    out[i+1,12] <- mean(leftover$steplen)
    out[i+1,2] <- eaten$ID_burst[1]
    out[i+1,3] <- eaten$sex[1]
    out[i+1,1] <- eaten$ID[1]
    out[i+1,13] <- leftover$x[length(leftover$x)]
    out[i+1,14] <- leftover$y[length(leftover$y)]
  }
  
  eatenss <- eaten[seq(1,nrow(eaten),by=subsam),]
  turns <- rep(NA,nss)
  for (k in 2:(nss)){
    x<-c(eatenss$x[k-1],eatenss$y[k-1])
    y<-c(eatenss$x[k],eatenss$y[k])
    z<-c(eatenss$x[k+1],eatenss$y[k+1])
    if (!is.na(x[1]) & !is.na(y[1]) & !is.na(z[1])) turns[k]<-turnAngle(x,y,z)
  }  
  
  out[1:nrow(eatenss),15] <- turns
  
  return(out)
}



split2list <- split(data, f = as.factor(data$ID_burst)) #split ID_bursts into a list of matricies

subsam <- 10 #subsample size
numvar <- 13 #number of variables that will be included in the data set produced


subsamplelist <- lapply(split2list, function(x) handleburst(x,subsam)) #apply function to each ID_burst

subsample <- do.call(rbind,subsamplelist) #reunite into final dataset
colnames(subsample) <- c("ID","ID_burst","sex","max.dep","dive.dur","surf.dur","dist2coast","dist2haulout","dist2fishingnet","temperature","salinty","steplen","x","y","turna")
subsample <- as.data.frame(subsample)
subsample <- subsample[complete.cases(subsample),]

subsample$ID <- factor(subsample$ID, c(1,12,15:21,2:3), labels = c(1:11))

write.csv(subsample, file = paste0("subsample_vit",subsam,".csv"))
