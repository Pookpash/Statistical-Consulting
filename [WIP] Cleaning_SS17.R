### load packages and set wd
require(maptools)

setwd("C:/Users/Pook/Documents/Robben")
data <- read.delim("~/Robben/seal_data_20170420.txt")

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

#wip :turning angles

#turning angles
trackAngle <- function(xy) {
  angles <- abs(c(trackAzimuth(xy), 0) -
                  c(0, rev(trackAzimuth(xy[nrow(xy):1, ]))))
  angles[is.na(angles)] <- 180
  angles[-c(1, length(angles))]
}

180 - trackAngle(as.matrix(data[1:1000,c(4,5)]))

data$turn.angle <- angles

#write to new csv-file to prevent repeating the steps above everytime
#one cleaned file with almost all columns
#one working file with the columns i think we will use
#one working file with only the grey seals
write.csv(data, file = "seal_data_cleaned.csv")

datax <- data[,c(1,2,8,10:13,15:20)]


write.csv(datax, file ="seal_work.csv")

datag <- datax[1:91830,] #needs to be adjusted after turning angles have been implemented


write.csv(datag, file ="greyseal_work.csv")

