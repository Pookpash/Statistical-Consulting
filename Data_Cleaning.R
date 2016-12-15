#2016_12_15_TW

### load packages and set wd
setwd("C:/Users/Admin/Documents/Robben") #change as needed
require(reshape)

### load original data
data <- read.delim("~/Robben/seal_data.txt")

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

# following code looks nicer but is too slow
#zerolist <- c("01","02","03","04","05","06","07","08","09")
#for (j in 1:(dim(data)[1]-500000)){ #remove "0"s where not necessary (-500000 to speed things up, if u want to check do "data[dim(data)[1]-500000,1]"
#        for (i in zerolist){
#                if (data[j,1]==i){
#                        data[j,1] <- gsub("0","",data[j,1])
#                }
#        }
#}

#second column
data[,2] <- gsub("Grey seal","1",data[,2])
data[,2] <- gsub("Harbour seal","0",data[,2])
data[,2] <- as.integer(data[,2])

#3rd column
data2 <- transform(data, datetime = colsplit(datetime, split = "\\ ", names = c("date", "time"))) 
#this produces a nested dataframe which can cause problems. therefore saved seperately and the new columns are added to the old df.
#I kept the original datetime column because we might sill need it.
date <- data2$datetime$date
time <- data2$datetime$time
data$date <- date
data$time <- time

#11th column                                                                                                 
data[,11] <- gsub("Male","0",data[,11])
data[,11] <- gsub("Female","1",data[,11])
data[,11] <- as.integer(data[,11])

#reorder dataframe
data <- data[,c(1,2,3,14,15,4,5,6,7,8,9,10,11,12,13)]

# Calculate distance in kilometers between two points
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
steplen <-earth.dist(long1=data[,7],lat1=data[,6],long2=data[,9],lat2=data[,8])
steplen <- steplen *1000 # convert to meters
data[,16] <- steplen
colnames(data)[colnames(data) == "V16"] <- "steplen"

data[,15] <- round(data[,15],2)

### create csv
#write to new csv-file to prevent repeating the steps above everytime
write.csv(data, file = "seal_data_cleaned.csv")

#test new csv-file
seal_clean <- read.csv("seal_data_cleaned.csv")
View(seal_clean) #looks good, "X" column can be dropped in further analysis(discussion)
