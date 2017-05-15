

setwd("C:/Users/Admin/Documents/Robben")
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

#turning angles

yberechnen <- function(w,x,y,z){
        if (x[1]-w[1]==0){ #prevent dividing by zero
                helpcalc1 <- 0.00000000000001
        } else {
                helpcalc1 <- x[1]-w[1]
        }
        
        m <- (x[2]-w[2])/helpcalc1
        b <- x[2]-m*x[1]
        
        if (z[1]-y[1]==0){ #prevent dividing by zero
                helpcalc2 <- 0.00000000000001
        } else {
                helpcalc2 <- z[1]-y[1]
        }
        mm <- (z[2]-y[2])/helpcalc2
        bb <- z[2]-y*x[1]
        if(mm==m){ #prevent problem induced above
                m<-m+0.000000000000000000000000001
        }
        yy <- c((b-bb)/(mm-m))
        yy[2] <- m*yy[1]+b
        return(yy)
}

#w startpunkt vorheriger
#x endpunkt vorher
#y start aktuell
#z end aktuell

y_intervec <- matrix(c(rep(0, 2*(length(data[,1])+1))),nrow=length(data[,1])+1)
        #correct size of vector specified beforehand 
        #to make code quicker since R does not 
        #allocate memory very well. Speed difference (~30times faster!)

for (i in 2: length(data[,1])){
        w <- c(data[i-1,4], data[i-1,5])
        x <- c(data[i-1,6], data[i-1,7])
        y <- c(data[i,4], data[i,5])
        z <- c(data[i,6], data[i,7])
        y_int <- yberechnen(w,x,y,z)
        y_intervec[i,1] <- y_int[1]
        y_intervec[i,2] <- y_int[2]
}

y_intervec <- as.matrix(y_intervec,ncol=2,byrow=T)

turnAngle <- function(x,y,z){
        v <- c(y[1]-x[1],y[2]-x[2])
        w <- c(z[1]-y[1],z[2]-y[2])
        angle <- atan2(w[2],w[1])-atan2(v[2],v[1])
        while(angle <= (-pi)){
                angle <- angle + 2*pi
        }
        while(angle > pi){
                angle <- angle -2*pi
        }
        return(angle)
}

turnang <- c(rep(0,length(data[,1])+1)) #correct size of vector specified beforehand 
                                        #to make code quicker since R does not 
                                        #allocate memory very well.
for (i in 2:length(data[,1])){
        x <- c(data[i-1,6], data[i-1,7])
        y <- y_intervec[i,]
        z <- c(data[i,4], data[i,5])
        turnang[i] <- turnAngle(x,y,z)
}

#put in on the dataframe

data$angle <- turnang

#get rid of wrong angles
ID_burst_vector <- match(unique(data$ID_burst), data$ID_burst)
data <- data[-ID_burst_vector,]

#x endpunkt vorheriger tauchgang, vec mit 2 elem
#y vorher berechnet
#z startpunkt aktueller tauchgang, vec 2 elem

#write to new csv-file to prevent repeating the steps above everytime
#one cleaned file with almost all columns
#one working file with the columns i think we will use
#one working file with only the grey seals
write.csv(data, file = "seal_data_cleaned.csv")

datax <- data[,c(1,2,8,10:13,15:20)]


write.csv(datax, file ="seal_work.csv")

grey_end <- match(unique(data$ID), data$ID)[12]
datag <- datax[1:(grey_end-1),] #needs to be adjusted after turning angles have been implemented


write.csv(datag, file ="greyseal_work.csv")


datah <- datax[grey_end:length(datax[,1]),]
write.csv(datag, file="harbseal_work.csv")
