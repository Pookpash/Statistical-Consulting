
data2 <- transform(data, datetime = colsplit(datetime, split = "\\ ", names = c("date", "time"))) 
#this produces a nested dataframe which can cause problems. therefore saved seperately and the new columns are added to the old df.
#I kept the original datetime column because we might sill need it.
date <- data2$datetime$date
time <- data2$datetime$time
data$date <- date
data$time <- time
