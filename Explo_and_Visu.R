#2016_10_30_TW

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#------------------------------------------------------------------------------#
# important: run the User defined functions at the bottom of the script first! #
#                 (after loading the packages of course!)                      #
#------------------------------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# table of contents                       #        
#                                         #
# 1. load data and packages               #
# 2. prepare data                         #
# 3. exploratory data analysis            #
#   3.1 some general numbers              #
#   3.2 visualize movement                #
#   3.3 histograms of different variables #
#   3.4 correlation matrices              #
#   3.5 ACF's
# Appendix (UDF's)                        #

### 1. load data and packages
setwd("C:/Users/Pook/Documents/Robben")
seal_clean <- read.csv("seal_data_cleaned.csv")
require(moveHMM)
require(reshape2)
require(ggplot2)
require(grid)

### 2. prepare data

#remove unecessary columns for visualization and exploration
cols.dont.want <- c("X", "datetime","date", "time")
seals <- seal_clean[, ! names(seal_clean) %in% cols.dont.want, drop = F]

#rename lat and lon to y and x to use in moveHMM
colnames(seals)[colnames(seals) == "start.lon"] <- "x"
colnames(seals)[colnames(seals) == "start.lat"] <- "y"

#create data for every seal
seal1 <- seals[which(seals$sealID=="1"),]
seal2 <- seals[which(seals$sealID=="2"),]
seal3 <- seals[which(seals$sealID=="3"),]
seal4 <- seals[which(seals$sealID=="4"),]
seal5 <- seals[which(seals$sealID=="5"),]
seal6 <- seals[which(seals$sealID=="6"),]
seal7 <- seals[which(seals$sealID=="7"),]
seal8 <- seals[which(seals$sealID=="8"),]
seal9 <- seals[which(seals$sealID=="9"),]
seal10 <- seals[which(seals$sealID=="10"),]
seal11 <- seals[which(seals$sealID=="11"),]
seal12 <- seals[which(seals$sealID=="12"),]
seal13 <- seals[which(seals$sealID=="13"),]
seal14 <- seals[which(seals$sealID=="14"),]
seal15 <- seals[which(seals$sealID=="15"),]
seal16 <- seals[which(seals$sealID=="16"),]
seal17 <- seals[which(seals$sealID=="17"),]
seal18 <- seals[which(seals$sealID=="18"),]
seal19 <- seals[which(seals$sealID=="19"),]
seal20 <- seals[which(seals$sealID=="20"),]
seal21 <- seals[which(seals$sealID=="21"),]

### 3. exploratory data analysis

## 3.1 some general numbers:
#
# male: 16, female: 5
# grey: 10, harbour: 21

## 3.2 visualize movement of all seals(first 1000 datapoints ~200 days)
plotSat(data=seal1[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal2[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal3[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal4[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal5[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal6[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal7[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal8[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal9[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal10[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal11[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal12[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal13[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal14[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal15[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal16[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal17[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal18[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal19[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal20[1:1000,], zoom= 10, col ="yellow")
plotSat(data=seal21[1:1000,], zoom= 10, col ="yellow")

## 3.3 histograms of different variables
par(mfrow=c(3,7))
seallist <- list(seal1,seal2,seal3,seal4,seal5,seal6,seal7,seal8,seal9,seal10,seal11,seal12,seal13,seal14,seal15,seal16,seal17,seal18,seal19,seal20,seal21)
# for smaller monitors use
# par(mfrow=c(3,3))
# seallist <- list(seal1,seal2,seal3,seal4,seal5,seal6,seal7,seal8,seal9)

# surf duration
for (i in seallist){
        hist(i[,7],freq=F)
}
# not very useful lets try to do two for different scales
for (i in seallist){
        hist(i[,7],freq=F,breaks=200)
}
# lets look more into the area with values <500
for (i in seallist){
        hist(i[,7],freq=F,breaks=200,xlim=c(0,500))
}

#better! now for values >500
for (i in seallist){
        hist(i[,7],freq=F,breaks=200,xlim=c(500,8000),ylim=c(0,0.0001))
}

# same for dive duration (different scales!!)
for (i in seallist){
        hist(i[,8],freq=F)
}

# same for max dep (different scales!!)
for (i in seallist){
        hist(i[,9],freq=F)
}

# same for dist2coast (different scales!!)
for (i in seallist){
        hist(i[,11],freq=F)
}

## 3.4 correlation matrices

#correlation matrix for all seals
cormat <- round(cor(seals[,7:13]),2)
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
#create map and plot it
cormatplot<- ggheatmap(mc=melted_cormat)
print(cormatplot)
fin_corplot(ggheatmap=cormatplot)

##correlation matrix seal 1 -4
seal1_red <- seal1[,7:13]               #reduce matrix 
seal1_red <- seal1_red[,-4]             #in addition to above we need to remove sex bc we are looking at an individual
cormat1 <- round(cor(seal1_red),2)
cormat1 <- reorder_cormat(cormat1)
upper_tri1 <- get_upper_tri(cormat1)
melted_cormat1 <- melt(upper_tri1, na.rm = TRUE)
cormatplot1<- ggheatmap(mc=melted_cormat1)
s1_cor<-fin_corplot(ggheatmap=cormatplot1)

seal2_red <- seal2[,7:13]
seal2_red <- seal2_red[,-4]
cormat2 <- round(cor(seal2_red),2)
cormat2 <- reorder_cormat(cormat2)
upper_tri2 <- get_upper_tri(cormat2)
melted_cormat2 <- melt(upper_tri2, na.rm = TRUE)
cormatplot2<- ggheatmap(mc=melted_cormat2)
s2_cor<-fin_corplot(ggheatmap=cormatplot2)

seal3_red <- seal3[,7:13]
seal3_red <- seal3_red[,-4]
cormat3 <- round(cor(seal3_red),2)
cormat3 <- reorder_cormat(cormat3)
upper_tri3 <- get_upper_tri(cormat3)
melted_cormat3 <- melt(upper_tri3, na.rm = TRUE)
cormatplot3<- ggheatmap(mc=melted_cormat3)
s3_cor<-fin_corplot(ggheatmap=cormatplot3)

seal4_red <- seal4[,7:13]
seal4_red <- seal4_red[,-4]
cormat4 <- round(cor(seal4_red),2)
cormat4 <- reorder_cormat(cormat4)
upper_tri4 <- get_upper_tri(cormat4)
melted_cormat4 <- melt(upper_tri4, na.rm = TRUE)
cormatplot4<- ggheatmap(mc=melted_cormat4)
s4_cor<-fin_corplot(ggheatmap=cormatplot4)

multiplot(s1_cor,s2_cor,s3_cor,s4_cor,cols = 2)

## 3.5 ACF's

# max dep

nslist <- list(seal1$max.dep,seal2$max.dep,seal3$max.dep,seal4$max.dep,seal5$max.dep,seal6$max.dep,seal7$max.dep,seal8$max.dep,
               seal9$max.dep,seal10$max.dep,seal11$max.dep,seal12$max.dep,seal13$max.dep,seal14$max.dep,seal15$max.dep,
               seal16$max.dep,seal17$max.dep,seal18$max.dep,seal19$max.dep,seal20$max.dep,seal21$max.dep)
acfs<-lapply(nslist, FUN=ggacf)
multiplot(acfs[[1]],acfs[[2]],acfs[[3]],acfs[[4]],acfs[[5]],acfs[[6]],acfs[[7]],acfs[[8]],acfs[[9]],acfs[[10]],acfs[[11]],
          acfs[[12]],acfs[[13]],acfs[[14]],acfs[[15]],acfs[[16]],acfs[[17]],acfs[[18]],acfs[[19]],acfs[[20]],acfs[[21]],cols = 7)

# for smaller monitors
# multiplot(acfs[[1]],acfs[[2]],acfs[[3]],acfs[[4]],acfs[[5]],acfs[[6]],acfs[[7]],acfs[[8]],acfs[[9]],cols=3)

# dive.dur
nslist.ddur <- list(seal1[,8],seal2[,8],seal3[,8],seal4[,8],seal5[,8],seal6[,8],seal7[,8],seal8[,8],seal9[,8],seal10[,8],seal11[,8],
                    seal12[,8],seal13[,8],seal14[,8],seal15[,8],seal16[,8],seal17[,8],seal18[,8],seal19[,8],seal20[,8],seal21[,8])
acfsd<-lapply(nslist, FUN=ggacf)
multiplot(acfsd[[1]],acfsd[[2]],acfsd[[3]],acfsd[[4]],acfsd[[5]],acfsd[[6]],acfsd[[7]],acfsd[[8]],acfsd[[9]],acfsd[[10]],acfsd[[11]],
          acfsd[[12]],acfsd[[13]],acfsd[[14]],acfsd[[15]],acfsd[[16]],acfsd[[17]],acfsd[[18]],acfsd[[19]],acfsd[[20]],acfsd[[21]],cols = 7)
# for smaller monitors
# multiplot(acfsd[[1]],acfsd[[2]],acfsd[[3]],acfsd[[4]],acfsd[[5]],acfsd[[6]],acfsd[[7]],acfsd[[8]],acfsd[[9]],cols=3)

#--------------------------------------#
# self made functions, run these first #
#--------------------------------------#

get_upper_tri <- function(cormat){
        cormat[lower.tri(cormat)]<- NA
        return(cormat)
}

reorder_cormat <- function(cormat){
        # Use correlation between variables as distance
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <-cormat[hc$order, hc$order]
}

#function for multiplotting, since ggplot does not support base multiplotting
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        require(grid)
        
        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)
        
        numPlots = length(plots)
        
        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }
        
        if (numPlots==1) {
                print(plots[[1]])
                
        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                        
                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}

# Create a ggheatmap
ggheatmap <- function(mc){
        ggplot(mc, aes(Var2, Var1, fill = value))+
                geom_tile(color = "white")+
                scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                     midpoint = 0, limit = c(-1,1), space = "Lab", 
                                     name="Pearson\nCorrelation") +
                theme_minimal()+ # minimal theme
                theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                                 size = 12, hjust = 1))+
                coord_fixed()
        
}

# add correlations and make plot nicer
fin_corplot<-function(ggheatmap){
        ggheatmap + 
                geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
                theme(
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.ticks = element_blank(),
                        legend.justification = c(1, 0),
                        legend.position = "right",
                        legend.direction = "vertical")+
                guides(fill = guide_colorbar(barwidth = 3, barheight = 15,
                                             title.position = "top", title.hjust = 0.5))
}

ggacf <- function(dvec){
        bacf <- acf(dvec, plot = FALSE)
        bacfdf <- with(bacf, data.frame(lag, acf))
        confint1 <- 1.96/sqrt(length(dvec))
        c <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf))
        c <- c + geom_hline(aes(yintercept = 0)) 
        c <- c + geom_segment(mapping = aes(xend = lag, yend = 0))   
        c <- c + geom_hline(aes(yintercept = confint1), color="blue", size = 1, linetype = "F1")
        c <- c + geom_hline(aes(yintercept = -confint1), color="blue", size = 1, linetype = "F1")
        c <- c + theme(text = element_text(size=17))
        return(c)
}