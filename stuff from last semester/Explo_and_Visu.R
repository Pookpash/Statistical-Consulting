
### 1. load data and packages
setwd("C:/Users/Pook/Documents/Robben")
seal_clean <- read.csv("seal_data_cleaned.csv")
require(moveHMM)
require(reshape2)
require(ggplot2)
require(grid)

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
