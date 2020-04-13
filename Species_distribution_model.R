
######################
# Charge libraries
######################
library(sp)
library(car)
library(nlme)
library(dismo)
library(gstat)
library(spdep)
library(Hmisc)
library(raster)
library(leaflet)
library(GGally)
library(ggplot2)
library(maptools)
library(corrplot)
library(sdmpredictors)
library(PresenceAbsence)

#####################################
### --- Read the data         --- ###
#####################################
data <- readRDS("data.rds")
str(data)

############################################
#  Check first 6 rows of the dataset      #
############################################
head(data)

############################################
#         How many rows and colums?        #
############################################
dim(data)

##########################NOTE##################################
#Coordinates are typically expressed as longitude and latitude,#
#but they could also be Easting and Northing in UTM or         #
#another planar coordinate reference system (map projection).  #
#Check always the format before analysis the data.             #
#################################################################

##############################################################
#Select the records that have longitude and latitude data
#############################################################
colnames(data)

###################################################################
#It is important to make maps to assure that the points are, at least
#roughly, in the right location.
##################################################################
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-80,70), ylim=c(-60,10), axes=TRUE,col="light yellow")

##############################################################
#Plot points
##############################################################
points(data$Lon, data$Lat, col="orange", pch=20, cex=0.75)
points(data$Lon, data$Lat, col="red", cex=0.75)

################################################################################
#Argentine hake is a species that occurs in the Southwest                     #
#Atlantic (off southern Brazil to Argentina and the Malvinas Islands).        #
#Do you see any errors on the map?                                            #
#There are a few records that map in the Indian Ocean. Any idea why that      #
#may have happened? It is a common mistake, missing minus signs.              #
#The coordinates are around (50.78, -31.06) but they should in Brazil         #
#around (-50.78, -31.06).                                                     #
################################################################################

click()#get coordinates
identify(data$Lon, data$Lat)#get position in the dataset

################################################################################
#So the latitude is probably correct, and erroneously copied to the longitude.
#In this case we eliminate these observations from the dataset.
################################################################################
data<-data[-c(259,260), ]

################################################################################
#Plot again to see if we eliminated them
################################################################################
plot(wrld_simpl, xlim=c(-80,70), ylim=c(-60,10), axes=TRUE,col="light yellow")
points(data$Lon, data$Lat, col="orange", pch=20, cex=0.75)
points(data$Lon, data$Lat, col="red", cex=0.75)

################################################################################
#which records are duplicates?
################################################################################
dups2 <- duplicated(data[, c("Species","Lon","Lat")]);dups2

################################################################################
#Number of duplicates
################################################################################
sum(dups2)

################################################################################
#Remove duplicates
################################################################################
data <- data[!dups2, ]

#################################################################################
#                               Environmental variables          ##
##############################################################################
#Get the file names and list layer codes for Bio-ORACLE and MARSPEC
list_layers(c("Bio-ORACLE","MARSPEC"))$layer_code

predictors<- load_layers(c("BO_calcite","BO_chlomean","BO_nitrate","BO_ph","BO_phosphate",
                           "BO_silicate","BO_salinity","BO_sstmean"))# datadir = tempdir())

names(predictors) <- c("calcite","chlomean","nitrate","ph","phos","silicate","salinity","sstmean")
plot(predictors)
predictors#check extent, resolution and projection

################################################################################
#                      Crop for your areas
################################################################################
ext<-extent(-70,-30,-55,-20)
predictors<-crop(predictors,ext)
plot(predictors)

################################################################################
#We can also make a plot of a single layer in a RasterStack, and plot data on top of it.
################################################################################

#tiff("Figure.tiff",3200,3500, res=300)#this is to save as a tiff file
par(mfrow=c(1,1))
plot(predictors, 8)
plot(wrld_simpl, col="grey",add=TRUE)
points(data$Lon, data$Lat, col='blue', pch=16)
#dev.off()

######################################################################################################################
###             Standardize predictors 
######################################################################################################################
predictors2<-scale(predictors)
round(apply(values(predictors2), 2, summary), 4)

################################################################################
#    How get variables values in presence records?
################################################################################
coords1=cbind(data$Lon, data$Lat)
colnames(coords1)<-c("x","y")
presvals <- extract(predictors2, coords1)
head(presvals)

##################################################################
#                   Generate Pseudo-absence
####################################################################
#Setting random seed to always create the same random set of points for this example
set.seed(500)
backgr <- randomPoints(predictors2, 200)
backgr <-as.data.frame(backgr)
head(backgr)

#Map psuedo-absences
#tiff("Figure.tiff",3500,3000, res=350)#this is to save as a tiff file
plot(wrld_simpl, xlim=c(-70,-30), ylim=c(-55,-20), axes=TRUE,col="light yellow")
points(data$Lon, data$Lat, col="orange", pch=20, cex=0.75)
points(backgr, col="black", pch=20, cex=0.75)
#dev.off()

#nNw we repeat the sampling, but limit the area of sampling using a spatial extent
e <- extent(-50, -30, -50, -30)
bg2 <- randomPoints(predictors2, 200, ext=e)

#Map again
#tiff("Figure.tiff",3500,3000, res=350)#this is to save as a tiff file
plot(wrld_simpl, xlim=c(-70,-30), ylim=c(-55,-20), axes=TRUE,col="light yellow")
points(data$Lon, data$Lat, col="orange", pch=20, cex=0.75)
plot(e, add=TRUE, col="red")
points(bg2, col="black", pch=20, cex=0.75)
#dev.off()

#############################################################
#Extract environmental values from absence locations
##############################################################
absvals <- extract(predictors2, backgr)
head(absvals)

#Make a new database 
coords<-as.data.frame(rbind(coords1,backgr))#generate a unique coordinates vector
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)));pb
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals), coords))
head(sdmdata)

###################################################
#    Remove NA
#################################################   
summary(sdmdata)
to.remove <- which(!complete.cases(sdmdata))
sdmdata <- sdmdata[-to.remove,]
summary(sdmdata)

################################################################################
#               Exploration of the dataset
################################################################################

################################################################################
#Step 1: Check correlation among explicative variables. 
#Variables highly correlated can NOT be used together in the model.
################################################################################

matrix<-rcorr(as.matrix(sdmdata[,c(1:9)]), type = "pearson")

# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

corrplot(matrix$r, type="lower", tl.col = "black",method="number",
         p.mat = matrix$P, sig.level = 0.05)

################################################################################
#Step 2: Check if the response variables and the explicative variables have linear relationship
################################################################################
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

g = ggpairs(sdmdata[,c(1:9)], lower = list(continuous = my_fn))
g

################################################################################
#Step 3: Check multicollinearity among variables
################################################################################
source("HighstatLib.r")
corvif(sdmdata[,c(2:9)])
corvif(sdmdata[,c(2,3,4,5,8)])

###############################################################################
#                  Generalized Linear Models - GLMs-
###############################################################################
#Complete model
m1<- glm(pb ~ calcite+chlomean+nitrate+ph+salinity, data=sdmdata,family="binomial")
summary(m1)
D1<-((m1$null.deviance-m1$deviance)/m1$null.deviance)*100;D1
vif(m1)## Evaluate Multi-collinearity

m2<- glm(pb ~ calcite+nitrate+salinity, data=sdmdata,family="binomial")
summary(m2)
D2<-((m2$null.deviance-m2$deviance)/m2$null.deviance)*100;D2
vif(m2)

#Prediction
pglm<-predict(predictors2,m2,type='response')

#map prediction
#tiff("Figure.tiff",3200,3500, res=350)#this is to save as a tiff file
par(mfrow=c(1,1))
plot(pglm)
plot(wrld_simpl, axes=TRUE,add=T,col="grey")

#plot points
points(sdmdata$x,sdmdata$y, col=sdmdata$pb+1,pch=20, cex=0.75)
#dev.off()

##Evaluate prediction#######
coords<-as.data.frame(cbind(sdmdata$x,sdmdata$y))
pp<-extract(pglm,coords)
pp<-as.data.frame(pp)
sdmdata<-cbind(sdmdata,pp)
sdmdata<-na.omit(sdmdata)
cor(sdmdata$pp,sdmdata$pb,method="spearman")

############## Create training datasets################
set.seed(600)
samp <- sample(nrow(sdmdata), round(0.80 * nrow(sdmdata)))
traindata <- sdmdata[samp,]
testdata <- sdmdata[-samp,]

m<- glm(pb ~   calcite+nitrate+salinity,data=traindata, family="binomial")
summary(m)

#Prediction
pglm1<-predict(predictors2,m,type='response')
plot(pglm1)
plot(wrld_simpl, add=TRUE, col="grey")
points(testdata$x,testdata$y,col=testdata$pb+1, pch=16,cex=.70)

eglm <- evaluate(testdata[testdata==1,], testdata[testdata==0,], m);eglm

#TSS  
TSSglm <- mean(eglm@TPR+eglm@TNR-1);TSSglm

#Plot a binary map
tr <- threshold(eglm,"spec_sens");tr

#tiff("Figure_binary.tiff",3200,3500, res=300)#this is to save as a tiff file
plot(pglm1 > tr, main="presence/absence")
plot(wrld_simpl, add=TRUE, col="grey")
points(traindata$x,traindata$y,col=traindata$pb+1, pch=16,cex=.70)
#dev.off()

#################SPATIAL CORRELATION###########################
#Calculate Moran's I values explicitly for a certain distance, 
#and to test for its significance: 
coords <-as.matrix(cbind(dat$x,dat$y))

#The function identifies neighbours of region points by Euclidean distance
#between lower (greater than) and upper (less than or equal to) bounds, 
#or with longlat = TRUE, by Great Circle distance in kilometers.
nb <- dnearneigh(as.matrix(coords), 1,8);nb

plot(wrld_simpl, xlim=c(-70,-30), ylim=c(-55,-20), axes=TRUE,col="light yellow")
plot(nb, coords, col='red', lwd=2, add=TRUE)
listw <- nb2listw(nb,style = "S")#turns neighbourhood object into a weighted list 

#W = row standardized (rows sum to 1), Emphasizes weakly connected points
#B = binary (0/1), Emphasizes strongly connected points
#C = global standardized (all links sum to n), Emphasizes strongly connected points
#S = variance stabilization (Tiefelsdorf et al. 1999), Try to balance

#Null hypothesis that there is no spatial autocorrelation 
#p-valor less than (0.05) the null hypothesis have to be rejected: 
#there is spatial autocorrelation
MoranI <- moran.test(residuals(m2), listw=listw, randomisation=FALSE); MoranI

#Instead of the approach above you should use Monte Carlo simulation. 
#The way it works that the values are randomly assigned to the matrxi, 
#and the Moran???s I is computed. This is repeated several times to establish 
#a distribution of expected values. The observed value of Moran???s I is then 
#compared with the simulated distribution to see how likely it is that the observed 
#values could be considered a random draw.
Moran_MC <- moran.mc(residuals(m2), listw=listw, nsim=100); Moran_MC


