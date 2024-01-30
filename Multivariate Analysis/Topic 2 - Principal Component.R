
#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&
#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&
# Topic 2 - Principal Component Analysis
#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&
#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&#&

##################################################################################################################
# Open R or R-studio.
##################################################################################################################

##########################################################################################################
# Select the number of digits in the outputs

options(digits=4)

##################################################################################################################
##################################################################################################################
# Principal Component Analysis
##################################################################################################################
##################################################################################################################

##################################################################################################################
# The NCI60 data set
##################################################################################################################

##################################################################################################################
# Install and load ISLR package

install.packages("ISLR")
library("ISLR")

# See the help page

?NCI60

# Load the data set into memory

data(NCI60)

# Remember that NCI60 is a list with two objects

summary(NCI60)

# The data matrix is the object called data

X <- NCI60$data

##################################################################################################################
# Sample sizes and dimensions of the NCI60 data set

n.NCI60 <- nrow(X)
n.NCI60
p.NCI60 <- ncol(X)
p.NCI60

##################################################################################################################
# PCA for the NCI60 data set
##################################################################################################################

##################################################################################################################
# There are several functions in R to perform PCA. Next, we use the function prcomp

PCS.NCI60 <- prcomp(X)

# Have a look at the outputs

names(PCS.NCI60)

##################################################################################################################
# The PC scores is are in x. In this case, the scores are a matrix of size 64x64

dim(PCS.NCI60$x)
head(PCS.NCI60$x)

##################################################################################################################
# Make a plot of the first two PC scores

plot(PCS.NCI60$x[,1:2],pch=19,col="deepskyblue2",main="First two PCs")

# Remember that one of our interests are in grouping cancer cell lines
# The first two PCs suggest the presence of at least two different groups
# This information is completely hidden if we consider the 6830 variables

# We can add the names of the cancer type

text(PCS.NCI60$x[,1:2],labels=NCI60$labs,pos=1,col="firebrick2",cex=0.7)

##################################################################################################################
# The eigenvalues of the sample covariance matrix of X, i.e., the variances of the PCs are the square of sdev
# As in this example n<p, only 64 eigenvalues can be different from 0. These are those that appear here.

PCS.NCI60$sdev^2

##################################################################################################################
# Have a look at these eigenvalues

install.packages("factoextra")
library("factoextra")

# Screeplot with the 64 eigenvalues

fviz_eig(PCS.NCI60,ncp=64,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")

# Screeplot with the first 30 eigenvalues

fviz_eig(PCS.NCI60,ncp=30,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")

##################################################################################################################
# How many PCs are important?

# Have a look at the proportion of explained variance and the cumulative proportion of explained variance

get_eigenvalue(PCS.NCI60)

# We need 21 PCs to have the 70% of the total variability of the data set
# Have a look at the sample mean of the eigenvalues

eval.NCI60 <- PCS.NCI60$sdev^2
mean(eval.NCI60)

# The number of eigenvalues larger than this sample mean is 

sum(eval.NCI60>mean(eval.NCI60))

# Then, the dimension of the data set is reduced from 6830 to either 21 or 17 that represents either the 0.30% 
# or the 0.24% of the number of variables in the original data set keeping around the 70% of the information inside

##################################################################################################################
# The loading matrix, i.e., the eigenvectors of the sample covariance matrix of X are given in rotation
# The output is a matrix of size 6830x64, so we only have a look at the first few rows

dim(PCS.NCI60$rotation)
head(PCS.NCI60$rotation)

# Note that only 64 eigenvectors (corresponding to 64 PCs) appear

##################################################################################################################
# Interpretation of the PCs

# The first eigenvector has 6830 values
# Each value represent the weight of the associated variable (the measurement on a particular gen)
# Have a look at the important variables (genes) in the first PC

plot(PCS.NCI60$rotation[,1],pch=20,col="deepskyblue2",main="Weights for the first PC")
abline(h=0)

# Have a look at the important variables (genes) in the second PC

plot(PCS.NCI60$rotation[,2],pch=20,col="deepskyblue2",main="Weights for the second PC")
abline(h=0)

# Have a look at the important variables (genes) in the first two PCs

plot(PCS.NCI60$rotation[,1:2],pch=19,col="deepskyblue2",main="Weights for the first two PCs")
abline(h=0,v=0)
text(PCS.NCI60$rotation[,1:2],labels=colnames(NCI60$data),pos = 1,col="firebrick2",cex=0.5)

# Those genes that appear to be outliers are those genes that creates the largest variability in
# the first two PCs, which are the most important ones
# Add a circle to point out the outliers. Note that the radius is arbitrary

install.packages("plotrix")
library("plotrix")

draw.circle(0,0,0.05,border="green2",lwd=3)

# These genes are those more related with the cancer cell lines

##################################################################################################################
# Make a plot of the first three PCs

pairs(PCS.NCI60$x[,1:3],col="deepskyblue2",pch=19,main="The first three PCs")

# Have a look at the important variables (genes) in the third PC

plot(PCS.NCI60$rotation[,3],pch=19,col="deepskyblue2",main="Weights for the third PC")
abline(h=0)

# Have a look at the important variables (genes) in the first three PCs

pairs(PCS.NCI60$rotation[,1:3],pch=19,col="deepskyblue2",main="Weights for the first three PCs")

##################################################################################################################
# PCA for the College data set
##################################################################################################################

##################################################################################################################
# Load the College data set

# College data set

?College

# Load data set

data(College)

# Have a look at the first rows

head(College)

# Size of the matrix

dim(College)

# We have 777 colleges and 18 variables

##################################################################################################################
# We define the data matrix with the quantitative variables and the indicator vector with the qualitative variable (Private)

X <- College[,2:18]
head(X)

Y <- College[,1]
head(Y)
table(Y)

##################################################################################################################
# Define sample size and the dimension

n.X <- nrow(X)
n.X
p.X <- ncol(X)
p.X

##################################################################################################################
# Plot the original variables

library(MASS)
colors.X <- c("deepskyblue2","firebrick2")[Y]
parcoord(X,col=colors.X,var.label=TRUE)

##################################################################################################################
# Two possibilities: Obtain PCs on the whole data set or on the two groups
# See what happens if we obtain PCs on the whole data set
# Use the sample correlation matrix

PCS.X <- prcomp(X,scale=TRUE)

##################################################################################################################
# PC scores

dim(PCS.X$x)
head(PCS.X$x)

##################################################################################################################
# Make a plot of the first two PCs

plot(PCS.X$x[,1:2],pch=19,col=colors.X)

# The first two PCs show the presence of the two different groups
# Indeed the relationship between the first and the second PC is totally different in terms of the indicator variable
# This information was not clear in a plot of the 17 variables

pairs(X,pch=20,col=colors.X)

##################################################################################################################
# The eigenvalues of the sample correlation matrix of X, i.e., the variances of the PCs are the square of sdev

PCS.X$sdev^2

##################################################################################################################
# Have a look at these eigenvalues

# Screeplot

fviz_eig(PCS.X,ncp=17,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")

##################################################################################################################
# How many PCs are important?

# Have a look at the proportion of explained variance and the cumulative proportion of explained variance

get_eigenvalue(PCS.X)

# We need 4 PCs to have the 70% of the total variability of the data set
# Check that the sample mean of the eigenvalues is 1

eval.X <- PCS.X$sdev^2
mean(eval.X)

# The number of eigenvalues larger than this sample mean is 

sum(eval.X > mean(eval.X))

# Then, the dimension of the data set is reduced from 17 to 4 that represents either the 23.52% of the number of 
# variables in the original data set keeping around the 70% of the information inside

##################################################################################################################
# The loading matrix, i.e., the eigenvectors of the sample correlation matrix of X are given in rotation.
# See the loadings of the first four PCs

dim(PCS.X$rotation)
PCS.X$rotation[,1:4]

##################################################################################################################
# Interpretation of the first PC: Weights for the first PC

plot(1:p.X,PCS.X$rotation[,1],pch=19,col="deepskyblue2",main="Weights for the first PC")
abline(h=0)
text(1:p.X,PCS.X$rotation[,1],labels=colnames(X),pos=1,col="firebrick2",cex=0.75)

##################################################################################################################
# Interpretation of the second PC: Weights for the second PC

plot(1:p.X,PCS.X$rotation[,2],pch=19,col="deepskyblue2",main="Weights for the second PC")
abline(h=0)
text(1:p.X,PCS.X$rotation[,2],labels=colnames(X),pos=1,col="firebrick2",cex=0.75)

##################################################################################################################
# Have a look at the important variables in the first two PCs
# Note the different groups in the data
# The radius is arbitrary

plot(PCS.X$rotation[,1:2],pch=19,col="deepskyblue2",main="Weights for the first two PCs")
abline(h=0,v=0)
text(PCS.X$rotation[,1:2],labels=colnames(X),pos=1,col="firebrick2",cex=0.75)
draw.circle(0,0,0.3,border="green2",lwd=3)

##################################################################################################################
# The biplot is an alternative way to plot points and the first two PCs together 
# However, it is only useful when the data set is not too large

biplot(PCS.X,col=c("deepskyblue2","firebrick2"),cex=c(0.5,0.8))

##################################################################################################################
# Plot the scores (the four PCs)

pairs(PCS.X$x[,1:4],col=colors.X,pch=19,main="The first four PCs")

# The PCs show that the variables have different behavior in terms of the two groups

##################################################################################################################
# Plot the correlations between the original data set and the scores
# This is useful to understand also which are the important variables in the data set in terms of variability

install.packages("corrplot")
library("corrplot")

corrplot(cor(X,PCS.X$x),is.corr=T)

# Reduce to only four PCs

corrplot(cor(X,PCS.X$x[,1:4]),is.corr=T)

##################################################################################################################
# Compress images with PCA
##################################################################################################################

##################################################################################################################
# Set working directory

setwd("C:/temp")

##################################################################################################################
# Install and load jpeg package

install.packages("jpeg")
require("jpeg")

##################################################################################################################
# Load image

image.jpg <- readJPEG('image.jpg')

##################################################################################################################
# Check that image.jpg is an array of dimension 354x630x3

is(image.jpg)
dim(image.jpg)
n <- ncol(image.jpg)
n
p <- nrow(image.jpg)
p

##################################################################################################################
# Therefore, there are three RGB matrices corresponding to Red, Green and Blue

image.jpg.red <- image.jpg[,,1]
image.jpg.green <- image.jpg[,,2]
image.jpg.blue <- image.jpg[,,3]

# Check that the dimension of the matrices is 354x630

dim(image.jpg.red)
dim(image.jpg.green)
dim(image.jpg.blue)

# See the first rows of the matrix for colour red and check that all the pixels goes from 0 to 1

head(image.jpg.red)
min(image.jpg.red)
max(image.jpg.red)

##################################################################################################################
# PCA on the RGB matrices

image.jpg.red.pca <- prcomp(image.jpg.red,center=FALSE)
image.jpg.green.pca <- prcomp(image.jpg.green,center=FALSE)
image.jpg.blue.pca <- prcomp(image.jpg.blue,center=FALSE)

##################################################################################################################
# Screeplots and eigenvalues

fviz_eig(image.jpg.red.pca,ncp=20,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")
get_eigenvalue(image.jpg.red.pca)

fviz_eig(image.jpg.green.pca,ncp=20,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")
get_eigenvalue(image.jpg.green.pca)

fviz_eig(image.jpg.blue.pca,ncp=20,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")
get_eigenvalue(image.jpg.blue.pca)

# With a few PCs we can recover more than the 90% of the total variability of the data matrices

##################################################################################################################
# Compression factors

r <- c(1,2,3,4,5,10,100)
comp.fact <- r/p + r/n
comp.fact

# The compression factors are smaller than 1 even if r is 100

##################################################################################################################
# Have a look at the images created for these values of r

# Install and load imager package

install.packages("imager")
library("imager")

# Create a list with the three PCs

image.jpg.pca <- list(image.jpg.red.pca,image.jpg.green.pca,image.jpg.blue.pca)

# Case r=1

r <- 1
image.pca.1 <- sapply(image.jpg.pca,function(ijp){comp.img <- ijp$x[,1:r] %*% t(ijp$rotation[,1:r])},simplify='array')
image.pca.1 <- as.cimg(image.pca.1,x=n,y=p,z=1,cc=3)
imager::save.image(image.pca.1,"C:/temp/image-1.jpg")

# Case r=2

r <- 2
image.pca.2 <- sapply(image.jpg.pca,function(ijp){comp.img <- ijp$x[,1:r] %*% t(ijp$rotation[,1:r])},simplify='array')
image.pca.2 <- as.cimg(image.pca.2,x=n,y=p,z=1,cc=3)
imager::save.image(image.pca.2,"C:/temp/image-2.jpeg")

# Case r=3

r <- 3
image.pca.3 <- sapply(image.jpg.pca,function(ijp){comp.img <- ijp$x[,1:r] %*% t(ijp$rotation[,1:r])},simplify='array')
image.pca.3 <- as.cimg(image.pca.3,x=n,y=p,z=1,cc=3)
imager::save.image(image.pca.3,"C:/temp/image-3.jpeg")

# Case r=4

r <- 4
image.pca.4 <- sapply(image.jpg.pca,function(ijp){comp.img <- ijp$x[,1:r] %*% t(ijp$rotation[,1:r])},simplify='array')
image.pca.4 <- as.cimg(image.pca.4,x=n,y=p,z=1,cc=3)
imager::save.image(image.pca.4,"C:/temp/image-4.jpeg")

# Case r=5

r <- 5
image.pca.5 <- sapply(image.jpg.pca,function(ijp){comp.img <- ijp$x[,1:r] %*% t(ijp$rotation[,1:r])},simplify='array')
image.pca.5 <- as.cimg(image.pca.5,x=n,y=p,z=1,cc=3)
imager::save.image(image.pca.5,"C:/temp/image-5.jpeg")

# Case r=10

r <- 10
image.pca.10 <- sapply(image.jpg.pca,function(ijp){comp.img <- ijp$x[,1:r] %*% t(ijp$rotation[,1:r])},simplify='array')
image.pca.10 <- as.cimg(image.pca.10,x=n,y=p,z=1,cc=3)
imager::save.image(image.pca.10,"C:/temp/image-10.jpeg")

# Case r=100

r <- 100
image.pca.100 <- sapply(image.jpg.pca,function(ijp){comp.img <- ijp$x[,1:r] %*% t(ijp$rotation[,1:r])},simplify='array')
image.pca.100 <- as.cimg(image.pca.100,x=n,y=p,z=1,cc=3)
imager::save.image(image.pca.100,"C:/temp/image-100.jpeg")