library(kableExtra)
library(calibrate)
library(e1071)
library(Matrix)
library(ggplot2)
library(GGally)
library(car)
library(dplyr)
library(corrplot)

names = c('age', 'bp','sg', 'al', 'su', 'rbc', 'pc', 'pcc', 'ba', 'bgr', 'bu', 'sc',
          'sod', 'pot', 'hemo', 'pcv', 'wc', 'rc', 'htn', 'dm', 'cad', 'appet',
          'pe', 'ane', 'class')
data =read.table("chronic_kidney_disease.csv", sep=';', na.strings='?', header=F, col.names=names)

data =na.omit(data)
head(data)
data
summary(data)


#Seleccionamos las variables numéricas
numericdata<- data[c(1,10:18)]
numericdata
head(numericdata)
summary(numericdata)


#DATA VISUALIZATION AND TRANSFORMATION (SIMMETRY AND STANDARIZATION)
# Función para personalizar la diagonal con histogramas
diagonal_hist <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(fill = "lightblue", color = "black", bins =5 ) +
    theme_minimal()
}

# Crear el gráfico de pares con histogramas en la diagonal usando GGally
ggpairs(numericdata, diag = list(continuous = diagonal_hist))
?ggpairs

#Hacemos boxplot para ver las simetrías de las variables
boxplot(numericdata)
?boxplot


#Miramos la simetría así que vemos la skewness
apply(numericdata, MARGIN=2, FUN=skewness)

#Las que tienen mucha skewness positiva las transformamos con el logaritmo
numericdata[,c(2,3,4,9)]=1/(numericdata[,c(2,3,4,9)])
boxplot(numericdata)

#Las que tienen skewness negativa las transformamos con el cuadrado 
numericdata[,c(5,7,8)]=(numericdata[,c(5,7,8)])^2
boxplot(numericdata)


apply(numericdata, MARGIN=2, FUN=skewness)


#Qué mala escala! Vamos a estandarizar las variables
numericdata<-as.data.frame(scale(numericdata,center=TRUE,scale=TRUE))
boxplot(numericdata)

ggpairs(numericdata, diag = list(continuous = diagonal_hist)) #vemos que ya son bastante simetricas todas

#Creamos matriz y pintamos
matriz=as.matrix(numericdata)
matriz
pairs(matriz)


#Vector de medias
means<- colMeans(numericdata)
means

#Matriz de covarianzas
S<-cov(numericdata)
S
rankS<-qr(S)$rank
rankS
#rank=10, así que tenemos 10 variables linearmente independientes

#Matrix correlaciones
R<-cor(numericdata)
R
rankR<-qr(R)$rank
rankR
heatmap(R,Colv=NA,Rowv=NA)

#Intercorrelation measures
intercorrelation <- function(X){
  q<-rep(0,6)
  p<-ncol(X)
  R=cor(X)
  q[1]<-(1-min(eigen(R)$value)/(max(eigen(R)$value)))^(p+2)
  q[2]<- 1-p/(sum(1/eigen(R)$value))
  q[3]<-1-sqrt(abs(det(R)))
  q[4]<-(max(eigen(R)$value)/p)^(3/2)
  q[5]<-(1-min(eigen(R)$value)/p)^5
  R_1=solve(R)
  q[6]<-sum((1-(1/diag(R_1)))/p) #este no se si esta bien definido
  q
}

intercorrelation(numericdata)
corrplot(R, method="color",  col = COL2('PiYG'), tl.col = "darkgreen", addCoef.col ='black',number.cex = 0.55)

#We see the relation along quantitative variables among the types of appet.
appet<-split(numericdata,data$appet)
poor<-appet$poor
good<-appet$good
interpoor<-intercorrelation(poor)
interpoor
intergood<-intercorrelation(good)
intergood

plot(1:6,interpoor,col='red',ylim=c(0.15,1))
lines(1:6,interpoor,col='red')
points(1:6,intergood,col='blue')
lines(1:6,intergood,col='blue')
legend(1,0.8,legend=c('poor','good'),col=c('red','blue'),lty=1:1,cex=1,bty='n')
#The red is poor and the blue is good (añadair esta aclaración en RMarkdown)

# Create the scatterplot matrix with diagonal histograms
ggpairs(numericdata, 
        diag = list(continuous = diagonal_hist),
        aes(color = data$appet))+
  scale_color_manual(values = c("red", "blue"))



############### We start with the PCA #############################

#Calculamos la matriz centrada
n<-nrow(matriz)
n
p<-ncol(matriz)
p

#Crea una nueva matriz utilizando el vector de medias. n es el número de filas
# en la matriz original y ncol(matriz) es el número de columnas. byrow = TRUE 
# indica que los valores del vector de medias se deben repetir por fila para
# formar la matriz.
matrizcentrada <- matriz-1/n*matrix(means, n, p, byrow =TRUE)
matrizcentrada

eigenlista<-eigen(R)
lambda<-eigenlista$values
lambda
V<-eigenlista$vectors
V
#Comprobamos que V es ortogonal
V %*% t(V)

##Dimensionality reduction 
#Variabilidad explicada
varexpm<-lambda*(100/sum(lambda))
varexpm
varexp<-diag(varexpm)
varexpacum<-cumsum(varexp)
varexpacum
typeof(varexpacum)
vectorcum<- varexpacum[c(10,20,30,40,50,60,70,80,90,100)]
vectorcum

#Scree-plot with % variability
plot(1:length(lambda), diag(varexp), type = "b", pch = 19, xlab = "Eigenvalue", ylab = "% variability", main = "Scree Plot")


#PCA a partir de R 2 componentes (54.70750%)
Y<-matriz %*% V
Y
Y[,1:2]

plot(Y[,1],Y[,2],xlab='Y1',ylab='Y2',main='PCA from R (54.70750%)',col='blue')
textxy(Y[,1],Y[,2],labs=1:156,cex=0.5,offset=1)


#Methods for discarding components
#1.Percentage of explained variability
#Scree-plot with % variability
plot(1:length(lambda), diag(varexp), type = "b", pch = 19, xlab = "Eigenvalue", ylab = "% variability", main = "Scree Plot")

#Kaiser's & Jollife's criteria
#Scree-plot with criteria
plot(1:length(lambda), lambda, type = "b", pch = 19, ylab = "Eigenvalue",xlab='', 
     main = "Scree-Plot with Kaiser's and Jollife's criteria")
abline(h = mean(lambda), col = "red")
abline(h=0.7*mean(lambda), col="blue") 
#Cattell's scree graph is used as a visual tool and we usually take the
#principal components with the steepest slopes
#We can not perform dimensionality test because the principal components are
#computed from the correlation matrix

#Summary of data to choose our principal components
summary=data.frame(rbind(lambda,diag(varexp),vectorcum),
                   row.names=c('eigenvalues','% variability','cumulated variability'))
colnames(summary)<-c('$Y_1$','$Y_2$','$Y_3$','$Y_4$','$Y_5$','$Y_6$','$Y_7$','$Y_8$','$Y_9$','$Y_10$')
summary
head(summary)
#Meter kable(summary) al hacer el RMarkdown
#Usar esta tabla para explicar con cuántos cpmponentes nos queadaríamos según cada criterio.


# Barplot pesos
names=c('PC1 weights','PC2 weights','PC3 weights','PC4 weights','PC5 weights')
barplot(V[1:5,],beside=TRUE, col = rainbow(5),legend.text = names, 
        args.legend = list(x = "topright", inset=c(-0,-0.2), y=0.7, bty = "n"),
        ylim=c(-1,1), axis.lty=1)
heatmap(V[,1:5],Colv=NA,Rowv=NA)

#Analysis of stability
PCAcovstability <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  H <- diag(n) - matrix(1, n, n) / n
  X <- H %*% X
  S <- cov(X)
  eigen_result <- eigen(S)
  V <- eigen_result$vectors[,order(-eigen_result$values)]
  
# control de los signos de los componentes
for (j in 1:p) {
  if (V[1, j] < 0) {
    V[, j] <- -V[, j]
    }
  }
  
Y <- X %*% V
  
# leave-one-out
v <- matrix(0, n, n)
for (i in 1:(n-1)) {
  Xi <- X[-i, , drop = FALSE]
  Si <- cov(Xi)
  eigen_result_i <- eigen(Si)
  Vi <- eigen_result_i$vectors[,order(-eigen_result_i$values)]    
  # control de los signos de los componentes
  for (j in 1:p) {
    if (Vi[1, j] < 0) {
      Vi[, j] <- -Vi[, j]
    }
  }
    
  Yi <- Xi %*% Vi
  y <- X[i, ] %*% Vi
  Yplusi <- rbind(Yi[1:(i-1), , drop = FALSE], y, Yi[(i):(n-1), , drop = FALSE])
    
  # Distancia euclidiana de cada unidad entre configuraciones
  for (k in 1:n) {
    v[k, i] <- sum((Yplusi[k, , drop = FALSE] - Y[k, , drop = FALSE])^2)
    }
}
  
# La i-ésima fila en la matriz v contiene la configuración para la i-ésima unidad
v <- sqrt(v)
  
# Estadísticas resumen para las filas de v
Rad <- matrix(0, n, 6)
for (i in 1:n) {
  Rad[i, 1] <- min(v[i, ])
  Rad[i, 2] <- max(v[i, ])
  Rad[i, 3] <- quantile(v[i, ], 0.50)
  Rad[i, 4] <- quantile(v[i, ], 0.75)
  Rad[i, 5] <- quantile(v[i, ], 0.90)
  Rad[i, 6] <- quantile(v[i, ], 0.95)
}
  
# Algunas estadísticas resumen por filas
rowsummary <- rbind(mean(Rad), sd(Rad), apply(Rad, 2, median), apply(Rad, 2, mad))
rowsummary<- as.matrix(rowsummary)
colnames(rowsummary)<- c("min", "max", "median", "75th-p", "90-th p", "95-th p")
rownames(rowsummary)<- c("mean", "median", "std", "mad")
return(rowsummary)
# PCA unit stability
# par(mfrow=c(2, 1))
# par(mar=c(4, 4, 2, 2))
print(boxplot(v, main='Distancia euclidiana de cada unidad entre configuraciones'))
  
lab <- sprintf('%3g', 1:n)
radius <- Rad[, 5]  # también podría ser Rad[, 6]
  
theta <- seq(0, 2 * pi, by=0.01)
  
plot(Y[, 1], Y[, 2], col='blue', pch='.', main='Representación PCA con regiones de 90% de estabilidad')
abline(h=0, v=0, col='black')  # línea base
text(Y[, 1], Y[, 2], labels=lab, cex=0.01)

  
for (i in 1:n) {
  circle <- matrix(0, length(theta), 2)
  for (j in 1:length(theta)) {
    circle[j, 1] <- Y[i, 1] + cos(theta[j]) * radius[i]
    circle[j, 2] <- Y[i, 2] + sin(theta[j]) * radius[i]
  }
  points(circle[, 1], circle[, 2], col='black', pch='.', cex=0.1)
}
}

PCAcorrstability <- function(X) {
  Z <- scale(X)
  rowsummary <- PCAcovstability(Z)
  return(rowsummary)
}

PCAcorrstability(numericdata)  

