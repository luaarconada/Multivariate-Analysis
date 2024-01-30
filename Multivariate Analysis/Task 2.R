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

#Hacemos boxplot para ver las simetrías de las variables
boxplot(numericdata)


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

ggpairs(numericdata, diag = list(continuous = diagonal_hist)) 
#vemos que ya son bastante simetricas todas


#Hacemos factor a las variables cualitativas
factdata=data.frame(lapply(data[c(2:9,19:25)],as.factor))
factdata
head(factdata)
summary(factdata)


#Creamos la dataframe de los datos final con los numéricos estandarizados y los
#categoricos puestos como factor
data=cbind(numericdata,factdata)

########################## Distancias #####################################

#Calcular la matriz de distancias D con la distancia de Gowell
gower_distance<- daisy(data, metric = "gower")

# Access the distance matrix
D2 <- as.matrix(gower_distance)

# Print the distance matrix
D2

#Compute the Gram matrix to see if it's positive semi-definite
gram<-function(d){
  n<-nrow(d)
  ones<-as.matrix(rep(1, n))
  ones_matrix<-ones%*%t(ones)
  identity<- diag(1, nrow = n, ncol = n)
  #Centering matrix
  H<-identity-1/n*ones_matrix
  #Gram matrix
  G<- -1/2* (H%*%d%*%H)
  G
}

G <- gram(D2)
lambda = eigen(G)$values
lambda
min(lambda) >= 0
#Es semidefinida positiva porque todos los autovalores son positivos

rankMatrix(G) #157
U <- t(eigen(G)$vectors)

prin_co <- U %*% sqrt(diag(lambda))
prin_co



# Esta función se utiliza dentro del programa sensitividad_Gower.R que realiza un análisis 
# de sensitividad de la representación MDS obtenida a partir de la distancia de Gower. 
#
# La función [Ycomplet, percentacum] = coorp_complet(D) devuelve las dos primeras coordenadas 
# principales del conjunto de datos Xdat.

coorp_complet <- function(D) {
  n <- nrow(D)
  
  # Comprobamos que D es euclidea (ie, B >= 0)
  H <- diag(n) - matrix(1, n, n) / n
  B <- -H %*% D %*% H / 2
  L <- eigen(B)$values
  m <- min(L)
  epsilon <- 1e-6
  
  if (abs(m) < epsilon) {
    # Hacemos la transformación non2euclid
    D1 <- non2euclid(D)
    B <- -H %*% D1 %*% H / 2
  }
  
  # Cálculo de las coordenadas principales (solo consideramos las tres primeras)
  svd_result <- svd(B)
  T1 <- Re(svd_result$u[, 1:3])  # 'Re' extracts the real part
  # Control de signos de las coordenadas (siempre el primer valor positivo)
  if (T1[1, 1] > 0) {
    T1[, 1] <- -T1[, 1]
  }
  if (T1[1, 2] > 0) {
    T1[, 2] <- -T1[, 2]
  }
  if (T1[1, 3] > 0) {
    T1[, 3] <- -T1[, 3]
  }
  Ycomplet <- T1 %*% sqrt(diag(svd_result$d[1:3]))
  
  percentacum <- sum(diag(svd_result$d[1:2])) / sum(diag(svd_result$d)) * 100  # en la representación en dim 2
  
  return(list(Ycomplet = Ycomplet, percentacum = percentacum))
}

# Function to simulate non2euclid
non2euclid <- function(D) {
  # Replace this function with your non2euclid implementation if available
  # For this example, just returning the input matrix
  return(D)
}

# Example usage:
# Replace the matrix 'D' with your actual matrix
result <- coorp_complet(D2)
print(result$Ycomplet)
print(result$percentacum)


