PCAcovstability <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  H <- diag(n) - matrix(1, n, n) / n
  X <- H %*% X
  S <- cov(X)
  eigen_result <- eigen(S)
  T <- eigen_result$vectors[,order(-eigen_result$values)]
  
  # control de los signos de los componentes
  for (j in 1:p) {
    if (T[1, j] < 0) {
      T[, j] <- -T[, j]
    }
  }
  
  Y <- X %*% T
  
  # leave-one-out
  v <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    Xi <- X[-i, , drop = FALSE]
    Si <- cov(Xi)
    eigen_result_i <- eigen(Si)
    Ti <- eigen_result_i$vectors[,order(-eigen_result_i$values)]
    
    # control de los signos de los componentes
    for (j in 1:p) {
      if (Ti[1, j] < 0) {
        Ti[, j] <- -Ti[, j]
      }
    }
    
    Yi <- Xi %*% Ti
    y <- X[i, ] %*% Ti
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
  
  # PCA unit stability
  boxplot(v, main='Distancia euclidiana de cada unidad entre configuraciones')
  
  lab <- sprintf('%3g', 1:n)
  radius <- Rad[, 5]  # también podría ser Rad[, 6]
  
  theta <- seq(0, 2 * pi, by=0.01)
  
  plot(Y[, 1], Y[, 2], col='blue', pch='.', main='Representación PCA con regiones de 90% de estabilidad')
  abline(h=0, v=0, col='black')  # línea base
  text(Y[, 1], Y[, 2], labels=lab, cex=0.01)
  
  for (i in 1:n) {s
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

PCAcorrstability(num_dataset)