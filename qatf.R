library(devtools)
install_github("zzh237/detrendr")
install_github("glmgen/genlasso")
library(detrendr)
library(genlasso)


rm(list = ls())
MSE <- function(a, b){
  len = length(a)
  return(norm(a - b, type="2")**2/len) 
}

## Need a model selection that penalizes overfitting for our univariate usage


lambda_list <- 10**seq(-1, 4.5, length.out=300)


tau <- 0.5
n <- 60
d <- 3



## prepare the inputs X, and the true function
x <- seq(1, n, 1)

# get sceneria 1, 2, 3
x <- x/n


# get scenerio 4
# half <- n/2
# x[1:half] <- 3*x[1:half]
# x <- replace(x, seq(half+1, n, 1), 3*n - 3*tail(x,half))
# x <- x/n

# get scenerio 5
# x<- cos(6*pi*x/n)

# get scenerio 6
# half <- n/2
# x[1:half] <- (0.25*(x[1:half]/n)**0.5 + 1.375)/3
# x <- replace(x, seq(half+1, n, 1), (7*(tail(x,half)/n)**0.5-2)/3)



# get the simulated heterogeneously-smooth data
x_j <- x # we could add permutation for different j.


# scenerio 6
# get_g_j <- function(x_j, j){
#   g_0 <- (x_j + 0.1)*(j/10)
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j**(x_j + 0.1)*(j/10)
#   return(g_0)
# }

# scenerio 5
# get_g_j <- function(x_j, j){
#   g_0 <- (x_j + 0.1)*(j/10)
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j*(x_j + 0.1)*(j/10)
#   return(g_0)
# }


# scenerio 4
# get_g_j <- function(x_j, j){
#   g_0 <- (x_j + 0.1)**(j/10)
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j*(x_j + 0.1)**(j/10)
#   return(g_0)
# }

# scenerio 1, 2, 3, Doppler-like
get_g_j <- function(x_j, j){
  g_0 <- sin(2*pi/(x_j + 0.1)**(j/10))
  g_0_n <- norm(g_0, type="2")/sqrt(n)
  a_j <- 1/g_0_n
  g_0 <- a_j*sin(2*pi/(x_j + 0.1)**(j/10))
  return(g_0)
}





# create y star and y 
y_list <- data.frame()
for (j in 1:d){
  y_j <- get_g_j(x_j, j) 
  y_list <- rbind(y_list, y_j)
}
#get the y_star
y_list <- unname(y_list)
y_star <- colSums(y_list)
#get the y_i
sce1 <- rnorm(n, 0, 1)
# sce2 <- rcauchy(n, 0, 1)
# e <- x**0.5/n**0.5
# te <- rt(n, 2)
# sce3 <- e*te
# sce4 <- rt(n,3)
# sce5 <- rcauchy(n, 0, 1)
y <- y_star + sce1

#scenario 6
# sce6 <-rt(n,2)
# sce6_9q <- qt(0.9, 2)
# y<-y_star*sce6
# y_9q <- y_star*sce6_9q
# plot(y_9q, type="l", col="black", ylab='0.9 quantile')

y_mean <- mean(y)

get_mse <- function(k, alpha = 0.0001, max_t = 50){
  zero_matrix <- matrix(0, nrow = n, ncol = d)
  trend_list <- as.data.frame(zero_matrix) 
  # we initialize all component functions to be 0
  
  trend_hat_prev <- rep(0, times = n)
  t <- 1
  
  repeat {
    for (j in 1:d){
      # calculate jth partial residual using components not equal to j
      resp <- y - t(rowSums(trend_list[, -j]))
      fit_matrix <- trendfilter(resp, ord=k)$fit 
      mses<-apply(fit_matrix,2,MSE,b=resp)
      print(mses)
      trend_list[, j] <- fit_matrix[,which.min(mses)]
    }
    
    trend_hat <- rowSums(trend_list)
    
    if (all((abs(trend_hat - trend_hat_prev) <= alpha))) {
      cat("terminating after ", t, "iterations from convergence\n")
      break
    } 
    else if (t >= max_t) {
      cat("terminating manually after ", t, "iterations\n")
      break
    }
    trend_hat_prev <- trend_hat
    cat("after ", t, " iterations MSE = ", MSE(y_star, trend_hat), "\n")
    t <- t + 1
  }
  return(MSE(y_star, trend_hat))
}



get_mse_q <- function(k, alpha = 0.0001, max_t = 50){
  zero_matrix <- matrix(0, nrow = n, ncol = d)
  q_trend_list <- as.data.frame(zero_matrix) 
  # we initialize all component functions to be 0
  
  q_trend_hat_prev <- rep(0, times = n)
  t <- 1
  
  repeat {
    for (j in 1:d){
      # calculate jth partial residual using components not equal to j
      resp <- y - as.numeric(t(rowSums(q_trend_list[, -j])))
      old_MSE <- Inf
      for (lambda in lambda_list) {
        fit <- get_trend(resp, tau, lambda, k)
        if (MSE(fit, resp) < old_MSE) {
          best_fit <- fit
          old_MSE <- MSE(fit, resp)
        }  
      }
      q_trend_list[, j] <- best_fit
    }
    
    q_trend_hat <- rowSums(q_trend_list)
    
    if (all(abs(q_trend_hat - q_trend_hat_prev) <= alpha)) {
      cat("terminating after ", t, "iterations from convergence\n")
      break
    } 
    else if (t >= max_t) {
      cat("terminating manually after ", t, "iterations\n")
      break
    }
    q_trend_hat_prev <- q_trend_hat
    cat("after ", t, " iterations MSE = ", MSE(y_star, q_trend_hat), "\n")
    t <- t + 1
  }
  return(MSE(y_star, q_trend_hat))
}


ATF2 = get_mse(k=2)
ATF2
QATF2 = get_mse_q(k=2)
QATF2
ATF3 = get_mse(k=3)
ATF3
QATF3 = get_mse_q(k=3)
QATF3


data = data.frame(n = n, Scenario = 6,
                  tau=tau, QATF2 = QATF2, QATF3 = QATF3,
                  ATF2 =ATF2, ATF3 =ATF3)

# write data to a sample.csv file
write.table(data, file = "sample2.csv", append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE) 



#> Using same lambda for all quantiles
plot(y, type="l", col="black") # this is the y with errors
plot(y_star, type="l", col="black") # this is the y_star, without errors


# lines(trend_hat~x, col="red")
# lines(q_trend_hat~x, col="blue")
# lines(trend[,1]~x, col="red")
# lines(trend[,2]~x, col="blue")



