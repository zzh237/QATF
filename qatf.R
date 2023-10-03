# install.packages("devtools")
library(devtools)
install_github("zzh237/detrendr")
install_github("glmgen/genlasso")
# not sure if this next one is necessary - Zhi can you check? 
# install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")

library(detrendr)
# trace(get_model, edit=TRUE)
library(glmgen)

rm(list = ls())

MSE <- function(a, b){
  len = length(a)
  return(norm(a - b, type="2")**2/len) 
}

tau <- c(0.5)
n <- 1000
d <- 10



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

# get scenario 6
# x<- x*0

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
# sum all of the columns to achieve the additive
y_star <- colSums(y_list)
#get the y_i
# get other scenerios error
sce1 <- rnorm(n, 0, 1)
sce2 <- rcauchy(n, 0, 1)
# e <- x**0.5/n**0.5
# te <- rt(n, 2)
# sce3 <- e*te
# sce4 <- rt(n,3)
# sce5 <- rcauchy(n, 0, 1)
y <- y_star + sce1

# get scenerio 6 error
# half <- n/2
# x[1:half] <- (0.25*(x[1:half]/n)**0.5 + 1.375)/3
# x <- replace(x, seq(half+1, n, 1), (7*(tail(x,half)/n)**0.5-2)/3)
# sce6 <-rt(n,2)
# sce6_9q <- qt(0.9, 2)
# y<-y_star+x*sce6
# y_9q <- y_star+x*sce6_9q
# plot(y_9q, type="l", col="black", ylab='0.9 quantile')



# Worth passing in x, y, and y_star? 
get_mse <- function(k, alpha = 10**-6, max_t = 50){
  lambda_list <- 10**seq(5, -10, length.out=100)
  y_mean <- mean(y)
  
  best_mse <- Inf 
  best_lambda <- Inf
  best_trend_hat <- rep(0, times = n)
  
  for (lambda in lambda_list) {
    zero_matrix <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
    trend_list <- as.data.frame(zero_matrix) 
    # we initialize all component functions to be 0
    
    trend_hat_prev <- rep(0, times = n)
    t <- 1
    repeat {
      for (j in 1:d){
        # calculate jth partial residual using components not equal to j
        resp <- y - y_mean - as.numeric(t(rowSums(trend_list[, -j])))
        trend_list[, j] <- as.numeric(trendfilter(x, resp, k = k, lambda = lambda)$beta)
      }
      trend_hat <- rowSums(trend_list)
      
      if (all((abs(trend_hat - trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      trend_hat_prev <- trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, trend_hat + y_mean) 
    cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_trend_hat <- trend_hat
    }
  }
  
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_trend_hat))
}


# Worth passing in x, y, and y_star? 

get_mse_q <- function(k, alpha = 10**-6, max_t = 50){
  lambda_list <- 10**seq(5, -10, length.out=100)
  
  best_mse <- Inf 
  best_lambda <- Inf
  best_q_trend_hat <- rep(0, times = n)
  
  for (lambda in lambda_list) {
    q_trend_list <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
    # we initialize all component functions to be 0
    
    q_trend_hat_prev <- rep(0, times = n)
    t <- 1
    repeat {
      for (j in 1:d){
        # calculate jth partial residual using components not equal to j
        resp <- y - as.numeric(t(rowSums(q_trend_list[, -j])))
        q_trend_list[, j] <- get_trend(resp, tau, lambda, k)
      }
      q_trend_hat <- rowSums(q_trend_list)
      
      if (all((abs(q_trend_hat - q_trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      q_trend_hat_prev <- q_trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, q_trend_hat) 
    cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")
    
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_q_trend_hat <- q_trend_hat
    }
  }
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_q_trend_hat))
}


ATF2 <- get_mse(2)
ATF3 <- get_mse(3)
QATF2 <- get_mse_q(2)
QATF3 <- get_mse_q(3)

cat(ATF2$MSE, "at lambda : ", ATF2$LAMBDA, "\n")
cat(ATF3$MSE, "at lambda : ", ATF3$LAMBDA, "\n")
cat(QATF2$MSE, "at lambda : ", QATF2$LAMBDA, "\n")
cat(QATF3$MSE, "at lambda : ", QATF3$LAMBDA, "\n")


data <- data.frame(n = n, Scenario = 1,
                  tau=tau, QATF2 = QATF2$MSE, QATF3 = QATF3$MSE,
                  ATF2 =ATF2$MSE, ATF3 =ATF3$MSE)

# write data to a sample.csv file
write.table(data, file = "~/QATF/sample2.csv", append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE) 



#> Using same lambda for all quantiles
plot(y, type="l", col="black") # this is the y with errors
plot(y_star, type="l", col="black") # this is the y_star, without errors


# lines(trend_hat~x, col="red")
# lines(q_trend_hat~x, col="blue")
# lines(trend[,1]~x, col="red")
# lines(trend[,2]~x, col="blue")


