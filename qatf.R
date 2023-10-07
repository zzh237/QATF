# install.packages("devtools")
library(devtools)
install_github("zzh237/detrendr")
install_github("glmgen/genlasso")
# not sure if this next one is necessary - Zhi can you check? 
# install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")

library(detrendr)
# trace(get_model, edit=TRUE)
library(glmgen)

# rm(list = ls())
# cum_data <- data.frame()
# read.csv("MSEs.csv")

MSE <- function(a, b){
  len = length(a)
  return(norm(a - b, type="2")**2/len) 
}

# These functions require inputs to be passed in
get_mse <- function(y, y_star, n, d, k, alpha = 10**-6, max_t = 50, prints = TRUE){
  lambda_list <- 10**seq(5, -10, length.out=50)
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
        trend_list[, j] <- as.numeric(trendfilter(seq(1, n, 1)/n, resp, k = k, lambda = lambda)$beta)
      }
      trend_hat <- rowSums(trend_list)
      
      if (all((abs(trend_hat - trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      trend_hat_prev <- trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, trend_hat + y_mean) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")}
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_trend_hat <- trend_hat
    }
  }
  
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_trend_hat))
}
get_mse_q <- function(y, y_star, n, d, tau, k, alpha = 10**-6, max_t = 50, prints = TRUE){
  lambda_list <- 10**seq(5, -10, length.out=50)
  
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
        # get_trend requires equally spaced points
      }
      q_trend_hat <- rowSums(q_trend_list)
      
      if (all((abs(q_trend_hat - q_trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      q_trend_hat_prev <- q_trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, q_trend_hat) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")}
    
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_q_trend_hat <- q_trend_hat
    }
  }
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_q_trend_hat))
}


scenario1 <- function(n, d, tau) {
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }

  x <- seq(1, n, 1)/n
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Doppler-like
    g_0 <- sin(2*pi/(x + 0.1)**(j/10))
    a_j <- 1 / (norm(g_0, type="2")/sqrt(n))
    b_j <- a_j*mean(g_0)
    y_j <- a_j*g_0 - b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  y <- y_star + rnorm(n, 0, 1)
  y_star_q <- y_star + qnorm(tau, 0, 1)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, type = "l", col = "black", ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario2 <- function(n, d, tau) {
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- seq(1, n, 1)/n
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Doppler-like
    g_0 <- sin(2*pi/(x + 0.1)**(j/10))
    a_j <- 1 / (norm(g_0, type="2")/sqrt(n))
    b_j <- a_j*mean(g_0)
    y_j <- a_j*g_0 - b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Cauchy Errors
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, type = "l", col = "black", ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario3 <- function(n, d, tau) {
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- seq(1, n, 1)/n
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Doppler-like
    g_0 <- sin(2*pi/(x + 0.1)**(j/10))
    a_j <- 1 / (norm(g_0, type="2")/sqrt(n))
    b_j <- a_j*mean(g_0)
    y_j <- a_j*g_0 - b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Heteroskedastic T errors
  e <- (x**.5)  # i**.5 / n**.5
  y <- y_star + e * rt(n, 2)
  y_star_q <- y_star + e * qt(tau, 2)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, type = "l", col = "black", ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}

# Haven't done these yet
scenario4 <- function(n, d, tau) {
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- seq(1, n, 1)/n
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Doppler-like
    g_0 <- sin(2*pi/(x + 0.1)**(j/10))
    a_j <- 1 / (norm(g_0, type="2")/sqrt(n))
    b_j <- a_j*mean(g_0)
    y_j <- a_j*g_0 - b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Heteroskedastic T errors
  e <- (x**.5)  # i**.5 / n**.5
  y <- y_star + e * rt(n, 2)
  y_star_q <- y_star + e * qt(tau, 2)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, type = "l", col = "black", ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario5 <- function(n, d, tau) {
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- seq(1, n, 1)/n
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Doppler-like
    g_0 <- sin(2*pi/(x + 0.1)**(j/10))
    a_j <- 1 / (norm(g_0, type="2")/sqrt(n))
    b_j <- a_j*mean(g_0)
    y_j <- a_j*g_0 - b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Heteroskedastic T errors
  e <- (x**.5)  # i**.5 / n**.5
  y <- y_star + e * rt(n, 2)
  y_star_q <- y_star + e * qt(tau, 2)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, type = "l", col = "black", ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario6 <- function(n, d, tau) {
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- seq(1, n, 1)/n
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Doppler-like
    g_0 <- sin(2*pi/(x + 0.1)**(j/10))
    a_j <- 1 / (norm(g_0, type="2")/sqrt(n))
    b_j <- a_j*mean(g_0)
    y_j <- a_j*g_0 - b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Heteroskedastic T errors
  e <- (x**.5)  # i**.5 / n**.5
  y <- y_star + e * rt(n, 2)
  y_star_q <- y_star + e * qt(tau, 2)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, type = "l", col = "black", ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
# Need to address scenario6 issues


run_custom_sce_simulations <- function(n, d, tau, sce, simulations = 1) {
  # This function only exists for us to make constructing the data frame quicker
  # Wrapper function is set to run with fixed k = 1 and 2, 
  # and append results to formatted data frame.

  ATF1_MSE <- 0
  ATF2_MSE <- 0
  QATF1_MSE <- 0
  QATF2_MSE <- 0
  for (i in 1:simulations) {
    if      (sce == 1) {vals <- scenario1(n, d, tau)}
    else if (sce == 2) {vals <- scenario2(n, d, tau)}
    else if (sce == 3) {vals <- scenario3(n, d, tau)}
    else if (sce == 4) {vals <- scenario4(n, d, tau)}
    else if (sce == 5) {vals <- scenario5(n, d, tau)}
    else if (sce == 6) {vals <- scenario6(n, d, tau)}
    else {stop("Only 6 scenarios at the time of this functions' construction")}
    
    # Currently preserving the full output in case we want cool plots
    # Also, only run get_mse for tau = 0.5
    if (tau == 0.5) {
      ATF1 <- get_mse(vals[[1]], vals[[2]], n, d, 1, prints = FALSE)
      ATF2 <- get_mse(vals[[1]], vals[[2]], n, d, 2, prints = FALSE)
      ATF1_MSE <- ATF1_MSE + ATF1$MSE
      ATF2_MSE <- ATF2_MSE + ATF2$MSE
      
      cat("mse for ATF1 was ", ATF1$MSE, "at lambda = ", ATF1$LAMBDA, "\n")
      cat("mse for ATF2 was ", ATF2$MSE, "at lambda = ", ATF2$LAMBDA, "\n")
    } else {
      cat("not running ATF1 or ATF2 since tau != 0.5")
    }
    
    QATF1 <- get_mse_q(vals[[1]], vals[[2]], n, d, tau, 1, prints = FALSE)
    QATF2 <- get_mse_q(vals[[1]], vals[[2]], n, d, tau, 2, prints = FALSE)
    QATF1_MSE <- QATF1_MSE + QATF1$MSE
    QATF2_MSE <- QATF2_MSE + QATF2$MSE
    
    cat("mse for QATF1 was ", QATF1$MSE, "at lambda = ", QATF1$LAMBDA, "\n")
    cat("mse for QATF2 was ", QATF2$MSE, "at lambda = ", QATF2$LAMBDA, "\n")
    if (simulations != 1) {cat("finished simulation ", i, "\n")}
  }
  
  if (tau == 0.5) {
    ATF1_MSE <- ATF1_MSE / simulations
    ATF2_MSE <- ATF2_MSE / simulations
  } else {
    ATF1_MSE <- NA
    ATF2_MSE <- NA
  }
  QATF1_MSE <- QATF1_MSE / simulations
  QATF2_MSE <- QATF2_MSE / simulations
  
  return(data.frame(n = n, Scenario = sce,
                    tau=tau, Simulations = simulations, 
                    QATF1 = QATF1_MSE, QATF2 = QATF2_MSE,
                    ATF1 = ATF1_MSE, ATF2 = ATF2_MSE))
  
}

cum_data <- data.frame()
cum_data <- rbind(cum_data, run_custom_sce_simulations(500, 10, 0.5, 1, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 1, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 1, simulations = 1))




vals <- scenario1(500, 10, 1, 0.5)
ATF1 <- get_mse(vals[[1]], vals[[2]], 500, 10, 1)
ATF2 <- get_mse(vals[[1]], vals[[2]], 500, 10, 2)
QATF1 <- get_mse_q(vals[[1]], vals[[2]], 500, 10, 0.5, 1)
QATF2 <- get_mse_q(vals[[1]], vals[[2]], 500, 10, 0.5, 2)

cat(ATF1$MSE, "at lambda : ", ATF1$LAMBDA, "\n")
cat(ATF2$MSE, "at lambda : ", ATF2$LAMBDA, "\n")
cat(QATF1$MSE, "at lambda : ", QATF1$LAMBDA, "\n")
cat(QATF2$MSE, "at lambda : ", QATF2$LAMBDA, "\n")


cum_data <- rbind(cum_data, data.frame(n = 500, Scenario = 1,
                                       tau=0.5, QATF1 = QATF1$MSE, QATF2 = QATF2$MSE,
                                       ATF1 = ATF1$MSE, ATF2 = ATF2$MSE))




## prepare the inputs X, and the true function
x <- seq(1, n, 1)/n



# get scenerio 4
# half <- n/2
# x[1:half] <- 3*x[1:half]
# x <- replace(x, seq(half+1, n, 1), 3*n - 3*tail(x,half))
# x <- x/n

# get scenerio 5
# x<- cos(6*pi*x/n)


# get the simulated heterogeneously-smooth data
x_j <- x # we could add permutation for different j.


# scenerio 6
# get_g_j <- function(x_j, j){
#   x_j <- 0
#   g_0 <- (0.1)*(j/10)
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j**(0.1)*(j/10)
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
# sce1 <- rnorm(n, 0, 1)
# sce2 <- rcauchy(n, 0, 1)
# e <- (x**.5)  # i**.5 / n**.5
# te <- rt(n, 2)
# sce3 <- e*te
# sce4 <- rt(n,3)
# sce5 <- rcauchy(n, 0, 1)
# y <- y_star + sce1

# get scenerio 6 error
e <- seq(1, n, 1)/n
e[1:n/4] <- (0.25*(e[1:n/4])**0.5 + 1.375)/3
e[(n/4+1):n] <- (7*(e[(n/4+1):n])**0.5 - 2)/3
sce6 <-rt(n,2)
sce6_9q <- qt(tau, 2)
y<-y_star+e*sce6

y_star <- y_star+e*sce6_9q
plot(y_star, type="l", col="black", ylab='0.9 quantile')



ATF2 <- get_mse(2)
ATF3 <- get_mse(3)
QATF2 <- get_mse_q(2)
QATF3 <- get_mse_q(3)

# cat(ATF2$MSE, "at lambda : ", ATF2$LAMBDA, "\n")
# cat(ATF3$MSE, "at lambda : ", ATF3$LAMBDA, "\n")
cat(QATF2$MSE, "at lambda : ", QATF2$LAMBDA, "\n")
cat(QATF3$MSE, "at lambda : ", QATF3$LAMBDA, "\n")


cum_data <- rbind(cum_data, data.frame(n = n, Scenario = 1,
                                       tau=tau, QATF2 = QATF2$MSE, QATF3 = QATF3$MSE,
                                       ATF2 = ATF2$MSE, ATF3 = ATF3$MSE))


# write data to a sample.csv file
write.csv(cum_data, file = "MSEs.csv")

cum_data_ten <- cum_data
cum_data_ten$QATF2 <- cum_data$QATF2 * 10
cum_data_ten$QATF3 <- cum_data$QATF3 * 10
cum_data_ten$ATF2 <- cum_data$ATF2 * 10
cum_data_ten$ATF3 <- cum_data$ATF3 * 10
write.csv(cum_data_ten, file = "10MSEs.csv")

#> Using same lambda for all quantiles
plot(y, type="l", col="black") # this is the y with errors
plot(y_star, type="l", col="black") # this is the y_star, without errors


# lines(trend_hat~x, col="red")
# lines(q_trend_hat~x, col="blue")
# lines(trend[,1]~x, col="red")
# lines(trend[,2]~x, col="blue")




