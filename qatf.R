# install.packages("devtools")
library(devtools)
install_github("zzh237/detrendr")
install_github("glmgen/genlasso")
# not sure if this next one is necessary - Zhi can you check? 
# install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
install.packages("fields")

library(detrendr)
# trace(get_model, edit=TRUE)
library(glmgen)
library(fields)


# rm(list = ls())
prior_results <- read.csv("MSEs.csv")

MSE <- function(a, b){
  len = length(a)
  return(norm(a - b, type="2")**2/len) 
}

# These functions require all inputs to be passed in
# Return an object with the best_mse, best_lambda, and best_fit
# optional parameters  allow for better control
get_mse <- function(y, y_star, n, d, k, alpha = 10**-6, max_t = 50, prints = TRUE){
  lambda_list <- 10**seq(8, -8, length.out=50)
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
get_mse_s <- function(y, y_star, n, d, tau, alpha = 10**-6, max_t = 50, prints = TRUE){
  lambda_list <- 10**seq(1, -15, length.out=50)
  # cubic splines of this form seem to prefer very small lambda values
  
  best_mse <- Inf 
  best_lambda <- Inf
  best_s_trend_hat <- rep(0, times = n)
  
  for (lambda in lambda_list) {
    s_trend_list <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
    # we initialize all component functions to be 0
    
    s_trend_hat_prev <- rep(0, times = n)
    t <- 1
    repeat {
      for (j in 1:d){
        # calculate jth partial residual using components not equal to j
        resp <- y - as.numeric(t(rowSums(s_trend_list[, -j])))
        s_trend_list[, j] <- qsreg(seq(1, n, 1)/n, resp, lam = lambda, alpha = tau)$fitted.values
      }
      s_trend_hat <- rowSums(s_trend_list)
      
      if (all((abs(s_trend_hat - s_trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      s_trend_hat_prev <- s_trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, s_trend_hat) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")}
    
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_s_trend_hat <- s_trend_hat
    }
  }
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_s_trend_hat))
}
get_mse_q <- function(y, y_star, n, d, tau, k, alpha = 10**-6, max_t = 50, prints = TRUE){
  lambda_list <- 10**seq(8, -8, length.out=50)
  
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

# Scenario functions wrap our scenarios
# construct plots withing
# and return the true quantile fit, y_star_q, and the data, y
# for tau = 0.5, y_star_q = y_star
scenario1 <- function(n, d, tau) {
  # Scenario 1
  # i <- 1:n
  # g_0(x) <- sin(2*pi/(x + 0.1)**(j/10))
  # x <- 1/n equally spaced
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i normal errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- seq(1, n, 1)/n
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Doppler-like
    g_0 <- sin(2*pi/(x + 0.1)**(j/10))
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type="2")/sqrt(n))
    y_j <- a_j*g_0 - a_j*b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Normal Errors
  y <- y_star + rnorm(n, 0, 1)
  y_star_q <- y_star + qnorm(tau, 0, 1)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, col = "black", pch = 19, cex = 0.5, ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario2 <- function(n, d, tau) {
  # Scenario 2
  # i <- 1:n
  # g_0(x) <- sin(2*pi/(x + 0.1)**(j/10))
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i cauchy errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  y_list <- data.frame()
  for (j in 1:d) {
    x_j <- sort(runif(n, 0, 1))
    # Doppler-like
    g_0 <- sin(2*pi/(x_j + 0.1)**(j/10))
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type="2")/sqrt(n))
    y_j <- a_j*g_0 - a_j*b_j
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
  plot(y, col = "black", pch = 19, cex = 0.5, ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario3 <- function(n, d, tau) {
  # Scenario 3
  # i <- 1:n
  # g_0(x) <- sin(2*pi/(x + 0.1)**(j/10))
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i heteroskedastic t(2) errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  y_list <- data.frame()
  for (j in 1:d) {
    x_j <- sort(runif(n, 0, 1))
    # Doppler-like
    g_0 <- sin(2*pi/(x_j + 0.1)**(j/10))
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type="2")/sqrt(n))
    y_j <- a_j*g_0 - a_j*b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Heteroskedastic T errors
  # function of index
  e <- (seq(1, n, 1)**.5)/(n**.5)  # i**.5 / n**.5
  y <- y_star + e * rt(n, 2)
  y_star_q <- y_star + e * qt(tau, 2)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, col = "black", pch = 19, cex = 0.5, ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario4 <- function(n, d, tau) {
  # Scenario 4
  # i <- 1:n
  # g_0(x) <- (x + 0.1)*(j/10)
  # x <- 3(i/n) for (1:n/2), 3(1 - i/n) for (n/2 + 1:n)
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i t(3) errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- seq(1, n, 1)/n
  x[1:(n/2)] <- 3*x[1:(n/2)]
  x[(n/2 + 1):n] <- 3*(1 - x[(n/2 + 1):n])
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Linear
    g_0 <- (x + 0.1)*(j/10)
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type="2")/sqrt(n))
    y_j <- a_j*g_0 - a_j*b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # T errors
  y <- y_star + rt(n, 3)
  y_star_q <- y_star + qt(tau, 3)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, col = "black", pch = 19, cex = 0.5, ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario5 <- function(n, d, tau) {
  # Scenario 5
  # i <- 1:n
  # g_0(x) <- (x + 0.1)*(j/10)
  # x <- cos(6*pi*(i/n))
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i cauchy errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x <- cos(6*pi*seq(1, n, 1)/n)
  
  y_list <- data.frame()
  for (j in 1:d) {
    # Linear 
    g_0 <- (x + 0.1)*(j/10)
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type="2")/sqrt(n))
    y_j <- a_j*g_0 - a_j*b_j
    y_list <- rbind(y_list, y_j)
  }
  
  y_list <- unname(y_list)
  y_star <- colSums(y_list)
  
  # Cauchy errors
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star, type = "l", col = "black", ylab = "true values")
  plot(y, col = "black", pch = 19, cex = 0.5, ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}
scenario6 <- function(n, d, tau) {
  # Scenario 1
  # i <- 1:n
  # x <- 1/n equally spaced
  # y <- vi*epsilon_i
  # vi complicated
  # epsilon_i t(2) errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  e <- seq(1, n, 1)/n
  e[1:n/4] <- (0.25*(e[1:n/4])**0.5 + 1.375)/3
  e[(n/4+1):n] <- (7*(e[(n/4+1):n])**0.5 - 2)/3
  y<-e*rt(n,2)
  y_star_q <- e*qt(tau, 2)
  
  par(mfrow = c(1, 2))
  # Plot the true signal and the data
  plot(y_star_q, type = "l", col = "black", ylab = paste(tau, "quantile"))
  plot(y, col = "black", pch = 19, cex = 0.5, ylab = "data")
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(y, y_star_q))
}


# This is a wrapper function, 
# calling the appropriate scenario
# the appropriate algorithm(s)
# and simulating! 
# The output is designed for building the data frame
run_custom_sce_simulations <- function(n, d, tau, sce, simulations = 1) {
  # This function only exists for us to make constructing the data frame quicker
  # Wrapper function is set to run with fixed k = 1 and 2, 
  # and append results to formatted data frame.
  
  ATF1_MSE <- 0
  ATF2_MSE <- 0
  QS_MSE <- 0
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
      cat("not running ATF1 or ATF2 since tau != 0.5\n")
    }
    
    QS <- get_mse_s(vals[[1]], vals[[2]], n, d, tau, prints = FALSE)
    QS_MSE <- QS_MSE + QS$MSE
    cat("mse for QS was ", QS$MSE, "at lambda = ", QS$LAMBDA, "\n")
    
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
  QS_MSE <- QS_MSE / simulations
  QATF1_MSE <- QATF1_MSE / simulations
  QATF2_MSE <- QATF2_MSE / simulations
  
  return(data.frame(n = n,
                    Scenario = sce,
                    tau = tau,
                    Simulations = simulations,
                    QATF1 = format(QATF1_MSE, scientific = FALSE, digits = 6),
                    QATF2 = format(QATF2_MSE, scientific = FALSE, digits = 6),
                    QS    = format(QS_MSE   , scientific = FALSE, digits = 6),
                    ATF1  = format(ATF1_MSE , scientific = FALSE, digits = 6),
                    ATF2  = format(ATF2_MSE , scientific = FALSE, digits = 6)))
  
}

# Test various scenarios
scenario1(500, 10, 0.5)
scenario2(500, 10, 0.5)
scenario3(500, 10, 0.5)
scenario4(500, 10, 0.5)
scenario5(500, 10, 0.5)
scenario6(500, 10, 0.9)
scenario6(500, 10, 0.1)

run_custom_sce_simulations( 500, 10, 0.5, 5, simulations = 1)


# Run these lines to build the data frame
# Adjust simulations to your liking
cum_data <- data.frame()
# Scenario 1
cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 1, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 1, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 1, simulations = 1))
# Scenario 2
cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 2, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 2, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 2, simulations = 1))
# Scenario 3
cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 3, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 3, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 3, simulations = 1))
# Scenario 4
cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 4, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 4, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 4, simulations = 1))
# Scenario 5
cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 5, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 5, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 5, simulations = 1))
# Scenario 6
cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.9, 6, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.9, 6, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.9, 6, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.1, 6, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.1, 6, simulations = 1))
cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.1, 6, simulations = 1))
write.csv(cum_data, file = "MSEswqs.csv")

# Test a single scenario, great for plots
vals <- scenario1(500, 10, 0.5)
ATF1 <- get_mse(vals[[1]], vals[[2]], 500, 10, 1)
ATF2 <- get_mse(vals[[1]], vals[[2]], 500, 10, 2)
QS <- get_mse_s(vals[[1]], vals[[2]], 500, 10, 0.5)
QATF1 <- get_mse_q(vals[[1]], vals[[2]], 500, 10, 0.5, 1)
QATF2 <- get_mse_q(vals[[1]], vals[[2]], 500, 10, 0.5, 2)

cat(ATF1$MSE, "at lambda : ", ATF1$LAMBDA, "\n")
cat(ATF2$MSE, "at lambda : ", ATF2$LAMBDA, "\n")
cat(QS$MSE, "at lambda : ", QS$LAMBDA, "\n")
cat(QATF1$MSE, "at lambda : ", QATF1$LAMBDA, "\n")
cat(QATF2$MSE, "at lambda : ", QATF2$LAMBDA, "\n")




# Extra plotting ideas

# lines(trend_hat~x, col="red")
# lines(q_trend_hat~x, col="blue")
# lines(trend[,1]~x, col="red")
# lines(trend[,2]~x, col="blue")
