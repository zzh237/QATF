# QATF simulation study
# sim_study.R

source("scenarios.R")
source("algorithms.R")

# This is a wrapper function, 
# calling the appropriate scenario
# the appropriate algorithm(s)
# and simulating specified times! 
# The output is designed for building the data frame

# method can be "min_lambda" (default) or "optimal" 
run_custom_sce_simulations <- function(n, d, tau, sce, simulations = 1, nlambdas = 50, 
                                       method = "min_lambda", save = FALSE) {
  if (!(method == "min_lambda" || method == "optimal")) {
    stop("method can be \"min_lambda\" (default) or \"optimal\"")
  }
  
  # This function only exists for us to make constructing the data frame quicker
  # Wrapper function is set to run with fixed k = 0 and 1, 
  # and append results to formatted data frame.
  
  ATF1_mat <- matrix(data = NA, nrow = simulations, ncol = nlambdas)
  ATF0_mat <- matrix(data = NA, nrow = simulations, ncol = nlambdas)
  QS_mat <- matrix(data = NA, nrow = simulations, ncol = nlambdas)
  QATF1_mat <- matrix(data = NA, nrow = simulations, ncol = nlambdas)
  QATF0_mat <- matrix(data = NA, nrow = simulations, ncol = nlambdas)
  
  ATF1_MSE <- 0
  ATF0_MSE <- 0
  QS_MSE <- 0
  QATF1_MSE <- 0
  QATF0_MSE <- 0
  
  # TODO: make sure to save the output!
  for (i in 1:simulations) {
    if      (sce == 0) {vals <- scenario0(n, d, tau)}
    else if (sce == 1) {vals <- scenario1(n, d, tau)}
    else if (sce == 2) {vals <- scenario2(n, d, tau)}
    else if (sce == 3) {vals <- scenario3(n, d, tau)}
    else if (sce == 4) {vals <- scenario4(n, d, tau)}
    else if (sce == 5) {vals <- scenario5(n, d, tau)}
    else if (sce == 6) {vals <- scenario6(n, d, tau)}
    else {stop("Only 0-6 scenarios at the time of this functions' construction")}
    
    if (save) { saveRDS(vals, file = paste0("sce_", sce, "_simulation_", i, ".txt"))}
    
    ### Additive Trend Filtering ###
    # Only run get_mse for tau = 0.5 and not scenario 3
    if (tau == 0.5 && sce != 3) {
      ATF1 <- get_mse(vals[[1]], vals[[2]], vals[[3]], n, d, 1, prints = FALSE)
      ATF0 <- get_mse(vals[[1]], vals[[2]], vals[[3]], n, d, 0, prints = FALSE)
      if (length(ATF1$MSES) != nlambdas) {stop("Inconsistent number of lambdas. Adjust nlambdas parameter (default 50)")}
      if (length(ATF0$MSES) != nlambdas) {stop("Inconsistent number of lambdas. Adjust nlambdas parameter (default 50)")}
      
      ATF1_mat[i, ] <- ATF1$MSES
      ATF0_mat[i, ] <- ATF0$MSES

      ATF1_MSE <- ATF1_MSE + min(ATF1$MSES)
      ATF0_MSE <- ATF0_MSE + min(ATF0$MSES)
      if (method == "optimal") {
        cat("For simulation ", i, " The best mse for ATF1 was ", min(ATF1$MSES), 
            " at lambda = ", ATF1$LAMBDAS[which.min(ATF1$MSES)], "\n", sep = "")
        cat("For simulation ", i, " The best mse for ATF0 was ", min(ATF0$MSES), 
            " at lambda = ", ATF0$LAMBDAS[which.min(ATF0$MSES)], "\n", sep = "")
      }
    } else {
      if (method == "optimal") {
        cat("not running ATF1 or ATF0\n")
      }
    }
    
    
    ### Quantile Smoothing Splines ###
    QS <- get_mse_s(vals[[1]], vals[[2]], vals[[3]], n, d, tau, prints = FALSE)
    if (length(QS$MSES) != nlambdas) {stop("Inconsistent number of lambdas. Adjust nlambdas parameter (default 50)")}
    QS_mat[i, ] <- QS$MSES
    
    QS_MSE <- QS_MSE + min(QS$MSES)
    if (method == "optimal") {
      cat("For simulation ", i, " The best mse for QS was ", min(QS$MSES), 
          " at lambda = ", QS$LAMBDAS[which.min(QS$MSES)], "\n", sep = "")    
    }
    
    
    ### Quantile Additive Trend Filtering ###
    QATF1 <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], n, d, tau, 1, prints = FALSE)
    QATF0 <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], n, d, tau, 0, prints = FALSE)
    if (length(QATF1$MSES) != nlambdas) {stop("Inconsistent number of lambdas. Adjust nlambdas parameter (default 50)")}
    if (length(QATF0$MSES) != nlambdas) {stop("Inconsistent number of lambdas. Adjust nlambdas parameter (default 50)")}
    
    QATF1_mat[i, ] <- QATF1$MSES
    QATF0_mat[i, ] <- QATF0$MSES
    
    QATF1_MSE <- QATF1_MSE + min(QATF1$MSES)
    QATF0_MSE <- QATF0_MSE + min(QATF1$MSES)
    if (method == "optimal") {
      cat("For simulation ", i, " The best mse for QATF1 was ", min(QATF1$MSES), 
          " at lambda = ", QATF1$LAMBDAS[which.min(QATF1$MSES)], "\n", sep = "")
      cat("For simulation ", i, " The best mse for QATF0 was ", min(QATF0$MSES), 
          " at lambda = ", QATF0$LAMBDAS[which.min(QATF0$MSES)], "\n", sep = "")
    }
    
    if (simulations != 1) {cat("finished simulation ", i, "\n")}
  }
  

  ### Optimal Handling ###
  if (method == "optimal") {
    if (tau == 0.5 && sce != 3) {
      ATF1_MSE <- ATF1_MSE / simulations
      ATF0_MSE <- ATF0_MSE / simulations
    } else {
      ATF1_MSE <- NA
      ATF0_MSE <- NA
    }
    QS_MSE <- QS_MSE / simulations
    QATF1_MSE <- QATF1_MSE / simulations
    QATF0_MSE <- QATF0_MSE / simulations
    
    return(data.frame(n = n,
                      Scenario = sce,
                      d = d,
                      tau = tau,
                      Simulations = simulations,
                      Method = method, 
                      QATF1 = format(QATF1_MSE, scientific = FALSE, digits = 6),
                      QATF0 = format(QATF0_MSE, scientific = FALSE, digits = 6),
                      QS    = format(QS_MSE   , scientific = FALSE, digits = 6),
                      ATF1  = format(ATF1_MSE , scientific = FALSE, digits = 6),
                      ATF0  = format(ATF0_MSE , scientific = FALSE, digits = 6)))
  }
  
  ### Min Lambda Handling ### 
  else {
    if (tau == 0.5 && sce != 3) {
      ATF1_MSE <- min(colMeans(ATF1_mat))
      ATF0_MSE <- min(colMeans(ATF0_mat))
      cat("Min lambda across ", simulations, " simulations for ATF1 was ", ATF1$LAMBDAS[which.min(colMeans(ATF1_mat))], "\n", sep = "")
      cat("Min lambda across ", simulations, " simulations for ATF0 was ", ATF0$LAMBDAS[which.min(colMeans(ATF0_mat))], "\n", sep = "")
      } else {
      ATF1_MSE <- NA
      ATF0_MSE <- NA
      cat("Did not run ATF1 or ATF0\n")
      }
    
    QS_MSE <- min(colMeans(QS_mat))
    cat("Min lambda across ", simulations, " simulations for QS was ", QS$LAMBDAS[which.min(colMeans(QS_mat))], "\n", sep = "")
    
    
    QATF1_MSE <- min(colMeans(QATF1_mat))
    QATF0_MSE <- min(colMeans(QATF0_mat))
    cat("Min lambda across ", simulations, " simulations for QATF1 was ", QATF1$LAMBDAS[which.min(colMeans(QATF1_mat))], "\n", sep = "")
    cat("Min lambda across ", simulations, " simulations for QATF0 was ", QATF0$LAMBDAS[which.min(colMeans(QATF1_mat))], "\n", sep = "")
    
    return(data.frame(n = n,
                      Scenario = sce,
                      d = d,
                      tau = tau,
                      Simulations = simulations,
                      Method = method, 
                      QATF1 = format(QATF1_MSE, scientific = FALSE, digits = 6),
                      QATF0 = format(QATF0_MSE, scientific = FALSE, digits = 6),
                      QS    = format(QS_MSE   , scientific = FALSE, digits = 6),
                      ATF1  = format(ATF1_MSE , scientific = FALSE, digits = 6),
                      ATF0  = format(ATF0_MSE , scientific = FALSE, digits = 6)))
  }

  
}

run_qgam_sce_simulations <- function(n, d, tau, sce, simulations = 1, nlambdas = 50, 
                                       method = "min_lambda", save = FALSE) {
  if (!(method == "min_lambda" || method == "optimal")) {
    stop("method can be \"min_lambda\" (default) or \"optimal\"")
  }
  
  QS_qgam_mat <- matrix(data = NA, nrow = simulations, ncol = nlambdas)
  QS_qgam_MSE <- 0
  
  QS_mat <- matrix(data = NA, nrow = simulations, ncol = nlambdas)
  QS_MSE <- 0
  
  # TODO: make sure to save the output!
  for (i in 1:simulations) {
    if      (sce == 0) {vals <- scenario0(n, d, tau)}
    else if (sce == 1) {vals <- scenario1(n, d, tau)}
    else if (sce == 2) {vals <- scenario2(n, d, tau)}
    else if (sce == 3) {vals <- scenario3(n, d, tau)}
    else if (sce == 4) {vals <- scenario4(n, d, tau)}
    else if (sce == 5) {vals <- scenario5(n, d, tau)}
    else if (sce == 6) {vals <- scenario6(n, d, tau)}
    else {stop("Only 0-6 scenarios at the time of this functions' construction")}
    
    if (save) { saveRDS(vals, file = paste0("sce_", sce, "_simulation_", i, ".txt"))}

    ### Quantile Smoothing Splines ###
    QS <- get_mse_s(vals[[1]], vals[[2]], vals[[3]], n, d, tau, prints = FALSE)
    QS_qgam <- get_mse_s_qgam(vals[[1]], vals[[2]], vals[[3]], n, d, tau, prints = FALSE)
    if (length(QS$MSES) != nlambdas) {stop("Inconsistent number of lambdas. Adjust nlambdas parameter (default 50)")}
    if (length(QS_qgam$MSES) != nlambdas) {stop("Inconsistent number of lambdas. Adjust nlambdas parameter (default 50)")}
    QS_mat[i, ] <- QS$MSES
    QS_qgam_mat[i, ] <- QS_qgam$MSES
    
    QS_MSE <- QS_MSE + min(QS$MSES)
    QS_qgam_MSE <- QS_qgam_MSE + min(QS_qgam$MSES)
    if (method == "optimal") {
      cat("For simulation ", i, " The best mse for QS was ", min(QS$MSES), 
          " at lambda = ", QS$LAMBDAS[which.min(QS$MSES)], "\n", sep = "")    
      cat("For simulation ", i, " The best mse for QS_qgam was ", min(QS_qgam$MSES), 
          " at lambda = ", QS_qgam$LAMBDAS[which.min(QS_qgam$MSES)], "\n", sep = "")  
    }
    if (simulations != 1) {cat("finished simulation ", i, "\n")}
  }
  
  
  ### Optimal Handling ###
  if (method == "optimal") {
    QS_MSE <- QS_MSE / simulations
    QS_qgam_MSE <- QS_qgam_MSE / simulations
    
    return(data.frame(n = n,
                      Scenario = sce,
                      d = d,
                      tau = tau,
                      Simulations = simulations,
                      Method = method, 
                      QS = format(QS_MSE   , scientific = FALSE, digits = 6),
                      QS_qgam = format(QS_qgam_MSE, scientific = FALSE, digits = 6)))
  }
  
  ### Min Lambda Handling ### 
  else {
    QS_MSE <- min(colMeans(QS_mat))
    QS_qgam_MSE <- min(colMeans(QS_qgam_mat))
    cat("Min lambda across ", simulations, " simulations for QS was ", QS$LAMBDAS[which.min(colMeans(QS_mat))], "\n", sep = "")
    cat("Min lambda across ", simulations, " simulations for QS_qgam was ", QS_qgam$LAMBDAS[which.min(colMeans(QS_qgam_mat))], "\n", sep = "")
    
    
    return(data.frame(n = n,
                      Scenario = sce,
                      d = d,
                      tau = tau,
                      Simulations = simulations,
                      Method = QS_qgam$LAMBDAS[which.min(colMeans(QS_qgam_mat))], 
                      QS = format(QS_MSE   , scientific = FALSE, digits = 6),
                      QS_qgam = format(QS_qgam_MSE, scientific = FALSE, digits = 6)))
  }
  
  
}

# run_qgam_sce_simulations(500, 10, 0.5, 1, simulations = 1)


# Construct Table: Main paper experiments
# Scenario 1, 2, 3, 4, 0, with tau 0.5
# Scenario 4 with tau 0.2, 0.8
{
  # Scenario 0
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 0, simulations = 10))
  write.csv(cum_data, file = "scenario0_05.csv")
  # Scenario 1
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 1, simulations = 10))
  write.csv(cum_data, file = "scenario1_05.csv")
  # # Scenario 2
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 2, simulations = 10))
  write.csv(cum_data, file = "scenario2_05.csv")
  # # Scenario 3
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 3, simulations = 10))
  write.csv(cum_data, file = "scenario3_05.csv")
  # # Scenario 4
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 4, simulations = 10))
  write.csv(cum_data, file = "scenario4_05.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.2, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.2, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.2, 4, simulations = 10))
  write.csv(cum_data, file = "scenario4_02.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.8, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.8, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.8, 4, simulations = 10))
  write.csv(cum_data, file = "scenario4_08.csv")
}

# Construct Table: Supplementary Experiments 
# Scenario 5, 6 with tau 0.5
# Scenario 1, 2, 3, 0, 5 with tau 0.2, 0.8
{
  # Scenario 0
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.2, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.2, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.2, 0, simulations = 10))
  write.csv(cum_data, file = "scenario0_02.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.8, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.8, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.8, 0, simulations = 10))
  write.csv(cum_data, file = "scenario0_08.csv")
  # Scenario 1
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.2, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.2, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.2, 1, simulations = 10))
  write.csv(cum_data, file = "scenario1_02.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.8, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.8, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.8, 1, simulations = 10))
  write.csv(cum_data, file = "scenario1_08.csv")
  # Scenario 2
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.2, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.2, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.2, 2, simulations = 10))
  write.csv(cum_data, file = "scenario2_02.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.8, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.8, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.8, 2, simulations = 10))
  write.csv(cum_data, file = "scenario2_08.csv")
  # Scenario 3
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.2, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.2, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.2, 3, simulations = 10))
  write.csv(cum_data, file = "scenario3_02.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.8, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.8, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.8, 3, simulations = 10))
  write.csv(cum_data, file = "scenario3_08.csv")
  # Scenario 5
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 5, simulations = 10))
  write.csv(cum_data, file = "scenario5_05.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.2, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.2, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.2, 5, simulations = 10))
  write.csv(cum_data, file = "scenario5_02.csv")
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.8, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.8, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.8, 5, simulations = 10))
  write.csv(cum_data, file = "scenario5_08.csv")
  # Scenario 6
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 6, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 6, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 6, simulations = 10))
  write.csv(cum_data, file = "scenario6_05.csv")
}

# Construct Table:qgam method only
# Scenario 1-5, tau 0.2, 0.8, 0.5
# Scenario 6, tau 0.5 only
{
  # Scenario 0
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.2, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.2, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.2, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.5, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.5, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.5, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.8, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.8, 0, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.8, 0, simulations = 10))
  write.csv(cum_data, file = "scenario0_qgam.csv")
  # Scenario 1
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.2, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.2, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.2, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.5, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.5, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.5, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.8, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.8, 1, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.8, 1, simulations = 10))
  write.csv(cum_data, file = "scenario1_qgam.csv")
  # Scenario 1
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.2, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.2, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.2, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.5, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.5, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.5, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.8, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.8, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.8, 2, simulations = 10))
  write.csv(cum_data, file = "scenario2_qgam.csv")
  # Scenario 3
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.2, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.2, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.2, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.5, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.5, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.5, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.8, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.8, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.8, 3, simulations = 10))
  write.csv(cum_data, file = "scenario3_qgam.csv")
  # Scenario 4
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.2, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.2, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.2, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.5, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.5, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.5, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.8, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.8, 4, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.8, 4, simulations = 10))
  write.csv(cum_data, file = "scenario4_qgam.csv")
  # Scenario 5
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.2, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.2, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.2, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.5, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.5, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.5, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.8, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.8, 5, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.8, 5, simulations = 10))
  write.csv(cum_data, file = "scenario5_qgam.csv")
  # Scenario 6
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_qgam_sce_simulations( 500, 10, 0.5, 6, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(1000, 10, 0.5, 6, simulations = 10))
  cum_data <- rbind(cum_data, run_qgam_sce_simulations(2500, 10, 0.5, 6, simulations = 10))
  write.csv(cum_data, file = "scenario6_qgam.csv")
}


############################################
# Next Code is for plots scenarios 5, 6, 7 # 
############################################

# Plot scenario plot 1, all competitors
{
  vals <- scenarioplot1(500, 3, 0.5)
  # vals[[4]] is y_list, vals[[1]] is x_list
  # out[[2]] is trend_list, out[[3]] is permutation matrix
  
  par(mfrow = c(3, 3))
  lambda_res <- get_mse(vals[[1]], vals[[2]], vals[[3]], 500, 3, 2)
  lambda <- lambda_res$LAMBDAS[which.min(lambda_res$MSES)]
  out <- fit_atf(vals[[1]], vals[[2]], 500, 3, 2, lambda)
  plot(vals[[4]][1, ][out[[3]][1, ]])
  lines(out[[2]][, 1][out[[3]][1, ]], col = "red")
  plot(vals[[4]][2, ][out[[3]][2, ]])
  lines(out[[2]][, 2][out[[3]][2, ]], col = "red")
  plot(vals[[4]][3, ][out[[3]][3, ]])
  lines(out[[2]][, 3][out[[3]][3, ]], col = "red")
  
  lambda_sres <- get_mse_s(vals[[1]], vals[[2]], vals[[3]], 500, 3, 0.5)
  lambdas <- lambda_sres$LAMBDAS[which.min(lambda_sres$MSES)]
  print(lambdas)
  outs <- fit_qass(vals[[1]], vals[[2]], 500, 3, 0.5, lambdas)
  plot(vals[[4]][1, ][outs[[3]][1, ]])
  lines(outs[[2]][, 1][outs[[3]][1, ]], col = "blue")
  plot(vals[[4]][2, ][outs[[3]][2, ]])
  lines(outs[[2]][, 2][outs[[3]][2, ]], col = "blue")
  plot(vals[[4]][3, ][outs[[3]][3, ]])
  lines(outs[[2]][, 3][outs[[3]][3, ]], col = "blue")
  
  lambda_qres <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], 500, 3, 0.5, 2)
  lambdaq <- lambda_qres$LAMBDAS[which.min(lambda_qres$MSES)]
  print(lambdaq)
  outq <- fit_qatf(vals[[1]], vals[[2]], 500, 3, 0.5, 2, lambdaq)
  plot(vals[[4]][1, ][outq[[3]][1, ]])
  lines(outq[[2]][, 1][outq[[3]][1, ]], col = "purple")
  plot(vals[[4]][2, ][outq[[3]][2, ]])
  lines(outq[[2]][, 2][outq[[3]][2, ]], col = "purple")
  plot(vals[[4]][3, ][outq[[3]][3, ]])
  lines(outq[[2]][, 3][outq[[3]][3, ]], col = "purple")
}

# Plot scenario plot 2, QATF only
{
  vals <- scenarioplot2(1000, 4, 0.5)
  source("algorithms.R")
  # vals[[4]] is y_list, vals[[1]] is x_list
  # out[[2]] is trend_list, out[[3]] is permutation matrix
  
  par(mfrow = c(1, 4))
  
  # Get the optimal lambda
  lambda_qres <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], 1000, 4, 0.5, 2)
  lambdaq <- lambda_qres$LAMBDAS[which.min(lambda_qres$MSES)]
  print(lambdaq)
  outq <- fit_qatf(vals[[1]], vals[[2]], 1000, 4, 0.5, 2, lambdaq)
  
  # Set up the plot layout (4 plots in 1 row)
  par(mfrow = c(1, 4), mar = c(5, 4, 4, 2) + 0.1) # Reset the margins for better spacing
  
  # Loop through each component (1 to 4) and plot with appropriate axis ticks and titles
  for (i in 1:4) {
    plot(vals[[4]][i, ][outq[[3]][i, ]], 
         xaxt = 's', yaxt = 's', # Auto-generate axis ticks
         main = paste("Component", i), # Set the main title dynamically
         xlab = "", ylab = "") # Optionally set axis labels
    lines(outq[[2]][, i][outq[[3]][i, ]], col = "blue", lwd = 2) # Add trend line
  }
}

# Plot scenario plot 3
{
  vals <- scenarioplot3(2500, 2, 0.5)
  # vals[[4]] is y_list, vals[[1]] is x_list
  # out[[2]] is trend_list, out[[3]] is permutation matrix
  
  lambda_qres <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], 2500, 2, 0.5, 2)
  lambdaq <- lambda_qres$LAMBDAS[which.min(lambda_qres$MSES)]
  outq <- fit_qatf(vals[[1]], vals[[2]], 2500, 2, 0.5, 2, lambdaq)
  plot <- plot_ly(x = ~vals[[1]][1, ], y = ~vals[[1]][2, ], z = ~vals[[3]], 
                  type = "scatter3d", mode = "markers", 
                  marker = list(size = 3, color = ~vals[[3]], colorscale = 'Viridis')) %>%
    layout(scene = list(
      xaxis = list(title = "x1"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "true y")
    ))
  
  print(plot)
  
  plot <- plot_ly(x = ~vals[[1]][1, ], y = ~vals[[1]][2, ], z = ~vals[[2]], 
                  type = "scatter3d", mode = "markers", 
                  marker = list(size = 3, color = ~vals[[2]], colorscale = 'Viridis')) %>%
    layout(scene = list(
      xaxis = list(title = "x1"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "noisy data")
    ))
  
  print(plot)
  
  plot <- plot_ly(x = ~vals[[1]][1, ], y = ~vals[[1]][2, ], z = ~outq[[1]], 
                  type = "scatter3d", mode = "markers", 
                  marker = list(size = 3, color = ~outq[[1]], colorscale = 'Viridis')) %>%
    layout(scene = list(
      xaxis = list(title = "x1"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "fit")
    ))
  
  print(plot)
  
}




############################################
# Next code is to generate component plots #
############################################

scenario0(500, 10, 0.5, TRUE)
scenario1(500, 10, 0.5, TRUE)
scenario2(500, 10, 0.5, TRUE)
scenario3(500, 10, 0.5, TRUE)
scenario4(500, 10, 0.5, TRUE)
scenario5(500, 10, 0.5, TRUE)
scenario6(500, 10, 0.5, TRUE)

scenarioplot1(500, 3, 0.5, TRUE)
scenarioplot2(500, 4, 0.5, TRUE)
scenarioplot3(500, 2, 0.5, TRUE)

# Run this in R terminal to remove all generated plots
# system("rm sce[0-9]_plot*")


