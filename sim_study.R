# QATF simulation study
# sim_study.R

source("scenarios.R")
source("algorithms.R")

# This is a wrapper function, 
# calling the appropriate scenario
# the appropriate algorithm(s)
# and simulating specified times! 
# The output is designed for building the data frame
run_custom_sce_simulations <- function(n, d, tau, sce, simulations = 1) {
  # This function only exists for us to make constructing the data frame quicker
  # Wrapper function is set to run with fixed k = 0 and 1, 
  # and append results to formatted data frame.
  
  ATF1_MSE <- 0
  ATF0_MSE <- 0
  QS_MSE <- 0
  QATF1_MSE <- 0
  QATF0_MSE <- 0
  for (i in 1:simulations) {
    if      (sce == 0) {vals <- scenario0(n, d, tau)}
    else if (sce == 1) {vals <- scenario1(n, d, tau)}
    else if (sce == 2) {vals <- scenario2(n, d, tau)}
    else if (sce == 3) {vals <- scenario3(n, d, tau)}
    else if (sce == 4) {vals <- scenario4(n, d, tau)}
    else if (sce == 5) {vals <- scenario5(n, d, tau)}
    else if (sce == 6) {vals <- scenario6(n, d, tau)}
    else if (sce == 7) {vals <- scenario6(n, d, tau)}
    else {stop("Only 0-7 scenarios at the time of this functions' construction")}
    
    # Currently preserving the full output in case we want cool plots
    # Also, only run get_mse for tau = 0.5 and not scenario 3
    if (tau == 0.5 && sce != 3) {
      ATF1 <- get_mse(vals[[1]], vals[[2]], vals[[3]], n, d, 1, prints = FALSE)
      ATF0 <- get_mse(vals[[1]], vals[[2]], vals[[3]], n, d, 0, prints = FALSE)
      
      ATF1_MSE <- ATF1_MSE + ATF1$MSE
      ATF0_MSE <- ATF0_MSE + ATF0$MSE
      
      cat("mse for ATF1 was ", ATF1$MSE, "at lambda = ", ATF1$LAMBDA, "\n")
      cat("mse for ATF0 was ", ATF0$MSE, "at lambda = ", ATF0$LAMBDA, "\n")
    } else {
      cat("not running ATF1 or ATF0\n")
    }
    
    QS <- get_mse_s(vals[[1]], vals[[2]], vals[[3]], n, d, tau, prints = FALSE)
    QS_MSE <- QS_MSE + QS$MSE
    cat("mse for QS was ", QS$MSE, "at lambda = ", QS$LAMBDA, "\n")
    
    QATF1 <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], n, d, tau, 1, prints = FALSE)
    QATF0 <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], n, d, tau, 0, prints = FALSE)
    QATF1_MSE <- QATF1_MSE + QATF1$MSE
    QATF0_MSE <- QATF0_MSE + QATF0$MSE
    
    cat("mse for QATF1 was ", QATF1$MSE, "at lambda = ", QATF1$LAMBDA, "\n")
    cat("mse for QATF0 was ", QATF0$MSE, "at lambda = ", QATF0$LAMBDA, "\n")
    if (simulations != 1) {cat("finished simulation ", i, "\n")}
  }
  
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
                    QATF1 = format(QATF1_MSE, scientific = FALSE, digits = 6),
                    QATF0 = format(QATF0_MSE, scientific = FALSE, digits = 6),
                    QS    = format(QS_MSE   , scientific = FALSE, digits = 6),
                    ATF1  = format(ATF1_MSE , scientific = FALSE, digits = 6),
                    ATF0  = format(ATF0_MSE , scientific = FALSE, digits = 6)))
  
}

# Construct Table
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
  # Scenario 2
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 2, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 2, simulations = 10))
  write.csv(cum_data, file = "scenario2_05.csv")
  # Scenario 3
  cum_data <- data.frame()
  cum_data <- rbind(cum_data, run_custom_sce_simulations( 500, 10, 0.5, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(1000, 10, 0.5, 3, simulations = 10))
  cum_data <- rbind(cum_data, run_custom_sce_simulations(2500, 10, 0.5, 3, simulations = 10))
  write.csv(cum_data, file = "scenario3_05.csv")
  # Scenario 4
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


############################################
# Next Code is for plots scenarios 5, 6, 7 # 
############################################

# Plot scenario 5
{
  vals <- scenario5(2000, 3, 0.5)
  # vals[[4]] is y_list, vals[[1]] is x_list
  # out[[2]] is trend_list, out[[3]] is permutation matrix
  
  par(mfrow = c(3, 3))
  lambda_atf <- get_mse(vals[[1]], vals[[2]], vals[[3]], 2000, 3, 2)
  out <- fit_atf(vals[[1]], vals[[2]], 2000, 3, 2, lambda_atf$LAMBDA)
  plot(vals[[4]][1, ][out[[3]][1, ]])
  lines(out[[2]][, 1][out[[3]][1, ]], col = "red")
  plot(vals[[4]][2, ][out[[3]][2, ]])
  lines(out[[2]][, 2][out[[3]][2, ]], col = "red")
  plot(vals[[4]][3, ][out[[3]][3, ]])
  lines(out[[2]][, 3][out[[3]][3, ]], col = "red")
  
  lambda_qass <- get_mse_s(vals[[1]], vals[[2]], vals[[3]], 2000, 3, 0.5)
  outs <- fit_qass(vals[[1]], vals[[2]], 2000, 3, 0.5, lambda_qass$LAMBDA)
  plot(vals[[4]][1, ][outs[[3]][1, ]])
  lines(outs[[2]][, 1][outs[[3]][1, ]], col = "blue")
  plot(vals[[4]][2, ][outs[[3]][2, ]])
  lines(outs[[2]][, 2][outs[[3]][2, ]], col = "blue")
  plot(vals[[4]][3, ][outs[[3]][3, ]])
  lines(outs[[2]][, 3][outs[[3]][3, ]], col = "blue")
  
  lambda_qatf <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], 2000, 3, 0.5, 2)
  outq <- fit_qatf(vals[[1]], vals[[2]], 2000, 3, 0.5, 2, lambda_qatf$LAMBDA)
  plot(vals[[4]][1, ][outq[[3]][1, ]])
  lines(outq[[2]][, 1][outq[[3]][1, ]], col = "purple")
  plot(vals[[4]][2, ][outq[[3]][2, ]])
  lines(outq[[2]][, 2][outq[[3]][2, ]], col = "purple")
  plot(vals[[4]][3, ][outq[[3]][3, ]])
  lines(outq[[2]][, 3][outq[[3]][3, ]], col = "purple")
}

# Plot scenario 6
{
  vals <- scenario6(1000, 4, 0.5)
  # vals[[4]] is y_list, vals[[1]] is x_list
  # out[[2]] is trend_list, out[[3]] is permutation matrix
  
  par(mfrow = c(1, 4))
  # lambda_atf <- get_mse(vals[[1]], vals[[2]], vals[[3]], 1000, 4, 2)
  # out <- fit_atf(vals[[1]], vals[[2]], 1000, 4, 2, lambda_atf$LAMBDA)
  # plot(vals[[4]][1, ][out[[3]][1, ]])
  # lines(out[[2]][, 1][out[[3]][1, ]], col = "red")
  # plot(vals[[4]][2, ][out[[3]][2, ]])
  # lines(out[[2]][, 2][out[[3]][2, ]], col = "red")
  # plot(vals[[4]][3, ][out[[3]][3, ]])
  # lines(out[[2]][, 3][out[[3]][3, ]], col = "red")
  # plot(vals[[4]][4, ][out[[3]][4, ]])
  # lines(out[[2]][, 4][out[[3]][4, ]], col = "red")
  #
  # lambda_qass <- get_mse_s(vals[[1]], vals[[2]], vals[[3]], 1000, 4, 0.5)
  # outs <- fit_qass(vals[[1]], vals[[2]], 1000, 4, 0.5, lambda_qass$LAMBDA)
  # plot(vals[[4]][1, ][outs[[3]][1, ]])
  # lines(outs[[2]][, 1][outs[[3]][1, ]], col = "blue")
  # plot(vals[[4]][2, ][outs[[3]][2, ]])
  # lines(outs[[2]][, 2][outs[[3]][2, ]], col = "blue")
  # plot(vals[[4]][3, ][outs[[3]][3, ]])
  # lines(outs[[2]][, 3][outs[[3]][3, ]], col = "blue")
  # plot(vals[[4]][4, ][outs[[3]][4, ]])
  # lines(outs[[2]][, 4][outs[[3]][4, ]], col = "blue")
  
  lambda_qatf <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], 1000, 4, 0.5, 2)
  outq <- fit_qatf(vals[[1]], vals[[2]], 1000, 4, 0.5, 2, lambda_qatf$LAMBDA)
  # Set the layout to have 4 plots in a 1x4 grid
  par(mfrow = c(1, 4), mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
  
  # Plot 1
  plot(vals[[4]][1, ][outq[[3]][1, ]], xaxt = 'n', yaxt = 'n', main = "", xlab = "", ylab = "")
  lines(outq[[2]][, 1][outq[[3]][1, ]], col = "blue", lwd = 2)
  
  # Plot 2
  plot(vals[[4]][2, ][outq[[3]][2, ]], xaxt = 'n', yaxt = 'n', main = "", xlab = "", ylab = "")
  lines(outq[[2]][, 2][outq[[3]][2, ]], col = "blue", lwd = 2)
  
  # Plot 3
  plot(vals[[4]][3, ][outq[[3]][3, ]], xaxt = 'n', yaxt = 'n', main = "", xlab = "", ylab = "")
  lines(outq[[2]][, 3][outq[[3]][3, ]], col = "blue", lwd = 2)
  
  # Plot 4
  plot(vals[[4]][4, ][outq[[3]][4, ]], xaxt = 'n', yaxt = 'n', main = "", xlab = "", ylab = "")
  lines(outq[[2]][, 4][outq[[3]][4, ]], col = "blue", lwd = 2)
  
  
}

# Plot scenario 7
{
  vals <- scenario7(2500, 2, 0.5)
  # vals[[4]] is y_list, vals[[1]] is x_list
  # out[[2]] is trend_list, out[[3]] is permutation matrix
  
  lambda_qatf <- get_mse_q(vals[[1]], vals[[2]], vals[[3]], 2500, 2, 0.5, 2)
  outq <- fit_qatf(vals[[1]], vals[[2]], 2500, 2, 0.5, 2, lambda_qatf$LAMBDA)
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

scenario5(500, 3, 0.5, TRUE)
scenario6(500, 4, 0.5, TRUE)
scenario7(500, 2, 0.5, TRUE)

# Run this in R terminal to remove all generated plots
# system("rm sce[0-9]_plot*")


