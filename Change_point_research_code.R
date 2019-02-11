#All libraries used throughout the simulations
library(grDevices)
library(gridExtra)
library(grid)
library(gtable)
library(sandwich)
library(tseries)
library(quantmod)
library(fGarch)
library(aTSA)
library(crypto)
library(xts)


###################################
#Original attempt at simulating Brownian motion. This method doesn't work well numerically because the sin function oscillates too much. 
###################################

w_process <- function(N, T_, beta){
  
  Z <- rnorm(N+1)
  const <- sqrt(2)/pi
  t_val <- seq(0, T_, by = (1/N)/T_)
  t_val_plus <- t_val + 1
  J <- c(1:N)
  S <- sin(pi*(J%*%t(t_val)))
  S_plus <- sin(pi*(J%*%t(t_val_plus)))
  f <- sqrt(T_)*(Z[1]*t_val + const*colSums(diag(Z[2:(N+1)]/J)%*%S))
  f_plus <- sqrt(T_)*(Z[1]*t_val_plus + const*colSums(diag(Z[2:(N+1)]/J)%*%S_plus))
  
  final_vec <- abs(f_plus - f)/((t_val_plus)^beta)
  
  return(final_vec)
}
#######################################################
#W_process redevelopment
#######################################################

W_process <- function(times, M, beta){
  
  set.seed(M)
  max_vector <- vector()
  for(i in 1:times){
    du = .0001
    N = floor((M + 1)/du)
    W <- c(0, cumsum(rnorm(N)/sqrt(M+1)))
    u <- seq(0, M+1, by = du)
    test_stat <- abs(W[2:(N+1)] - W[1:N])/((u[1:N]+1)^beta)
    max_vector[i] <- max(test_stat)
  }

  return(max_vector)
}


##########################################
#Final method used for Brownian motion simulation
#times = number of times to run the simulation
#M = endpoint of the interval over which to generate the random quantities
#beta = value of beta for the test denominator
##########################################
W_process_mod <- function(times, M, beta){
  
  
  max_vector <- vector()
  for(i in 1:times){
    set.seed(i)
    du = .0001
    N = floor((M + 1)/du)
    #W <- c(0,rwiener(end = M+1, frequency = 1/du))
    W <- c(0, cumsum(rnorm(N)))*(sqrt((M+1)/N))
    #W <- c(0, cumsum(rnorm(N)/sqrt(N)))*(sqrt(M))
    u <- seq(0, M, by = du)
    start <- 1 + 1/du
    end <- (M+1)/du + 1
    end2 <- M/du + 1
    test_stat <- abs(W[start:end] - W[1:end2])
    final <- test_stat/((u+1)^beta)
    max_vector[i] <- max(final)
  }
  
  return(max_vector)
}



#Generates the vector of maximums for a given number of iterations and a given beta value
test_statistic = function(N, beta, endpoint){
  

  maximums <- vector()
  
  for(i in 1:(N+1)){
    value_vec <- w_process(N, endpoint, beta)
    maximums[i] <- max(value_vec)
  }
  
  return(maximums)
  
}


#############################################
#
#############################################
run_simulation <- function(N, Beta, endpoint){
  
  max_values <- W_process_mod(N, endpoint, Beta)
  
  cdf <- ecdf(max_values)
  
  #put these plots side by side
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/cdf", endpoint, "b", Beta, ".pdf", sep = "")
  pdf(file_name)

  plot(cdf, main = bquote("Empirical CDF When " ~ beta == .(Beta) ~ ", " ~ End == .(endpoint)))
  dev.off()
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/pdf", endpoint, "b", Beta, ".pdf", sep = "")
  pdf(file_name)
  plot(density(max_values), main = bquote("Density Estimate When " ~ beta == .(Beta) ~ ", " ~ End == .(endpoint)), xlab = bquote(x))
  dev.off()
  
  #obtain critical values
  
  cv <- unname(quantile(max_values, c(.90, .95, .99)))
  
  return(cv)
}



#############################################
#
#############################################
critical_value_table <- function(N, endpoint){
  
  
  cv1 <- run_simulation(N, 1, endpoint)
  cv2 <- run_simulation(N, 2, endpoint)
  cv3 <- run_simulation(N, 3, endpoint)
  cv4 <- run_simulation(N, 4, endpoint)
  cv5 <- run_simulation(N, 5, endpoint)
  cv6 <- run_simulation(N, 6, endpoint)
  cv7 <- run_simulation(N, 7, endpoint)
  cv8 <- run_simulation(N, 8, endpoint)
  
  critical_values <- rbind(cv1, cv2, cv3, cv4, cv5, cv6, cv7, cv8)
  #input into data frame with beta values
  df <- data.frame(c("1", "2", "3", "4", "5", "6", "7", "8"), critical_values)
  
  names(df) <- c("Beta", "90%", "95%", "99%" )
  
  g <- tableGrob(df, rows = NULL)
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, l = 1, r = ncol(g))
  grid.newpage()
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/cv", endpoint,".jpeg", sep = "")
  jpeg(file = file_name)
  grid.draw(g)
  dev.off()
  
}

############Upper Bound at 1###################

critical_value_table(5000, 1)

############Upper Bound at 2###################

critical_value_table(5000, 2)

############Upper Bound at 3###################

critical_value_table(5000, 3)

############Upper Bound at 4###################
critical_value_table(5000, 4)
############Upper Bound at 5###################
critical_value_table(5000, 5)

############Upper Bound at 5###################
critical_value_table(5000, 6)
critical_value_table(5000, 10)
critical_value_table(500, 20)
critical_value_table(500, 40)
critical_value_table(500, 80)
critical_value_table(1000, 100)
critical_value_table(200, 200)
critical_value_table(500, 400)
critical_value_table(500, 1000)




unname(quantile(W_process_mod(1000, 100, 3/4), c(.90, .95, .99)))
##########################################Analysis of Type 1 Error################################################

#####Choose M, h such that h < < M, assume that mean is 0
c_alphas_final_34 <- c(2.069211, 2.316345, 2.924225)
c_alphas_final_1 <- c(1.961971, 2.236345, 2.725577)
c_alphas_final_2 <- c(1.813968, 2.105014, 2.625778)
c_alphas_final_3 <- c(1.756301, 2.046612, 2.586148)
c_alphas_final_4 <- c(1.735013, 2.010628, 2.568760)


#############################################
#Calculates the empirical test size using a standard normal distribution
#M
#t
#N
#beta
#c_alpha
#############################################
calculate_alpha <- function(M, t, N, beta, c_alpha){
 
  max_vector <- vector() 
  h = floor(t^(1/2))
  for(i in 1:N){
    set.seed(i)
    training <- rnorm(M)
    window_vec <- vector()
    M_mean <- mean(training)
    test <- rnorm(t)
    diff_vector <- vector()
    sd_estimate <- sd(training)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+ 1/2})
      diff_vector[k] <- abs(M_mean - window_vec[k])/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }

  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}
#############################################
#
#############################################
sample_test_statistic_stop_time <- function(train, test, beta, c_alpha){
  
  x_bar <- mean(train)
  t <- length(test)
  window <- floor(t^(1/2))
  #window <- 2
  window_vec <- vector()
  diff_vector <- vector()
  
  stop_obs <- 0
  
  sd_estimate <- sd(train)
  
  for(k in 1:(t-window)){
    sum =0 
    
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    #g_T <- sd_estimate/(sqrt(window*t))
    if(abs(x_bar - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- k
      break
    }
  }
  
  
  return(stop_obs)
}

#############################################
#
#############################################
simulate_tau_sn <- function(M, T_, N, beta, c_alphas){
  sizes <- vector()
  stops <- vector()
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j+1)
      train <- rnorm(M)
      test <- rnorm(T_)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  return(sizes)
}

#Simulations

#Beta = 3/4
calculate_alpha(100, 100, 5000, 3/4,c_alphas_final_34[1])
calculate_alpha(100, 100, 5000, 3/4,c_alphas_final_34[2])
calculate_alpha(100, 100, 5000, 3/4,c_alphas_final_34[3])
calculate_alpha(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha(1000, 100, 5000, 1,c_alphas_final_1[3])

simulate_tau_sn(1000, 100, 5000, 1,c_alphas_final_1)

calculate_alpha(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha(100, 1000, 5000, 1,c_alphas_final_1[3])
simulate_tau_sn(100, 1000, 5000, 1,c_alphas_final_1)
calculate_alpha(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 1
calculate_alpha(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha(100, 100, 5000, 1,c_alphas_final_1[3])
simulate_tau_sn(100, 100, 5000, 1,c_alphas_final_1)
calculate_alpha(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha(1000, 100, 5000, 1,c_alphas_final_1[3])
simulate_tau_sn(1000, 100, 5000, 1,c_alphas_final_1)
calculate_alpha(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha(1000, 1000,5000, 2,c_alphas_final_2[1])
calculate_alpha(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha(1000, 1000, 5000, 2,c_alphas_final_2[3])

#Beta = 3
calculate_alpha(100, 100, 5000, 3,c_alphas_final_3[1])
calculate_alpha(100, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha(100, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha(1000, 100,5000,3,c_alphas_final_3[1])
calculate_alpha(1000, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha(1000, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha(100, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha(100, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha(100, 1000, 5000, 3,c_alphas_final_3[3])
calculate_alpha(1000, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha(1000, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha(1000, 1000, 5000, 3,c_alphas_final_3[3])

#Beta = 4
calculate_alpha(100, 100, 5000, 4,c_alphas_final_4[1])
calculate_alpha(100, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha(100, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha(1000, 100, 5000,4,c_alphas_final_4[1])
calculate_alpha(1000, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha(1000, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha(100, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha(100, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha(100, 1000, 5000, 4,c_alphas_final_4[3])
calculate_alpha(1000, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha(1000, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha(1000, 1000, 5000, 4,c_alphas_final_4[3])


#########iid Exponential############

calculate_alpha_exp <- function(M, t, N, beta, c_alpha){
  max_vector <- vector() 
  h = floor(t^(1/2))
  for(i in 1:N){
    set.seed(i)
    training <- rexp(M, 2)
    window_vec <- vector()
    M_mean <- mean(training)
    test <- rexp(t, 2)
    diff_vector <- vector()
    sd_estimate <- sd(training)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      #g_T <- ((sd_estimate)*(h + k)^{beta})/(h^{2*beta + 1/2})
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}
#Simulations

#Beta = 1
calculate_alpha_exp(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_exp(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_exp(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_exp(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_exp(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_exp(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_exp(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_exp(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_exp(100, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_exp(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_exp(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_exp(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_exp(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_exp(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_exp(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_exp(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_exp(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_exp(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_exp(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_exp(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_exp(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_exp(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_exp(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_exp(1000, 1000, 5000, 2,c_alphas_final_2[3])

#Beta = 3
calculate_alpha_exp(100, 100, 5000, 3,c_alphas_final_3[1])
calculate_alpha_exp(100, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_exp(100, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_exp(1000, 100, 5000,3,c_alphas_final_3[1])
calculate_alpha_exp(1000, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_exp(1000, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_exp(100, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_exp(100, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_exp(100, 1000, 5000, 3,c_alphas_final_3[3])
calculate_alpha_exp(1000, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_exp(1000, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_exp(1000, 1000, 5000, 3,c_alphas_final_3[3])

#Beta = 4
calculate_alpha_exp(100, 100, 5000, 4,c_alphas_final_4[1])
calculate_alpha_exp(100, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_exp(100, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_exp(1000, 100, 5000,4,c_alphas_final_4[1])
calculate_alpha_exp(1000, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_exp(1000, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_exp(100, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_exp(100, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_exp(100, 1000, 5000, 4,c_alphas_final_4[3])
calculate_alpha_exp(1000, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_exp(1000, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_exp(1000, 1000, 5000, 4,c_alphas_final_4[3])


###########iid T-Distribution###################
calculate_alpha_t_dist <- function(M,t, N, beta, c_alpha, df = 5){
  h = floor(t^(1/2)/beta)
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    training <- rt(M, df)
    window_vec <- vector()
    M_mean <- mean(training)
    test <- rt(t,df)
    diff_vector <- vector()
    sd_estimate <- sd(training)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- ((sd_estimate)*(h + k)^{beta})/(h^{beta + 1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}

#Beta = 1
calculate_alpha_t_dist(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_t_dist(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_t_dist(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_t_dist(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_t_dist(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_t_dist(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_t_dist(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_t_dist(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_t_dist(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_t_dist(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_t_dist(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_t_dist(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_t_dist(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_t_dist(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_t_dist(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_t_dist(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_t_dist(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_t_dist(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_t_dist(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_t_dist(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_t_dist(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_t_dist(1000, 1000,5000, 2,c_alphas_final_2[1])
calculate_alpha_t_dist(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_t_dist(1000, 1000, 5000, 2,c_alphas_final_2[3])

#Beta = 3
calculate_alpha_t_dist(100, 100, 5000, 3,c_alphas_final_3[1])
calculate_alpha_t_dist(100, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_t_dist(100, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_t_dist(1000, 100,5000,3,c_alphas_final_3[1])
calculate_alpha_t_dist(1000, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_t_dist(1000, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_t_dist(100, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_t_dist(100, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_t_dist(100, 1000, 5000, 3,c_alphas_final_3[3])
calculate_alpha_t_dist(1000, 1000,5000, 3,c_alphas_final_3[1])
calculate_alpha_t_dist(1000, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_t_dist(1000, 1000, 5000, 3,c_alphas_final_3[3])

#Beta = 4
calculate_alpha_t_dist(100, 100, 5000, 4,c_alphas_final_4[1])
calculate_alpha_t_dist(100, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_t_dist(100, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_t_dist(1000, 100,5000,4,c_alphas_final_4[1])
calculate_alpha_t_dist(1000, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_t_dist(1000, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_t_dist(100, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_t_dist(100, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_t_dist(100, 1000, 5000, 4,c_alphas_final_4[3])
calculate_alpha_t_dist(1000, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_t_dist(1000, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_t_dist(1000, 1000, 5000, 4,c_alphas_final_4[3])


###########iid Uniform###################
calculate_alpha_unif <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    training <- runif(M, max = 5)
    window_vec <- vector()
    M_mean <- mean(training)
    test <- runif(t, max = 5)
    diff_vector <- vector()
    sd_estimate <- sd(training)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      #g_T <- ((1/sqrt(12))*(h + k)^{beta}*log(t))/(sqrt(h))
      g_T <- ((sd_estimate)*(h + k)^{beta})/(h^{beta + 1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  alpha <- vector()
  for(i in 1:length(c_alpha)){
    alpha[i] <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  }
  
  return(alpha)
}

#Beta = 1
calculate_alpha_unif(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_unif(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_unif(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_unif(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_unif(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_unif(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_unif(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_unif(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_unif(100, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_unif(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_unif(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_unif(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_unif(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_unif(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_unif(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_unif(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_unif(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_unif(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_unif(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_unif(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_unif(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_unif(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_unif(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_unif(1000, 1000, 5000, 2,c_alphas_final_2[3])

#Beta = 3
calculate_alpha_unif(100, 100, 5000, 3,c_alphas_final_3[1])
calculate_alpha_unif(100, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_unif(100, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_unif(1000, 100, 5000,3,c_alphas_final_3[1])
calculate_alpha_unif(1000, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_unif(1000, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_unif(100, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_unif(100, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_unif(100, 1000, 5000, 3,c_alphas_final_3[3])
calculate_alpha_unif(1000, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_unif(1000, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_unif(1000, 1000, 5000, 3,c_alphas_final_3[3])

#Beta = 4
calculate_alpha_unif(100, 100, 5000, 4,c_alphas_final_4[1])
calculate_alpha_unif(100, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_unif(100, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_unif(1000, 100, 5000,4,c_alphas_final_4[1])
calculate_alpha_unif(1000, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_unif(1000, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_unif(100, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_unif(100, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_unif(100, 1000, 5000, 4,c_alphas_final_4[3])
calculate_alpha_unif(1000, 1000,5000, 4,c_alphas_final_4[1])
calculate_alpha_unif(1000, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_unif(1000, 1000, 5000, 4,c_alphas_final_4[3])


##############iid Gamma#######################
calculate_alpha_gamm <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    training <- rgamma(M, 1)
    window_vec <- vector()
    M_mean <- mean(training)
    test <- rgamma(t, 1)
    diff_vector <- vector()
    sd_estimate <- sd(training)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- abs(M_mean - window_vec[k])/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  alpha <- vector()
  for(i in 1:length(c_alpha)){
    alpha[i] <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  }
  
  return(alpha)
}

#Beta = 1
calculate_alpha_gamm(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_gamm(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_gamm(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_gamm(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_gamm(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_gamm(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_gamm(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_gamm(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_gamm(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_gamm(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_gamm(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_gamm(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_gamm(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_gamm(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_gamm(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_gamm(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_gamm(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_gamm(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_gamm(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_gamm(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_gamm(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_gamm(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_gamm(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_gamm(1000, 1000, 5000, 2,c_alphas_final_2[3])

#Beta = 3
calculate_alpha_gamm(100, 100, 5000, 3,c_alphas_final_3[1])
calculate_alpha_gamm(100, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_gamm(100, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_gamm(1000, 100,5000,3,c_alphas_final_3[1])
calculate_alpha_gamm(1000, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_gamm(1000, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_gamm(100, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_gamm(100, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_gamm(100, 1000, 5000, 3,c_alphas_final_3[3])
calculate_alpha_gamm(1000, 1000,5000, 3,c_alphas_final_3[1])
calculate_alpha_gamm(1000, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_gamm(1000, 1000, 5000, 3,c_alphas_final_3[3])

#Beta = 4
calculate_alpha_gamm(100, 100, 5000, 4,c_alphas_final_4[1])
calculate_alpha_gamm(100, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_gamm(100, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_gamm(1000, 100,5000,4,c_alphas_final_4[1])
calculate_alpha_gamm(1000, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_gamm(1000, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_gamm(100, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_gamm(100, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_gamm(100, 1000, 5000, 4,c_alphas_final_4[3])
calculate_alpha_gamm(1000, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_gamm(1000, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_gamm(1000, 1000, 5000, 4,c_alphas_final_4[3])

##############iid F-Distribution#######################
calculate_alpha_F <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2)/beta)
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i+1)
    training <- rf(M, 6, 7)
    window_vec <- vector()
    M_mean <- mean(training)
    test <- rf(t, 6, 7)
    diff_vector <- vector()
    sd_estimate <- sd(training)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j] 
      }
      window_vec[k] <- sum/h
      
          g_T <- (sd_estimate*(h + k)^{beta})/(h^{beta + 1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  alpha <- vector()
  for(i in 1:length(c_alpha)){
    alpha[i] <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  }

  
  return(alpha)
}

#Beta = 3/4
calculate_alpha_F(100, 100, 5000, 3/4,c_alphas_final_34[1])
calculate_alpha_F(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_F(1000, 100,5000, 1,c_alphas_final_1[1])
calculate_alpha_F(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_F(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_F(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 1
calculate_alpha_F(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_F(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_F(1000, 100,5000, 1,c_alphas_final_1[1])
calculate_alpha_F(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_F(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_F(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_F(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_F(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_F(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_F(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_F(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_F(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_F(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_F(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_F(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_F(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_F(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_F(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_F(1000, 1000, 5000, 2,c_alphas_final_2[3])

#Beta = 3
calculate_alpha_F(100, 100, 5000, 3,c_alphas_final_3[1])
calculate_alpha_F(100, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_F(100, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_F(1000, 100, 5000,3,c_alphas_final_3[1])
calculate_alpha_F(1000, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_F(1000, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_F(100, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_F(100, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_F(100, 1000, 5000, 3,c_alphas_final_3[3])
calculate_alpha_F(1000, 1000,5000, 3,c_alphas_final_3[1])
calculate_alpha_F(1000, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_F(1000, 1000, 5000, 3,c_alphas_final_3[3])

#Beta = 4
calculate_alpha_F(100, 100, 5000, 4,c_alphas_final_4[1])
calculate_alpha_F(100, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_F(100, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_F(1000, 100,5000,4,c_alphas_final_4[1])
calculate_alpha_F(1000, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_F(1000, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_F(100, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_F(100, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_F(100, 1000, 5000, 4,c_alphas_final_4[3])
calculate_alpha_F(1000, 1000,5000, 4,c_alphas_final_4[1])
calculate_alpha_F(1000, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_F(1000, 1000, 5000, 4,c_alphas_final_4[3])


##############iid Chi-Square Distribution#######################
calculate_alpha_Chi <- function(M, t, N, beta, c_alpha, df = 5){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    training <- rchisq(M, df)
    window_vec <- vector()
    M_mean <- mean(training)
    test <- rchisq(t, df)
    diff_vector <- vector()
    sd_estimate <- sd(training)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1
calculate_alpha_Chi(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_Chi(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_Chi(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_Chi(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_Chi(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_Chi(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_Chi(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_Chi(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_Chi(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_Chi(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_Chi(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_Chi(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_Chi(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_Chi(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_Chi(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_Chi(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_Chi(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_Chi(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_Chi(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_Chi(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_Chi(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_Chi(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_Chi(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_Chi(1000, 1000, 5000, 2,c_alphas_final_2[3])

#Beta = 3
calculate_alpha_Chi(100, 100, 5000, 3,c_alphas_final_3[1])
calculate_alpha_Chi(100, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_Chi(100, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_Chi(1000, 100, 5000,3,c_alphas_final_3[1])
calculate_alpha_Chi(1000, 100, 5000, 3,c_alphas_final_3[2])
calculate_alpha_Chi(1000, 100, 5000, 3,c_alphas_final_3[3])
calculate_alpha_Chi(100, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_Chi(100, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_Chi(100, 1000, 5000, 3,c_alphas_final_3[3])
calculate_alpha_Chi(1000, 1000, 5000, 3,c_alphas_final_3[1])
calculate_alpha_Chi(1000, 1000, 5000, 3,c_alphas_final_3[2])
calculate_alpha_Chi(1000, 1000, 5000, 3,c_alphas_final_3[3])

#Beta = 4
calculate_alpha_Chi(100, 100, 5000, 4,c_alphas_final_4[1])
calculate_alpha_Chi(100, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_Chi(100, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_Chi(1000, 100,5000,4,c_alphas_final_4[1])
calculate_alpha_Chi(1000, 100, 5000, 4,c_alphas_final_4[2])
calculate_alpha_Chi(1000, 100, 5000, 4,c_alphas_final_4[3])
calculate_alpha_Chi(100, 1000, 5000, 4,c_alphas_final_4[1])
calculate_alpha_Chi(100, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_Chi(100, 1000, 5000, 4,c_alphas_final_4[3])
calculate_alpha_Chi(1000, 1000,5000, 4,c_alphas_final_4[1])
calculate_alpha_Chi(1000, 1000, 5000, 4,c_alphas_final_4[2])
calculate_alpha_Chi(1000, 1000, 5000, 4,c_alphas_final_4[3])



##################Function of M#########################
#Beta = 2, T = 100, h = 2
M_Vec <- c(50:100)
M_sim <- vector()
for(i in 1:length(M_Vec-1)){
M_sim[i] <- calculate_alpha(M_Vec[i], 2, 800, 1000, 2)
}

plot(M_Vec, M_sim, main = bquote(alpha ~ "values as function of " ~ M), xlab = bquote(M), ylab = bquote(alpha))


##################Function of T#########################
#Beta = 2, M = 100, h = 2
T_Vec <- c(50:100)
T_sim <- vector()
for(i in 1:length(T_Vec-1)){
  T_sim[i] <- calculate_alpha(100, 2, T_Vec[i], 1000, 2)
}

plot(T_Vec, T_sim, main = bquote(alpha ~ "values as function of " ~ T), xlab = bquote(T), ylab = bquote(alpha))


##################Function of h#########################
#Beta = 2, T = 100, M = 100
h_Vec <- c(2:50)
h_sim <- vector()
for(i in 1:length(h_Vec-1)){
  h_sim[i] <- calculate_alpha(100, h_Vec[i], 100 , 1000, 2)
}

plot(h_Vec, h_sim, main = bquote(alpha ~ "values as function of " ~ h), xlab = bquote(h), ylab = bquote(alpha))



#Simulate Test statistics for different types of models


#############################################
#
#############################################
sample_test_statistic <- function(train, test, window){
  x_bar <- mean(train)
  t <- length(test)
  window_vec <- vector()
  diff_vector <- vector()
  
  sd_estimate <- sd(train)
  
  for(k in 1:(t-window)){
    sum =0 
    for(j in k:(k+window)){
      sum = sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^2)/(window^{2 + 2})
    
    diff_vector[k] <- abs(x_bar - window_vec[k])/(g_T)
  }
  
  maximum <- max(diff_vector)
  
  return(maximum)
}



###########2 Normals with one mean 0 and one with different mean################

normal_zero_sample <- rnorm(100)

normal_one_sample <- rnorm(100, 3)

sd(normal_zero_sample)

sample_test_statistic(normal_zero_sample, normal_one_sample, 2)


###########1 Normals with one mean 0 and one exponential################

normal_zero_sample <- rnorm(100, 0)

exp_2 <- rexp(100, 3)

exp_1 <- rexp(100, 1)

sample_test_statistic(exp_1, exp_2, 2)



###########Difference in Uniforms################


unif_1 <- runif(100)

unif_2 <- runif(100, 0, 2)

sample_test_statistic(unif_1, unif_2, 2)




###########1 Normals with one mean 0 and one Chi-Square################

normal_zero_sample <- rnorm(100, 0)

chi_1 <- rchisq(100, 1)

sample_test_statistic(normal_zero_sample, chi_1, 2)

plot(exp_1)

############AR(1), AR(2)#################

ar_1 <- arima.sim(model=list(ar = c(.9)), n = 100)


ar_2 <-  arima.sim(model= list(ar = c(.9)), n = 100)

plot(ar_2)

sample_test_statistic(ar_1, ar_2, 2)

############MA(1), MAd(2)#################

ma_1 <- 1+ arima.sim(model=list(ma = c(.5)), n = 100)


ma_2 <-  arima.sim(model= list(ma = c(.5)), n = 200)

sample_test_statistic(ma_1, ma_2, 2)
#############GARCH##############
ar_1 <- garch.sim(model=list(ar = c(.9, -.2)), n = 100)


ar_2 <- garch.sim(model= list(ar = c(.5)), n = 100)


###########ARMA#########################

arma_1 <- 3 + arima.sim(model=list(ar = c(.9, -.2),ma=c(-.7,.1)), n = 100)


arma_2 <- arima.sim(model= list(ar = c(.5),ma=c(-.7,.1)), n = 100)

sample_test_statistic(arma_1, arma_2, 2)

#############Changes in Variance###############

###########2 Normals with one mean 0 and one with different variance################

normal_zero_sample <- rnorm(100)

normal_one_sample <- rnorm(100, 0, 5)


sample_test_statistic(normal_zero_sample, normal_one_sample, 2)

###########################Optimal Stopping Time#######################

#############################################
#
#############################################
sample_test_statistic_stop_time <- function(train, test, beta, c_alpha){
  
  x_bar <- mean(train)
  t <- length(test)
  window <- floor(t^(1/2))
  #window <- 2
  window_vec <- vector()
  diff_vector <- vector()
  
  stop_obs <- 0
  
  sd_estimate <- sd(train)
  
  for(k in 1:(t-window)){
    sum =0 
    
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    #g_T <- sd_estimate/(sqrt(window*t))
    if(abs(x_bar - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- k
      break
    }
  }
  
  
  return(stop_obs)
}

#############################################
#Perform sequential monitoring procedure
#M
#beta
#c_alphas
#N
#############################################
simulate_tau <- function(M, beta, c_alphas, N){
  sizes <- vector()
  stops <- vector()
  T_ <- 100
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j)
      train <- rnorm(M)
      test <- rnorm(T_)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  kbar <- mean(stops)
  return(c(sizes, kbar))
}
#Under the null
simulate_tau(100, 1, c_alphas_final_1, 5000)
simulate_tau(100, 2, c_alphas_final_2, 5000)
simulate_tau(100, 3, c_alphas_final_3, 5000)
simulate_tau(100, 4, c_alphas_final_4, 5000)


#Under the alternative

simulate_tau_alt <- function(M, mean_test, beta, c_alphas, N){
  sizes <- vector()
  stops <- vector()
  T_ <- 100
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j)
      train <- rnorm(M)
      test <- rnorm(T_, mean_test)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  return(sizes)
}

#Beta = 1
simulate_tau_alt(100, 1, 1, c_alphas_final_1, 5000)
simulate_tau_alt(100, 1, 2, c_alphas_final_2, 5000)
simulate_tau_alt(100, 1, 3, c_alphas_final_3, 5000)
simulate_tau_alt(100, 1, 4, c_alphas_final_4, 5000)

#Beta = 2
simulate_tau_alt(100, 2, 1, c_alphas_final_1, 5000)
simulate_tau_alt(100, 2, 2, c_alphas_final_2, 5000)
simulate_tau_alt(100, 2, 3, c_alphas_final_3, 5000)
simulate_tau_alt(100, 2, 4, c_alphas_final_4, 5000)

#Beta = 3
simulate_tau_alt(100, 3, 1, c_alphas_final_1, 5000)
simulate_tau_alt(100, 3, 2, c_alphas_final_2, 5000)
simulate_tau_alt(100, 3, 3, c_alphas_final_3, 5000)
simulate_tau_alt(100, 3, 4, c_alphas_final_4, 5000)

#Beta = 4
simulate_tau_alt(100, 4, 1, c_alphas_final_1, 5000)
simulate_tau_alt(100, 4, 2, c_alphas_final_2, 5000)
simulate_tau_alt(100, 4, 3, c_alphas_final_3, 5000)
simulate_tau_alt(100, 4, 4, c_alphas_final_4, 5000)

#############################################
#
#############################################
simulate_tau_chi <- function(M, beta, c_alphas, N){
  sizes <- vector()
  stops <- vector()
  T_ <- 100
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j)
      train <- rnorm(M)
      test <- rchisq(T_, 2)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  return(sizes)
}


#Beta = 1
simulate_tau_chi(100, 1, c_alphas_final_1, 5000)
simulate_tau_chi(100, 2, c_alphas_final_2, 5000)
simulate_tau_chi(100, 3, c_alphas_final_3, 5000)
simulate_tau_chi(100, 4, c_alphas_final_4, 5000)



#############################################
#
#############################################
simulate_tau_gamm <- function(M, beta, c_alphas, N){
  sizes <- vector()
  stops <- vector()
  T_ <- 100
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j)
      train <- rnorm(M)
      test <- rgamma(T_, shape = 5,rate = 5)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  return(c(sizes))
}



simulate_tau_gamm(100, 1, c_alphas_final_1, 5000)
simulate_tau_gamm(100, 2, c_alphas_final_2, 5000)
simulate_tau_gamm(100, 3, c_alphas_final_3, 5000)
simulate_tau_gamm(100, 4, c_alphas_final_4, 5000)


#############################################
#
#############################################
simulate_tau_exp <- function(M, beta, c_alphas, N){
  sizes <- vector()
  stops <- vector()
  T_ <- 100
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j)
      train <- rnorm(M)
      test <- rexp(T_, rate = 1/5)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  return(sizes)
}



simulate_tau_exp(100, 1, c_alphas_final_1, 5000)
simulate_tau_exp(100, 2, c_alphas_final_2, 5000)
simulate_tau_exp(100, 3, c_alphas_final_3, 5000)
simulate_tau_exp(100, 4, c_alphas_final_4, 5000)


#############################################
#
#############################################
simulate_tau_unif <- function(M, beta, c_alphas, N){
  sizes <- vector()
  stops <- vector()
  T_ <- 100
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j)
      train <- rnorm(M)
      test <- runif(T_, max = 5)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  kbar <- mean(stops)
  return(c(sizes, kbar))
}


#############################################
#
#############################################
simulate_tau_unif(100, 1, c_alphas_final_1, 5000)
simulate_tau_unif(100, 2, c_alphas_final_2, 5000)
simulate_tau_unif(100, 3, c_alphas_final_3, 5000)
simulate_tau_unif(100, 4, c_alphas_final_4, 5000)

simulate_tau_pois <- function(M, beta, c_alphas, N){
  sizes <- vector()
  stops <- vector()
  T_ <- 100
  for(i in 1:length(c_alphas)){
    for(j in 1:N){
      set.seed(j)
      train <- rnorm(M)
      test <- rpois(T_, 1)
      stops[j] <- sample_test_statistic_stop_time(train, test, beta, c_alphas[i])
    }
    sizes[i] <- length(stops[stops != 0])/length(stops)
  }
  kbar <- mean(stops)
  return(c(sizes, kbar))
}


#############################################
#
#############################################
simulate_tau_pois(100, 1, c_alphas_final_1, 5000)
simulate_tau_pois(100, 2, c_alphas_final_2, 5000)
simulate_tau_pois(100, 3, c_alphas_final_3, 5000)
simulate_tau_pois(100, 4, c_alphas_final_4, 5000)



sample_test_statistic_stop_time(normal_zero_sample, normal_one_sample, 2, 1.68)

####################################Comparing with pdf of Maximums over Absolute Value of Normals###########################  

maximums <- vector()
for(i in 1:1000){
  sample <- rnorm(1000)
  sd_estimate <- sd(sample)
  gauss <- abs(sample)
  maximums[i] <- max(gauss/(sd_estimate*sqrt(1000)))
}



plot(density(maximums), main = bquote("Density Estimate"), xlab = "x")


cv3 <- unname(quantile(maximums, c(.90, .95, .99))) 




tau_probability <- function (N, M, beta, c){
  max_vector <- W_process_mod(N, M, beta)/c
  
  prob <- length(max_vector[max_vector <= 1])/length(max_vector)
  
  return(prob)
}

tau_probability(1000, 100, 2, 3)


#Dependent Data


#AR(1)

set.seed(101)

#############################################
#
#############################################
calculate_alpha_AR <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  phi <- 0.2
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar1_train = 1+arima.sim(model=list(ar=phi),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = 1+arima.sim(model=list(ar=phi),n = t, innov=err_test)
    diff_vector <- vector()
    #sd_estimate <- 1/(1-0.75)
    sd_estimate <- 1/(1-phi)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1
calculate_alpha_AR(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_AR(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_AR(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR(100, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_AR(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_AR(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_AR(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_AR(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_AR(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(1000, 1000, 5000, 2,c_alphas_final_2[3])



#############################################
#
#############################################
#MA(1)
calculate_alpha_MA <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  theta <- 0.3
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar1_train = arima.sim(model=list(ma=theta),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ma=theta),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- (1+theta)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}

#Beta = 3/4
calculate_alpha_MA(100, 100, 5000, 3/4,c_alphas_final_34[1])
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 1
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_MA(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_MA(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_MA(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_MA(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_MA(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_MA(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_MA(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_MA(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_MA(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_MA(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_MA(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_MA(1000, 1000, 5000, 2,c_alphas_final_2[3])

#############################################
#
#############################################
#AR(2)
calculate_alpha_AR_2 <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar2_train = arima.sim(model=list(ar=c(0.75, 0.55)),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ar=c(0.75, 0.55)),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- 1/(1-0.75-0.55)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1
calculate_alpha_AR_2(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR_2(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR_2(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_AR_2(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR_2(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR_2(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_AR_2(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR_2(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR_2(1000, 1000,5000, 1,c_alphas_final_1[3])
calculate_alpha_AR_2(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR_2(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR_2(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_AR(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_AR(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_AR(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_AR(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR(1000, 1000, 5000, 2,c_alphas_final_2[3])

#############################################
#
#############################################
calculate_alpha_MA_2 <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar2_train = arima.sim(model=list(ma=c(0.75, 0.55)),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ma=c(0.75, 0.55)),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- (1+0.75+0.55)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}

#############################################
#
#############################################
#Beta = 1
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_MA(100, 100, 2000, 2,c_alphas_final_2[1])
calculate_alpha_MA(100, 100, 2000, 2,c_alphas_final_2[2])
calculate_alpha_MA(100, 100, 2000, 2,c_alphas_final_2[3])
calculate_alpha_MA(1000, 100, 2000, 2,c_alphas_final_2[1])
calculate_alpha_MA(1000, 100, 2000, 2,c_alphas_final_2[2])
calculate_alpha_MA(1000, 100, 2000, 2,c_alphas_final_2[3])
calculate_alpha_MA(100, 1000, 2000, 2,c_alphas_final_2[1])
calculate_alpha_MA(100, 1000, 2000, 2,c_alphas_final_2[2])
calculate_alpha_MA(100, 1000, 2000, 2,c_alphas_final_2[3])
calculate_alpha_MA(1000, 1000, 2000, 2,c_alphas_final_2[1])
calculate_alpha_MA(1000, 1000, 2000, 2,c_alphas_final_2[2])
calculate_alpha_MA(1000, 1000, 2000, 2,c_alphas_final_2[3])

#ARMA(1,1)
#############################################
#
#############################################
calculate_alpha_ARMA <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_arma_train = arima.sim(model=list(ma=0.3, ar=0.2),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ma=0.3, ar=0.2),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- (1+0.3)/(1 - 0.2)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1
calculate_alpha_ARMA(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARMA(1000, 100,5000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARMA(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARMA(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA(1000, 1000, 5000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_ARMA(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA(1000, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA(1000, 1000, 5000, 2,c_alphas_final_2[3])

#ARMA(2,2)
#############################################
#
#############################################

calculate_alpha_ARMA_2 <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_arma_train = arima.sim(model=list(ma=c(0.75, 0.6), ar=c(0.75, 0.45)),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ma=c(0.75, 0.6), ar=c(0.75, 0.45)),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- (1+0.75+0.6)/(1 - 0.75 - 0.45)
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1
calculate_alpha_ARMA_2(100, 100, 2000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA_2(100, 100, 2000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA_2(100, 100, 2000, 1,c_alphas_final_1[3])
calculate_alpha_ARMA_2(1000, 100, 2000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA_2(1000, 100, 2000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA_2(1000, 100, 2000, 1,c_alphas_final_1[3])
calculate_alpha_ARMA_2(100, 1000, 2000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA_2(100, 1000, 2000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA_2(1000, 1000, 2000, 1,c_alphas_final_1[3])
calculate_alpha_ARMA_2(1000, 1000, 2000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA_2(1000, 1000, 2000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA_2(1000, 1000, 2000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_ARMA_2(100, 100, 2000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA_2(100, 100, 2000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA_2(100, 100, 2000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA_2(1000, 100, 2000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA_2(1000, 100, 2000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA_2(1000, 100, 2000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA_2(100, 1000, 2000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA_2(100, 1000, 2000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA_2(100, 1000, 2000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA_2(1000, 1000, 2000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA_2(1000, 1000, 2000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA_2(1000, 1000, 2000, 2,c_alphas_final_2[3])


#############################################
#
#############################################
calculate_alpha_ARCH <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    alpha_0 <- 0.25
    omega <- 0.3
    spec <- garchSpec(model = list(omega = omega, alpha = c(alpha_0), beta = 0))
    y_arch_train <- garchSim(spec, n = M)
    diff_vector <- vector()
    window_vec <- vector()
    M_mean <- mean(y_arch_train)
    test <- garchSim(spec, n = t)
    sd_estimate <- sqrt((omega)/(1 - alpha_0))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}

#Beta 1

calculate_alpha_ARCH(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARCH(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARCH(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARCH(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARCH(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARCH(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARCH(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARCH(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARCH(100, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARCH(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARCH(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARCH(1000, 1000, 5000, 1,c_alphas_final_1[3])


#Beta 2

calculate_alpha_ARCH(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARCH(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARCH(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARCH(1000, 100, 5000,2,c_alphas_final_2[1])
calculate_alpha_ARCH(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARCH(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARCH(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARCH(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARCH(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARCH(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARCH(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARCH(1000, 1000, 5000, 2,c_alphas_final_2[3])


#Do same simulations for ARMA and ARIMA


calculate_alpha_GARCH <- function(M, t, N, beta, c_alpha){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    alpha_0 <- 0.25
    omega <- 0.3
    beta_0 <- 0.1
    spec <- garchSpec(model = list(omega = omega, alpha = c(alpha_0), beta = beta_0))
    y_garch_train <- garchSim(spec, n = M)
    diff_vector <- vector()
    window_vec <- vector()
    M_mean <- mean(y_garch_train)
    test <- garchSim(spec, n = t)
    sd_estimate <- sqrt((omega)/(1 - alpha_0 - beta_0))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}

#Beta 1

calculate_alpha_GARCH(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_GARCH(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_GARCH(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_GARCH(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_GARCH(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_GARCH(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_GARCH(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_GARCH(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_GARCH(100, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_GARCH(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_GARCH(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_GARCH(1000, 1000, 5000, 1,c_alphas_final_1[3])


#Beta 2

calculate_alpha_GARCH(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_GARCH(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_GARCH(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_GARCH(1000, 100, 5000,2,c_alphas_final_2[1])
calculate_alpha_GARCH(1000, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_GARCH(1000, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_GARCH(100, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_GARCH(100, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_GARCH(100, 1000, 5000, 2,c_alphas_final_2[3])
calculate_alpha_GARCH(1000, 1000, 5000, 2,c_alphas_final_2[1])
calculate_alpha_GARCH(1000, 1000, 5000, 2,c_alphas_final_2[2])
calculate_alpha_GARCH(1000, 1000, 5000, 2,c_alphas_final_2[3])




########################################

#############################################
#
#############################################

calculate_alpha_AR_lrv <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  phi <- 0.2
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar1_train = arima.sim(model=list(ar=phi),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ar=phi),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}



#Beta = 1 
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Flat-Top")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Flat-Top")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Flat-Top")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Parzen")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Parzen")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Parzen")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Quadratic Spectral")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Quadratic Spectral")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Quadratic Spectral")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Tukey-Hanning")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Tukey-Hanning")
calculate_alpha_AR_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Tukey-Hanning")



#Beta = 2
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Flat-Top")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Flat-Top")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Flat-Top")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Parzen")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Parzen")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Parzen")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Quadratic Spectral")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Quadratic Spectral")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Quadratic Spectral")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Tukey-Hanning")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Tukey-Hanning")
calculate_alpha_AR_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Tukey-Hanning")

#############################################
#
#############################################
calculate_alpha_AR_lrv_2 <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar1_train = arima.sim(model=list(ar=c(0.75, -0.35)),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ar=c(0.75, -0.35)),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Parzen")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Parzen")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Parzen")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Flat-Top")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Flat-Top")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Flat-Top")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Quadratic Spectral")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Quadratic Spectral")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Quadratic Spectral")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Truncated")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Truncated")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Truncated")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Tukey-Hanning")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Tukey-Hanning")
calculate_alpha_AR_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Tukey-Hanning")
calculate_alpha_AR_lrv_2(1000, 100, 1000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv_2(1000, 100, 1000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv_2(1000, 100, 1000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv_2(100, 1000, 1000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv_2(100, 1000, 1000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv_2(1000, 1000, 1000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv_2(1000, 1000, 2000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv_2(1000, 1000, 2000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv_2(1000, 1000, 2000, 1,c_alphas_final_1[3])


#############################################
#
#############################################
#MA(1)
calculate_alpha_MA_lrv <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar1_train = arima.sim(model=list(ma=0.3),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ma=0.3),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1 
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Parzen")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Parzen")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Parzen")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Flat-Top")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Flat-Top")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Flat-Top")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Quadratic Spectral")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Quadratic Spectral")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Quadratic Spectral")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Tukey-Hanning")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Tukey-Hanning")
calculate_alpha_MA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Tukey-Hanning")

#Beta = 2 
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Parzen")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Parzen")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Parzen")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Flat-Top")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Flat-Top")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Flat-Top")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Quadratic Spectral")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Quadratic Spectral")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Quadratic Spectral")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Truncated")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Truncated")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Truncated")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Tukey-Hanning")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Tukey-Hanning")
calculate_alpha_MA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Tukey-Hanning")
calculate_alpha_MA_lrv(1000, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA_lrv(1000, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA_lrv(1000, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA_lrv(100, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA_lrv(100, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA_lrv(1000, 1000, 5000, 1,c_alphas_final_1[3])
calculate_alpha_MA_lrv(1000, 1000, 5000, 1,c_alphas_final_1[1])
calculate_alpha_MA_lrv(1000, 1000, 5000, 1,c_alphas_final_1[2])
calculate_alpha_MA_lrv(1000, 1000, 5000, 1,c_alphas_final_1[3])


#############################################
#
#############################################
######MA(2)
calculate_alpha_MA_lrv_2 <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar1_train = arima.sim(model=list(ma=c(0.75, 0.55)),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ma=c(0.75, 0.55)),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1])
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2])
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3])
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Parzen")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Parzen")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Parzen")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Flat-Top")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Flat-Top")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Flat-Top")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Quadratic Spectral")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Quadratic Spectral")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Quadratic Spectral")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Truncated")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Truncated")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Truncated")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[1], "Tukey-Hanning")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[2], "Tukey-Hanning")
calculate_alpha_MA_lrv_2(100, 100, 1000, 1,c_alphas_final_1[3], "Tukey-Hanning")
calculate_alpha_AR_lrv(1000, 100, 1000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv(1000, 100, 1000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv(1000, 100, 1000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv(100, 1000, 1000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv(100, 1000, 1000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv(1000, 1000, 1000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv(1000, 1000, 2000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv(1000, 1000, 2000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv(1000, 1000, 2000, 1,c_alphas_final_1[3])

#Beta = 2
calculate_alpha_MA_lrv_2(1000, 100, 2000, 2,c_alphas_final_2[1])
calculate_alpha_MA_lrv_2(1000, 100, 2000, 2,c_alphas_final_2[2])
calculate_alpha_MA_lrv_2(1000, 100, 2000, 2,c_alphas_final_2[3])


######ARMA(1,1)
calculate_alpha_ARMA_lrv <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err <- rnorm(M)
    y_ar1_train <- arima.sim(model=list(ar=0.2, ma=0.3),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ar=0.2, ma=0.3),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


#Beta = 1 
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Tukey-Hanning")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Tukey-Hanning")
calculate_alpha_ARMA_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Tukey-Hanning")



#Beta = 2
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Tukey-Hanning")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Tukey-Hanning")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Tukey-Hanning")







#############################################
#
#############################################
calculate_alpha_ARMA_lrv2 <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err = rnorm(M)
    y_ar1_train = arima.sim(model=list(ar=c(0.75, -0.35), ma=c(0.75, 0.55)),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ar=c(0.75, -0.35), ma=c(0.75, 0.55)),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}



calculate_alpha_ARMA_lrv2(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARMA_lrv2(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARMA_lrv2(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Parzen")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Flat-Top")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Quadratic Spectral")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Tukey-Hanning")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Tukey-Hanning")
calculate_alpha_ARMA_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Tukey-Hanning")
calculate_alpha_AR_lrv(1000, 100, 1000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv(1000, 100, 1000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv(100, 1000, 1000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv(100, 1000, 1000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv(1000, 1000, 1000, 1,c_alphas_final_1[3])
calculate_alpha_AR_lrv(1000, 1000, 2000, 1,c_alphas_final_1[1])
calculate_alpha_AR_lrv(1000, 1000, 2000, 1,c_alphas_final_1[2])
calculate_alpha_AR_lrv(1000, 1000, 2000, 1,c_alphas_final_1[3])






#############################################
#
#############################################
######ARMA(1,1)
calculate_alpha_ARCH_lrv <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    alpha_0 <- 0.25
    omega <- 0.3
    spec <- garchSpec(model = list(omega = omega, alpha = c(alpha_0), beta = 0))
    y_arch_train <- garchSim(spec, n = M)
    diff_vector <- vector()
    window_vec <- vector()
    M_mean <- mean(y_arch_train)
    test <- garchSim(spec, n = t)
    sd_estimate <- sqrt(asymp_var(y_arch_train^2, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}

#Beta = 1
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Parzen")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Parzen")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Parzen")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Flat-Top")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Flat-Top")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Flat-Top")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Quadratic Spectral")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Quadratic Spectral")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Quadratic Spectral")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Tukey-Hanning")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Tukey-Hanning")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Tukey-Hanning")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "None")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "None")
calculate_alpha_ARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "None")




#Beta = 2
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Parzen")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Parzen")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Parzen")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Flat-Top")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Flat-Top")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Flat-Top")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Quadratic Spectral")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Quadratic Spectral")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Quadratic Spectral")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Tukey-Hanning")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Tukey-Hanning")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Tukey-Hanning")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "None")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "None")
calculate_alpha_ARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "None")




######ARMA(1,1)
calculate_alpha_GARCH_lrv <- function(M, t, N, beta, c_alpha, kernel = "Bartlett"){
  h = floor(M^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    alpha_0 <- 0.25
    omega <- 0.3
    beta_0 <- 0.1
    spec <- garchSpec(model = list(omega = omega, alpha = c(alpha_0), beta = beta_0))
    y_arch_train <- garchSim(spec, n = M)
    diff_vector <- vector()
    window_vec <- vector()
    M_mean <- mean(y_arch_train)
    test <- garchSim(spec, n = t)
    sd_estimate <- sqrt(asymp_var(y_arch_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}



#Beta = 1
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1])
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2])
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3])
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Parzen")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Parzen")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Parzen")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Flat-Top")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Flat-Top")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Flat-Top")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Quadratic Spectral")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Quadratic Spectral")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Quadratic Spectral")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "Tukey-Hanning")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "Tukey-Hanning")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "Tukey-Hanning")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[1], "None")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[2], "None")
calculate_alpha_GARCH_lrv(100, 100, 5000, 1,c_alphas_final_1[3], "None")


#Beta = 2
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1])
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2])
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3])
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Parzen")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Parzen")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Parzen")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Flat-Top")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Flat-Top")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Flat-Top")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Quadratic Spectral")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Quadratic Spectral")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Quadratic Spectral")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "Tukey-Hanning")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "Tukey-Hanning")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "Tukey-Hanning")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[1], "None")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[2], "None")
calculate_alpha_GARCH_lrv(100, 100, 5000, 2,c_alphas_final_2[3], "None")





################################
#Long run variance estimation methods
################################


asymp_var <- function(x, ker = c("Bartlett", "Parzen", "Tukey-Hanning", "Quadratic Spectral", "Flat-Top")){
  T_ <- length(x)
  M_T <- floor(4*(T_/100)^{2/9})
  auto <- acf(x, lag.max = M_T + 1,type = c("covariance"), plot = FALSE)
  covar <- c(auto$acf)
  gamma_0 <- covar[1]
  gamma_sum <- vector()
  gamma_sum[1] <- 0
  for(i in 2:(M_T + 1)){
    index <- i/(M_T + 1)
    if(ker == "Flat-Top"){gamma_sum[i] <- flat_top_ft(index)*covar[i] } else if(ker == "None") { break }else{
    gamma_sum[i] <- kweights(index, kernel = ker, normalize = FALSE)*covar[i]
    }
  }
  lrv <- gamma_0 + 2*sum(gamma_sum)
  
  return(lrv)
  
}

flat_top_ft <- function(x){
  if(abs(x)<= 1/2){
    return(1)
  }
  else{
    2*(1-abs(x))
  }
}



##########################################################
#Time Series Stopping Times
##########################################################

###############################################
#ARCH/GARCH
###############################################



plot_arch_sample_statistic_stop_time <- function(M, t,beta, c_alpha){
  alpha_0 <- 0.25
  omega <- 0.3
  train_spec <- garchSpec(model = list(omega = omega, alpha = c(alpha_0), beta = 0))
  y_arch_train <- garchSim(train_spec, n = M)
  M_mean <- mean(y_arch_train)
  test_spec <- garchSpec(model = list(omega = 1, alpha = c(alpha_0), beta = 0))
  test <- garchSim(test_spec, n = t)
  window <- floor(t^(1/2))
  window_vec <- vector()
  diff_vector <- vector()
  
  stop_obs <- 0
  
  sd_estimate <- sqrt(asymp_var(y_arch_train, ker = "None"))
  
  for(k in 1:(t-window)){
    sum =0 
    
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    if(abs(M_mean - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- k
      break
    }
  }
  
  k_val <- M + k
  title <- bquote("ARCH(1) Stopping Time,  " ~ omega == .(omega) ~ ", " ~ alpha == .(alpha_0) ~", Test " ~ omega == .(1))
  plot.ts(c(y_arch_train, test), main = title)
  abline(v = k_val, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/change_means/", "arch_stoptime1.pdf", sep = "")
  pdf(file_name)
  plot.ts(c(y_arch_train, test), ylab = "ARCH(1)", main = title, xlab = "Index")
  abline(v = k_val, col = "red")
  dev.off()
  
  
}


#############################################
#
#############################################
plot_garch_sample_statistic_stop_time <- function(M, t,beta, c_alpha){
  alpha_0 <- 0.25
  omega <- 0.3
  beta_0 = 0.1
  omega_2 = 1
  train_spec <- garchSpec(model = list(omega = omega, alpha = c(alpha_0), beta = beta_0))
  y_arch_train <- garchSim(train_spec, n = M)
  M_mean <- mean(y_arch_train)
  test_spec <- garchSpec(model = list(omega = 1, alpha = c(alpha_0), beta = beta_0))
  test <- garchSim(test_spec, n = t)
  window <- floor(t^(1/2))
  window_vec <- vector()
  diff_vector <- vector()
  
  stop_obs <- 0
  
  sd_estimate <- sqrt(asymp_var(y_arch_train, ker = "None"))
  
  for(k in 1:(t-window)){
    sum =0 
    
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    if(abs(M_mean - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- k
      break
    }
  }
  
  k_val <- M + k
  title <- bquote("GARCH(1,1),  " ~ omega == .(omega) ~ ", " ~ alpha == .(alpha_0) ~", Test " ~ omega == .(omega_2))
  plot.ts(c(y_arch_train, test), main = title)
  abline(v = k_val, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/change_means/", "garch_stoptime1.pdf", sep = "")
  pdf(file_name)
  plot.ts(c(y_arch_train, test), ylab = "GARCH(1,1)", main = title)
  abline(v = k_val, col = "red")
  dev.off()
  
  
}

#############################################
#
#############################################
plot_AR_sample_statistic_stop_time <- function(M, t,beta, c_alpha){
  phi <- 0.2
  window <- floor(t^(1/2))
  err = rnorm(M)
  y_ar1_train = arima.sim(model=list(ar=phi),n = M,innov=err)
  window_vec <- vector()
  M_mean <- mean(y_ar1_train)
  phi_2 <- 0.3
    err_test = rnorm(t)
    test = arima.sim(model=list(ar=phi_2),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = "Parzen"))
  
  for(k in 1:(t-window)){
    sum =0 
    
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    if(abs(M_mean - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- k
      break
    }
  }
  
  k_val <- M + k
  title <- bquote("AR(1) Stopping Time,  Train " ~ phi == .(phi) ~ ", Test " ~ phi == .(phi_2))
  plot.ts(c(y_arch_train, test), main = title)
  abline(v = k_val, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/change_means/", "ar1_stoptime1.pdf", sep = "")
  pdf(file_name)
  plot.ts(c(y_arch_train, test), ylab = "AR(1)", xlab = "Index", main = title)
  abline(v = k_val, col = "red")
  dev.off()
  
  
}
#############################################
#
#############################################
plot_MA_sample_statistic_stop_time <- function(M, t,beta, c_alpha){
  theta <- 0.2
  window <- floor(t^(1/2))
  err = rnorm(M)
  y_ma1_train = arima.sim(model=list(ma=theta),n = M,innov=err)
  window_vec <- vector()
  M_mean <- mean(y_ma1_train)
  theta_2 <- 0.3
  err_test = rnorm(t)
  test = arima.sim(model=list(ma=theta_2),n = t, innov=err_test)
  diff_vector <- vector()
  sd_estimate <- sqrt(asymp_var(y_ma1_train, ker = "Parzen"))
  
  for(k in 1:(t-window)){
    sum =0 
    
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    if(abs(M_mean - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- k
      break
    }
  }
  
  k_val <- M + k
  title <- bquote("MA(1) Stopping Time,  Train " ~ theta == .(theta) ~ ", Test " ~ theta == .(theta_2))
  plot.ts(c(y_arch_train, test), main = title)
  abline(v = k_val, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/change_means/", "ma1_stoptime1.pdf", sep = "")
  pdf(file_name)
  
  
  plot.ts(c(y_arch_train, test), ylab = "MA(1)", main = title, xlab = "Index")
  abline(v = k_val, col = "red")
  dev.off()
  
  
}

#############################################
#
#############################################
plot_ARMA_sample_statistic_stop_time <- function(M, t,beta, c_alpha){
  phi <- 0.2
  theta <- 0.3
  window <- floor(t^(1/2))
  err = rnorm(M)
  y_arma1_train = arima.sim(model=list(ar=phi, ma = theta),n = M,innov=err)
  window_vec <- vector()
  M_mean <- mean(y_arma1_train)
  phi_2 <- 0.4
  theta_2 <- 0.5
  err_test = rnorm(t)
  test = arima.sim(model=list(ar=phi_2, ma = theta_2),n = t, innov=err_test)
  diff_vector <- vector()
  sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = "Parzen"))
  
  for(k in 1:(t-window)){
    sum =0 
    
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    if(abs(M_mean - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- k
      break
    }
  }
  
  k_val <- M + k
  title <- bquote("ARMA(1) Stopping Time," ~ phi == .(phi) ~ " "~ theta == .(theta) ~ ", Test " ~ phi == .(phi_2) ~ ", "~ theta == .(theta_2))
  plot.ts(c(y_arch_train, test), main = title)
  abline(v = k_val, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/change_means/", "arma1_stoptime1.pdf", sep = "")
  pdf(file_name)
  plot.ts(c(y_arch_train, test), ylab = "ARMA(1)", main = title, xlab = "Index")
  abline(v = k_val, col = "red")
  dev.off()
  
  
}

plot_arch_sample_statistic_stop_time(100, 100, 1, c_alphas_final_1[1])
plot_garch_sample_statistic_stop_time(100, 100, 1, c_alphas_final_1[1])
plot_AR_sample_statistic_stop_time(100, 100, 1, c_alphas_final_1[1])
plot_MA_sample_statistic_stop_time(100, 100, 1, c_alphas_final_1[1])
plot_ARMA_sample_statistic_stop_time(100, 100, 1, c_alphas_final_1[1])


###############################################
#GARCH
###############################################
plot_garch_stopping_time <- function(M, t, beta, c_alpha, kernel = "Bartlett"){
  h = floor(t^(1/2))
  max_vector <- vector() 
  for(i in 1:N){
    set.seed(i)
    err <- rnorm(M)
    y_ar1_train <- arima.sim(model=list(ar=0.2, ma=0.3),n = M,innov=err)
    window_vec <- vector()
    M_mean <- mean(y_ar1_train)
    err_test = rnorm(t)
    test = arima.sim(model=list(ar=0.2, ma=0.3),n = t, innov=err_test)
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(y_ar1_train, ker = kernel))
    for(k in 1:(t-h)){
      sum =0 
      for(j in k:(k+h)){
        sum = sum + test[j]
      }
      window_vec[k] <- sum/h
      
      
      g_T <- (sd_estimate*(h + k)^beta)/(h^{beta+1/2})
      
      diff_vector[k] <- (abs(M_mean - window_vec[k]))/(g_T)
    }
    
    max_vector[i] <- max(diff_vector)
  }
  
  alpha <- length(max_vector[max_vector >= c_alpha])/length(max_vector)
  
  return(alpha)
}


###############################################
#
###############################################
stock_sample_test_statistic_stop_time <- function(series_xts, M, beta, c_alpha, name, tick, ker, class){
  series <- as.ts(series_xts)
  stationary.test(series)
  series_df <- data.frame(index(series_xts), coredata(series_xts), stringsAsFactors = FALSE)
  colnames(series_df) <- c("Date", "Log Return")
  T_ <- length(series)-M
  h <- floor(sqrt(T_))
  sample <- series[1:M]
  test <- series[(M+1):length(series)]
  window <- h
  x_bar <- mean(sample)
  window_vec <- vector()
  diff_vector <- vector()
  stop_obs <- 0
  sd_estimate <- sqrt(asymp_var(sample, ker = ker))
  for(k in 1:(T_-h)){
    sum <- 0 
      
    for(j in k:(k+window)){
      sum <- sum + test[j]
    }
    window_vec[k] <- sum/window
      
      
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    diff <- abs(x_bar - window_vec[k])
    crit <- (g_T*c_alpha)
    if(diff >= crit){
      stop_obs <- M + k
      break
    }
  }
  
  plot.ts(series)
  abline(v = stop_obs, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name, "_stoptime.pdf", sep = "")

  title <- paste(tick, " Log ", class," Returns", sep="")
  pdf(file_name)

  
  plot(series, main = title, ylab = paste("Log ", class, " Return"))
  abline(v = stop_obs, col = "red")
  dev.off()
   
  #print(M)
  #print(T_)
  #print(k)
  #print(series_df[stop_obs, "Date"])
  
  return(c(stop_obs, series_df[stop_obs, "Date"]))
}

#################################################################################################################
#Stock Returns
#################################################################################################################



###############################################
#plot the time series of a given stock 
#series = time series of stock price
#tick = ticker symbol of the stock
#type = type of return being shown
###############################################
plot_stock_series<- function(series, tick, type){
  
  weekly_df <- data.frame(index(series), coredata(series), stringsAsFactors = FALSE)
  colnames(weekly_df) <- c("Date", "Log Return")
  
  file_name2 <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", tick, "_", type,  "_returns.pdf", sep = "")
  pdf(file_name2)
  title <- paste(tick, " Log ", type, " Returns", sep = "")
  plot(weekly_df, type = "l", main = title)
  dev.off()
  
}


calculate_returns_and_plots <- function(name, sample_pct = 0.10, weekly_sample_size = 100, toFile = FALSE, stop_time = FALSE, ker = "None"){

  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "2005-01-01", to = as.Date("2018-12-31"), env = NULL)
  weekly_returns <- weeklyReturn(Stock, subset = '2005::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '2005::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '2000::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '2000::', type = "log")
  
    plot_stock_series(weekly_returns, name, type = "Weekly")
    plot_stock_series(monthly_returns, name, type = "Monthly")
    plot_stock_series(quarterly_returns, name, type = "Quarterly")
    plot_stock_series(yearly_returns, name, type = "Annual")
  
    price_df <- data.frame(index(Stock), as.ts(Stock), stringsAsFactors = FALSE)
    colnames(price_df) <- c("Date", "Price (USD)")
    
    
    years <- unique(format(as.Date(price_df[,1]), "%Y"))
    
    
    
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    dev.off()

  if(toFile & stop_time){
    #Weekly Returns
    M <- weekly_sample_size
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 24
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly") 
    
    
    
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    #axis(1, at = years,labels = years)
    abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "blue")
    dev.off()

   
    M <- ceiling(sample_pct*length(quarterly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), tick = name,  ker, "Quarterly") 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), tick = name, ker, "Annual") 

    
    

    
  }

  
  
  
}


calculate_returns_and_plots_daily <- function(name, sample_pct = 0.20, weekly_sample_size = 100, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "2008-01-01", to = as.Date("2008-12-31"), env = NULL)
  daily_returns <- dailyReturn(Stock, subset = '2008::', type = "log")

  plot_stock_series(daily_returns, name, type = "daily")

  
  price_df <- data.frame(index(Stock), as.ts(Stock), stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  
  years <- unique(format(as.Date(price_df[,1]), "%Y"))
  
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_daily_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    #Weekly Returns
    M <- ceiling(sample_pct*length(daily_returns))
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_daily", sep = ""), tick = name, ker, "Daily") 
    #monthly returns

    print(stop_time)
    print(M)
    print(T_)
    print(window)
    print(price_df[stop_time[1],1])
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_daily_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    #axis(1, at = years,labels = years)
    abline(v = stop_time[2], col = "red")
    dev.off()
    
    
  }
  
  
  
  
}


calculate_returns_and_plots_alpha <- function(name, sample_pct = 0.10, weekly_sample_size = 100, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "2005-01-01", to = as.Date("2018-12-31"), env = NULL)
  weekly_returns <- weeklyReturn(Stock, subset = '2005::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '2005::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '2000::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '2000::', type = "log")
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  plot_stock_series(quarterly_returns, name, type = "Quarterly")
  plot_stock_series(yearly_returns, name, type = "Annual")
  
  price_df <- data.frame(index(Stock), as.ts(Stock), stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    #Weekly Returns
    M <- weekly_sample_size
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 24
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly") 
    
    #alpha = 0.05
    #Weekly Returns
    M <- weekly_sample_size
    stop_time_w_05 <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[2], paste(name, "_weekly_05", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 24
    stop_time_month_05 <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[2], paste(name, "_monthly_05", sep = ""), tick = name, ker, "Monthly") 
    
    #alpha = 0.05
    #Weekly Returns
    M <- weekly_sample_size
    stop_time_w_01 <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[2], paste(name, "_weekly_01", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 24
    stop_time_month_01 <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[2], paste(name, "_monthly_01", sep = ""), tick = name, ker, "Monthly") 
    
    
    
    
    
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime_full.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "blue")
    abline(v = stop_time_w_05[2], col = "green")
    abline(v = stop_time_month_05[2], col = "orange")
    abline(v = stop_time_w_01[2], col = "purple")
    abline(v = stop_time_month_01[2], col = "cyan")
    dev.off()
    
    
    #M <- ceiling(sample_pct*length(quarterly_returns))
    #stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), tick = name,  ker, "Quarterly") 
    
    #M <- ceiling(sample_pct*length(yearly_returns))
    #stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), tick = name, ker, "Annual") 
    
    
    
    
    
  }
  
  
  
  
}


calculate_returns_and_plots_dcom <- function(name, sample_pct = 0.10, weekly_sample_size = 104, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "1997-01-01", to = as.Date("2001-12-31"), env = NULL)
  weekly_returns <- weeklyReturn(Stock, subset = '1998::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '1998::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '1997::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '1997::', type = "log")
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  plot_stock_series(quarterly_returns, name, type = "Quarterly")
  plot_stock_series(yearly_returns, name, type = "Annual")
  
  price_df <- data.frame(index(Stock), as.ts(Stock), stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    
    #Weekly Returns
    M <- weekly_sample_size
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 24
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly") 
    
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "blue")
    dev.off()
    
    
    M <- ceiling(sample_pct*length(quarterly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), tick = name,  ker, "Quarterly") 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), tick = name, ker, "Annual") 
    
    
    
    
  }
  
  
  
  
}




calculate_returns_and_plots_neut <- function(name, sample_pct = 0.10, weekly_sample_size = 100, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "2005-01-01", to = as.Date("2018-12-31"), env = NULL)
  weekly_returns <- weeklyReturn(Stock, subset = '2005::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '2005::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '2005::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '2005::', type = "log")
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  plot_stock_series(quarterly_returns, name, type = "Quarterly")
  plot_stock_series(yearly_returns, name, type = "Annual")
  
  price_df <- data.frame(index(Stock), as.ts(Stock), stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    
    #Weekly Returns
    M <- weekly_sample_size
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 24
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly") 
    
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "blue")
    dev.off()
    
    
    M <- ceiling(sample_pct*length(quarterly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), tick = name,  ker, "Quarterly") 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), tick = name, ker, "Annual") 
    
    
    
    
    
  }
  
  
  
  
}



calculate_returns_and_plots_up <- function(name, sample_pct = 0.10, weekly_sample_size = 100, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "2005-01-01", to = as.Date("2018-12-31"), env = NULL)
  weekly_returns <- weeklyReturn(Stock, subset = '2014::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '2014::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '2014::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '2014::', type = "log")
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  plot_stock_series(quarterly_returns, name, type = "Quarterly")
  plot_stock_series(yearly_returns, name, type = "Annual")
  
  price_df <- data.frame(index(Stock), as.ts(Stock), stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    
    #Weekly Returns
    M <- weekly_sample_size
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 12
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly") 
    
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "blue")
    dev.off()
    
    
    M <- ceiling(sample_pct*length(quarterly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), tick = name,  ker, "Quarterly") 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), tick = name, ker, "Annual") 
    
    
    
  }
  
  
  
  
}

calculate_returns_and_plots_jap <- function(name, weeks = 30, months = 8, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols("^N225", from = "1900-01-01", to = as.Date("1992-12-31"), env = NULL)
  Stock <- na.omit(Stock)
  weekly_returns <- weeklyReturn(Stock, subset = '1988::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '1988::', type = "log")
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  
  price_df <- data.frame(index(Stock), as.ts(Stock), stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    
    #Weekly Returns
    M <- weeks
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- months
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly") 
    
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "blue")
    dev.off()
    
    
  }
  
  
  
  
}
calculate_returns_and_plots_SP <- function(name, sample_pct = 0.10, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- na.omit( 
    getSymbols(
      "SP500",
      src = "FRED",
      from = "1949-12-31",
      auto.assign = FALSE
    )
  )
  weekly_returns <- weeklyReturn(Stock, subset = '2005::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '2000::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '2000::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '2000::', type = "log")
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  plot_stock_series(quarterly_returns, name, type = "Quarterly")
  plot_stock_series(yearly_returns, name, type = "Annual")
  
  
  
  if(toFile & stop_time){
    
    M <- ceiling(sample_pct*length(weekly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), ker) 
    
    M <- ceiling(sample_pct*length(monthly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(monthly_returns, M, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), ker) 
    
    M <- ceiling(sample_pct*length(quarterly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), ker) 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), ker) 
    
    
    
    
    
  }
  
  
  
  
}



file_name1 <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", "AAPL",  "_returns.pdf", sep = "")
pdf(file_name1)
plot.xts(weekly_returns, ylab = "Log Return")
dev.off()


options("getSymbols.warning4.0"=FALSE)
Stock <- getSymbols("JPM", from = "1900-01-01", to = as.Date("2018-12-31"), env = NULL)
Stock <- na.omit(Stock)

plot_stock_price(Stock, "JPM")
weekly_returns <- weeklyReturn(Stock, subset = '2005::')
stock_sample_test_statistic_stop_time(weekly_returns, 73, 1, c_alphas_final_1[1], "JPM", ker = "None")
plot(monthly_returns)
monthly_returns <- monthlyReturn(Stock, subset = '2000::')
quarterly_returns <- quarterlyReturn(Stock, subset = '2000::', type = "log")
yearly_returns <- annualReturn(Stock, subset = '1990::', type = "log")

calculate_returns_and_plots("AAPL", sample_pct = 0.50, toFile = TRUE, stop_time = TRUE, ker = "Flat-Top")


calculate_returns_and_plots("BP", sample_pct = 0.25, toFile = TRUE, stop_time = TRUE, ker = "Flat-Top")


calculate_returns_and_plots("BAC", toFile = TRUE, stop_time = TRUE, ker = "Flat-Top")


calculate_returns_and_plots("PYS", toFile = TRUE, stop_time = TRUE, ker = "Flat-Top")

calculate_returns_and_plots_alpha("AIG", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("AIG", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("BAC", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("JPM", sample_pct = 0.10, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("C", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("AXP", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("GS", sample_pct = 0.1,toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("WMT", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("DLTR", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("GE", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("ZION", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots("NVDA", toFile = TRUE, stop_time = TRUE, ker = "None")

calculate_returns_and_plots_dcom("QCOM", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_dcom("CSCO", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_dcom("MU", sample_pct = 0.1, weekly_sample_size = 100, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_dcom("MSFT", sample_pct = 0.1, weekly_sample_size = 100, toFile = TRUE, stop_time = TRUE, ker = "None")

calculate_returns_and_plots_neut("PG", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_neut("CLX", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_neut("JNJ", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_neut("AZO", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")

calculate_returns_and_plots_up("GPRO", sample_pct = 0.1, weekly_sample_size = 24, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_up("P", sample_pct = 0.1, toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_returns_and_plots_up("FIT", toFile = TRUE, weekly_sample_size = 24, stop_time = TRUE, ker = "None")


calculate_returns_and_plots_jap("Nikkei", months = 12, weeks = 52, toFile = TRUE, stop_time = TRUE, ker = "None")

calculate_returns_and_plots_daily("AAPL", sample_pct = .30, toFile = TRUE, stop_time = TRUE, ker = "Quadratic Spectral")

########################################
#Bitcoin
########################################




calculate_returns_and_plots_crypto <- function(name, sample_pct = 0.10, weekly_sample_size = 104, toFile = FALSE, stop_time = FALSE, ker = "None"){

  Stock <- read.csv("data/btc.csv", header = TRUE)
  Stock_clean <- na.omit(Stock)
  coin_xts <- as.xts(Stock_clean$Close, order.by = as.Date(as.character(Stock_clean$Date), format = "%m/%d/%y"))
  weekly_returns <- na.omit(weeklyReturn(coin_xts, type = "log"))
  monthly_returns <- na.omit(monthlyReturn(coin_xts, type = "log"))
  quarterly_returns <- na.omit(quarterlyReturn(coin_xts, type = "log"))
  yearly_returns <- na.omit(annualReturn(coin_xts, type = "log"))
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  plot_stock_series(quarterly_returns, name, type = "Quarterly")
  plot_stock_series(yearly_returns, name, type = "Annual")
  
  price_df <- data.frame(as.Date(as.character(Stock_clean$Date),format = "%m/%d/%y"), Stock_clean$Close, stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    
    #Weekly Returns
    M <- weekly_sample_size
    stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    #monthly returns
    M_mo <- 24
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly") 
    
    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Stock Price", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "blue")
    dev.off()
    
    
    M <- ceiling(sample_pct*length(quarterly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), tick = name,  ker, "Quarterly") 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), tick = name, ker, "Annual") 
    
    
    
    
  }
  
  
  
  
}


calculate_returns_and_plots_crypto("BTC", toFile = TRUE, stop_time = TRUE, ker = "None")



##################################################
#External data sets
##################################################


calculate_returns_and_plots_stoptime <- function(path, name, weeks = 100, months = 12 , toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  Stock <- read.csv(path, header = TRUE)
  Stock_clean <- na.omit(Stock)
  coin_xts <- as.xts(Stock_clean$Close, order.by = as.Date(as.character(Stock_clean$Date), format = "%m/%d/%Y"))
  #weekly_returns <- na.omit(weeklyReturn(coin_xts, type = "log"))
  monthly_returns <- na.omit(monthlyReturn(coin_xts, type = "log"))
  #quarterly_returns <- na.omit(quarterlyReturn(coin_xts, type = "log"))
  #yearly_returns <- na.omit(annualReturn(coin_xts, type = "log"))
  
  #plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  #plot_stock_series(quarterly_returns, name, type = "Quarterly")
  #plot_stock_series(yearly_returns, name, type = "Annual")
  
  price_df <- data.frame(as.Date(as.character(Stock_clean$Date),format = "%m/%d/%Y"), Stock_clean$Close, stringsAsFactors = FALSE)
  colnames(price_df) <- c("Date", "Price (USD)")
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price.pdf", sep = "")
  pdf(file_name)
  title <- paste(name, " Stock Price", sep = "")
  plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
  dev.off()
  
  if(toFile & stop_time){
    
    # #Weekly Returns
    # M <- weeks
    # stop_time <- stock_sample_test_statistic_stop_time(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), tick = name, ker, "Weekly") 
    # #monthly returns
    M_mo <- months
    stop_time_month <- stock_sample_test_statistic_stop_time(monthly_returns, M_mo, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""), tick = name, ker, "Monthly")

    #Construct Plot
    file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", name,"_price_stoptime.pdf", sep = "")
    pdf(file_name)
    title <- paste(name, " Price Index", sep = "")
    plot(price_df[,1], price_df[,2], type = "l", xlab = "Date", ylab = "Price (USD)", main = title)
    #abline(v = stop_time[2], col = "red")
    abline(v = stop_time_month[2], col = "red")
    dev.off()
    
    
    # M <- ceiling(sample_pct*length(quarterly_returns))
    # stop_time <- stock_sample_test_statistic_stop_time(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""), tick = name,  ker, "Quarterly") 
    # 
    # M <- ceiling(sample_pct*length(yearly_returns))
    # stop_time <- stock_sample_test_statistic_stop_time(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""), tick = name, ker, "Annual") 
    # 
    # 
    # 
    
  }
  
  
  
  
}
##################################################
#Perform sequential procedure on great depression data
##################################################
calculate_returns_and_plots_stoptime("data/DJIA_data_1929.csv", "DJIA", toFile = TRUE, stop_time = TRUE, ker = "None")



##############################################################################################
#All Stop Times
##############################################################################################


stock_sample_stop_times <- function(series_xts, M, beta, c_alpha, file_name, stock, ker){
  series <- as.ts(series_xts)
  series_df <- data.frame(index(series_xts), coredata(series_xts), stringsAsFactors = FALSE)
  colnames(series_df) <- c("Date", "Log Return")
  T_ <- length(series)-M
  h <- floor(T_^(1/2))
  sample <- series[1:M]
  test <- series[(M+1):length(series)]
  window <- h
  x_bar <- mean(sample)
  window_vec <- vector()
  diff_vector <- vector()
  stop_obs <- vector()
  sd_estimate <- sqrt(asymp_var(sample, ker = ker))
  for(k in 1:(T_-h)){
    sum <- 0 
    
    for(j in k:(k+window)){
      sum <- sum + series[j]
    }
    window_vec[k] <- sum/window
    
    
    g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
    if(abs(x_bar - window_vec[k]) >= (g_T*c_alpha)){
      stop_obs <- c(stop_obs, k)
    }
  }
  
  plot.ts(series)
  abline(v = stop_obs, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", file_name, "_stoptimes.pdf", sep = "")
  if(length(stop_obs) != 0){
    title <- paste(stock, " M = ", M, ", Stops = ", length(stop_obs), sep = "")
  }
  else{
    title <- paste("K = 0"," M = ", M,  sep = "")
  }
  pdf(file_name)
  
  
  plot(series, main = title, ylab = "Log Return")
  abline(v = stop_obs, col = "red")
  dev.off()
  
  
  return(series_df[stop_obs, "Date"])
}



############################################
#Calculates and plots all stop times that are found with the sequential monitoring procedure.
############################################

calculate_all_stop_times <- function(name, sample_pct = 0.10, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "1900-01-01", to = as.Date("2018-12-31"), env = NULL)
  Stock <- na.omit(Stock)
  weekly_returns <- weeklyReturn(Stock, subset = '2005::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '2000::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '2000::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '1990::', type = "log")
  
  #plot_stock_series(weekly_returns, name, type = "Weekly")
  #plot_stock_series(monthly_returns, name, type = "Monthly")
  #plot_stock_series(quarterly_returns, name, type = "Quarterly")
  #plot_stock_series(yearly_returns, name, type = "Annual")
  
  
  
  if(toFile & stop_time){
    
    M <- ceiling(sample_pct*length(weekly_returns))
    week_stop_time <- stock_sample_stop_times(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), name, ker) 
    
    M <- ceiling(sample_pct*length(monthly_returns))
    month_stop_time <- stock_sample_stop_times(monthly_returns, M, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""),name, ker) 
    
    M <- ceiling(sample_pct*length(quarterly_returns))
    quarter_stop_time <- stock_sample_stop_times(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""),name,ker) 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    year_stop_time <- stock_sample_stop_times(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""),name, ker) 
    
    
    
    
    
  }
  
  
  ###Construct table of stop times
  
}

calculate_all_stop_times("AIG", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_all_stop_times("BAC", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_all_stop_times("JPM", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_all_stop_times("C", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_all_stop_times("AXP", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_all_stop_times("GS", toFile = TRUE, stop_time = TRUE, ker = "None")



############################################
#Rolling window approach
############################################

stock_sample_rolling_stop_times <- function(series_xts, M, beta, c_alpha, file_name, stock, ker){
  series <- as.ts(series_xts)
  series_df <- data.frame(index(series_xts), coredata(series_xts), stringsAsFactors = FALSE)
  colnames(series_df) <- c("Date", "Log Return")
  T_ <- 100
  h <- floor(T_^(1/2))
  sample <- series[1:M]
  window <- h
  x_bar <- mean(sample)
  window_vec <- vector()
  diff_vector <- vector()
  stop_obs <- vector()
  sd_estimate <- sqrt(asymp_var(sample, ker = ker))
  while(M +T_ < length(series)){
    sample <- series[1:M]
    x_bar <- mean(sample)
    window_vec <- vector()
    diff_vector <- vector()
    sd_estimate <- sqrt(asymp_var(sample, ker = ker))
  for(k in (M+1):(M+T_-window)){
    sum <- 0 
    
    for(j in k:(k+window)){
      sum <- sum + series[j]
      window_vec[k] <- sum/window
      
      
      g_T <- (sd_estimate*(window + k)^beta)/(window^{beta+ 1/2})
      #print(abs(x_bar - window_vec[k]))
      #print((g_T*c_alpha))
      if(abs(x_bar - window_vec[k]) >= (g_T*c_alpha)){
        stop_obs <- c(stop_obs, k)
      }
    }
    
  }
    M <- M + T_
  }
  
  plot.ts(series)
  abline(v = stop_obs, col = "red")
  
  
  file_name <- paste("/Users/Jon/Library/Mobile Documents/com~apple~CloudDocs/Class Work 2018/Thesis/plots/stocks/", file_name, "_rolling_stoptimes.pdf", sep = "")
  if(length(stop_obs) != 0){
    title <- paste(stock, ", Stops =  ", length(stop_obs), sep = "")
  }
  else{
    title <- paste("K = 0"," Stops = ", 0,  sep = "")
  }
  pdf(file_name)
  
  
  plot(series, main = title, ylab = "Log Return")
  abline(v = stop_obs, col = "red")
  dev.off()
  
  
  return(series_df[stop_obs, "Date"])
}




calculate_sequential_stop_times <- function(name, sample_pct = 0.10, toFile = FALSE, stop_time = FALSE, ker = "None"){
  
  options("getSymbols.warning4.0"=FALSE)
  Stock <- getSymbols(name, from = "1900-01-01", to = as.Date("2011-12-31"), env = NULL)
  Stock <- na.omit(Stock)
  weekly_returns <- weeklyReturn(Stock, subset = '2000::', type = "log")
  monthly_returns <- monthlyReturn(Stock, subset = '2000::', type = "log")
  quarterly_returns <- quarterlyReturn(Stock, subset = '2000::', type = "log")
  yearly_returns <- annualReturn(Stock, subset = '1990::', type = "log")
  
  plot_stock_series(weekly_returns, name, type = "Weekly")
  plot_stock_series(monthly_returns, name, type = "Monthly")
  plot_stock_series(quarterly_returns, name, type = "Quarterly")
  plot_stock_series(yearly_returns, name, type = "Annual")
  
  
  
  if(toFile & stop_time){
    
    M <- ceiling(sample_pct*length(weekly_returns))
    week_stop_time <- stock_sample_rolling_stop_times(weekly_returns, M, 1, c_alphas_final_1[1], paste(name, "_weekly", sep = ""), name, ker) 
    
    M <- ceiling(sample_pct*length(monthly_returns))
    month_stop_time <- stock_sample_rolling_stop_times(monthly_returns, M, 1, c_alphas_final_1[1], paste(name, "_monthly", sep = ""),name, ker) 
    
    M <- ceiling(sample_pct*length(quarterly_returns))
    quarter_stop_time <- stock_sample_rolling_stop_times(quarterly_returns, M, 1, c_alphas_final_1[1], paste(name, "_quarterly", sep = ""),name,ker) 
    
    M <- ceiling(sample_pct*length(yearly_returns))
    year_stop_time <- stock_sample_rolling_stop_times(yearly_returns, M, 1, c_alphas_final_1[1], paste(name, "_annual", sep = ""),name, ker) 
    
    
    
    
    
  }
  
  
  ###Construct table of stop times
  
}


calculate_sequential_stop_times("AAPL", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_sequential_stop_times("AIG", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_sequential_stop_times("BAC", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_sequential_stop_times("JPM", toFile = TRUE, stop_time = TRUE, ker = "Flat-Top")
calculate_sequential_stop_times("C", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_sequential_stop_times("AXP", toFile = TRUE, stop_time = TRUE, ker = "None")
calculate_sequential_stop_times("GS", toFile = TRUE, stop_time = TRUE, ker = "None")


citation(package = "quantmod")
citation(package = "fGarch")
