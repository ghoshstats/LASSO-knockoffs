generate_nonzero_vector <- function(p, num_nonzero=30, min_value = -2, max_value = 3) {
  # Initialize vector with zeros
  vec <- numeric(num_nonzero)
  # Loop until we have filled num_nonzero non-zero entries
  while (sum(vec != 0) < num_nonzero) {
    
    # Generate random values from U(min_value, max_value)
    rand_values <- runif(num_nonzero, min = min_value, max = max_value)
    
    # Assign non-zero random values to vec
    vec[1:30] <- rand_values
    
    # Check if any zero values were assigned (should not happen due to runif properties)
    if (any(vec == 0)) {
      vec[vec == 0] <- runif(sum(vec == 0), min = min_value, max = max_value)
    }
  }
  vec2<-rep(0,p-num_nonzero)
  return(c(vec, vec2))
}

generateOmega <- function(q, type = "AR1", params = NULL) {
  Omega <- matrix(0, nrow = q, ncol = q)
  
  if (type == "AR1") {
    rho <- params$rho
    Sigma <- rho ^ (abs(row(Omega) - col(Omega)))
    Omega <- solve(Sigma)
    
  } else if (type == "AR2") {
    rho <- params$rho
    Omega <- rho ^ (abs(row(Omega) - col(Omega))) * (abs(row(Omega) - col(Omega)) <= 2)
    Sigma <- solve(Omega)
    return(list(Sigma = Sigma, Omega = Omega))
    
  } else if (type == "BlockDiagonal") {
    rho <- params$rho
    # Assuming two blocks for simplicity
    Sigma <- matrix(0, nrow = q, ncol = q)
    block_size <- q / 2
    Sigma[1:block_size, 1:block_size] <- rho
    Sigma[(block_size+1):q, (block_size+1):q] <- rho
    diag(Sigma) <- 1
    Omega <- solve(Sigma)
    
  } else if (type == "StarGraph") {
    rho <- params$rho
    Omega <- diag(1, q, q)
    Omega[2:q, 1] <- rho
    Omega[1, 2:q] <- rho
    Sigma <- solve(Omega)
    return(list(Sigma = Sigma, Omega = Omega))
    
  } else if (type == "SmallWorld") {
    rhos <- params$rhos
    Omega <- diag(rhos[1], q, q)
    Omega[abs(row(Omega)-col(Omega))==1] <- rhos[2]
    Omega[1, q] <- Omega[q, 1] <- rhos[3]
    Sigma <- solve(Omega)
    return(list(Sigma = Sigma, Omega = Omega))
    
  } else if (type == "TreeNetwork") {
    rhos <- params$rhos
    Omega <- matrix(rhos[2], q, q)
    diag(Omega) <- rhos[1]
    Sigma <- solve(Omega)
    return(list(Sigma = Sigma, Omega = Omega))
    
  } else {
    stop("Unknown type")
  }
  
  # For types that compute Sigma within the conditional
  return(list(Sigma = solve(Omega), Omega = Omega))
}

mean_FDR_TPR <- function(n, N, p, type, params = NULL, FDR) {
  beta <- generate_nonzero_vector(p)  # Assuming this function generates a non-zero vector of length p
  FDR_values <- numeric(N)  # Initialize an empty numeric vector to store FDR values
  TPR_values <- numeric(N)  # Initialize an empty numeric vector to store TPR values
  
  for (i in 1:N) {
    U <- generateOmega(p, type, params)[["Sigma"]]  # Generate covariance matrix or similar
    X <- mvtnorm::rmvnorm(n, rep(0, p), U)  # Generate multivariate normal data
    
    # Generate the response from a linear model
    eta <- X %*% beta
    pi <- 1 / (1 + exp(-eta))
    y.sample <- rbinom(n, 1, pi)
    
    # Apply knockoff procedure
    library(knockoff)
    result <- knockoff.filter(X, y.sample, fdr = FDR, offset = 0)
    
    # Compute False Discovery Rate (FDR)
    indices_outside_first_30 <- sum(!(result$selected %in% 1:30))
    FDR_value <- indices_outside_first_30 / length(result$selected)
    FDR_values[i] <- FDR_value  # Store FDR value for this iteration
    
    # Compute True Positive Rate (TPR)
    true_positives <- sum(result$selected %in% 1:30)  # Count selected indices within 1 to 30
    TPR_value <- true_positives / 30  # Assuming 30 as the number of true signals
    TPR_values[i] <- TPR_value  # Store TPR value for this iteration
  }
  
  mean_FDR <- mean(FDR_values,na.rm = TRUE)  # Compute the mean FDR over all iterations
  mean_TPR <- mean(TPR_values,na.rm = TRUE)  # Compute the mean TPR over all iterations
  
  result <- data.frame(FDR = mean_FDR, TPR = mean_TPR)
  return(result)
}

#0.01
l=c(0.01,0.05,0.10,0.15,0.2)
results_list <- list()
for (i in 1:2){
  k<-mean_FDR_TPR(200,100,400,type="AR1",params=list(rho=0.8),FDR=l[i])
  results_list[[as.character(i)]] <- k
}
results_list
results_df <- do.call(rbind, results_list)
results_df
write.csv(results_df,file="FDR_AR1_LR.csv", row.names = TRUE) 
results_df
read.csv("FDR_AR1_LR.csv")


l=c(0.01,0.05,0.10,0.15,0.2)
results_list <- list()
for (i in 1:length(l)){
  k<-mean_FDR_TPR(200,100,400,type="AR1",params=list(rho=0.8),FDR=l[i])
  results_list[[as.character(i)]] <- k
}
results_df <- do.call(rbind, results_list)
write.csv(results_df,file="FDR_AR1_LR.csv", row.names = TRUE) 
results_df
read.csv("FDR_AR1_LR.csv")

l=c(0.01,0.05,0.10,0.15,0.2)
results_list <- list()
for (i in 1:length(l)) {
  k<-mean_FDR_TPR(200,100,400,type="BlockDiagonal",params=list(rho=0.8),FDR=l[i])
  results_list[[as.character(i)]] <- k
}
results_df <- do.call(rbind, results_list)
write.csv(results_df,file="FDR_BlockDiagonal_LR.csv", row.names = TRUE)

l=c(0.01,0.05,0.10,0.15,0.2)
results_list <- list()
for (i in 1:length(l)) {
  k<-mean_FDR_TPR(200,100,400,type="SmallWorld",params=list(rhos=c(0.8,0.7,0.6)),FDR=l[i])
  results_list[[as.character(i)]] <- k
}
results_df <- do.call(rbind, results_list)
write.csv(results_df,file="FDR_SmallWorld_LR.csv", row.names = TRUE)

#0.05,0.15
l=c(0.05,0.15)
results_list <- list()
for (i in 1:length(l)) {
  k<-mean_FDR_TPR(200,100,400,type="StarGraph",params=list(rho=0.8),FDR=l[i])
  results_list[[as.character(i)]] <- k
}
results_list
results_df <- do.call(rbind, results_list)
results_df
write.csv(results_df,file="FDR_StarGraph_LR.csv", row.names = TRUE)


l=c(0.01,0.05,0.10,0.15,0.2)
results_list <- list()
for (i in 1:length(l)) {
  k<-mean_FDR_TPR(200,100,400,type="StarGraph",params=list(rho=0.8),FDR=l[i])
  results_list[[as.character(i)]] <- k
}
results_df <- do.call(rbind, results_list)
write.csv(results_df,file="FDR_StarGraph_LR.csv", row.names = TRUE)