normal_lkh <- function(x, mu, sigma){
  scale <- (1 / (sqrt(2 * pi * sigma)))
  quad <- (x - mu)^2 / (2 * sigma)
  likelihood <- scale * exp(-quad)
  return(as.numeric(likelihood))
}



mix_normal_lkh <- function(z, s, sigma, Pi, lambda){
  lkh_bind <- function(x){
    return(normal_lkh(x[1], 0, x[2]^2+sigma[k]^2))
  }
  p <- length(Pi)
  N <- length(z)
  lkh_all <- NULL
  for (k in 1:p){
    lkh <- apply(cbind(z, s), 1, lkh_bind)
    lkh_temp <- tibble(i=seq(1, N, 1), lkh=lkh, Pi=Pi[k], 
                       lambda=lambda[k], k=k-1) %>% relocate(k)
    lkh_all <- rbind(lkh_all, lkh_temp)
  }
  lkh_all <- lkh_all %>% mutate(lkh_Pi = lkh * Pi) 
  return(lkh_all)
}



log_lkh <- function(lkh){
  log_lkh_idv <- lkh %>% group_by(i) %>% summarize(log_sum = log(sum(lkh_Pi))) 
  return(sum(log_lkh_idv$log_sum))
}



post_mean <- function(z, s, Pi, sigma){
  Pi[Pi < 1e-10] <- 0
  nume <- z * Pi * sigma^2
  demo <- sigma^2 + s^2
  pos_mean <- sum(nume / demo)
  return(pos_mean)
}



post_mean_bind <- function(x, Pi, sigma){
  return(post_mean(x[1], x[2], Pi, sigma^2))
}



ash_em <- function(z, s, iter=500, verbose=FALSE, return='list'){

  N <- length(z)

  # determine the grid of sigma
  sigma_min <- min(s) / 10
  if (max(z^2 - s^2) > 0){
    sigma_max <- 2 * sqrt(max(z^2 - s^2))
  }else{
    sigma_max = 8 * sigma_min
  }
  sigma_grid <- sigma_min * sqrt(2)^(seq(0,100,1))
  sigma_grid <- sigma_grid[sigma_grid <= sigma_max]
  sigma_grid <- c(sigma_grid, sigma_max)
  sigma_grid <- c(0, sigma_grid)
  K <- length(sigma_grid) - 1
  
  # initialize Pi and Lambda
  Pi <- rep(1/N, K)
  Pi <- c(1-sum(Pi), Pi)
  Pi <- abs(Pi) / sum(abs(Pi)) 
  lambda <- c(10, rep(1, K))

  # create vector/matrix to store updates of parameters
  Pi_mat <- matrix(0, nrow=iter, ncol=K+1)
  Pi_mat[1,] <- Pi
  log_lkh_list <- rep(0, iter)
  penalty <- rep(0, iter)

  # EM algorithm
  for (i in 2:iter){
    lkh <- mix_normal_lkh(z=z, Pi=Pi_mat[i-1, ], s=s, sigma=sigma_grid, lambda=lambda)
    if (i == 2){
      log_lkh_list[1] <- log_lkh(lkh)
      penalty[1] <- sum((lambda-1) * log(Pi_mat[1, ]))
    }

    # update Pi
    nume_Pi <- lkh %>% group_by(i) %>% mutate(demo = sum(lkh_Pi)) %>%
      ungroup %>% mutate(gamma = lkh_Pi / demo) %>% group_by(k) %>%
      mutate(a = sum(gamma) + lambda - 1) %>% summarize(a = mean(a))
    demo_Pi <- nume_Pi$a %>% sum()
    Pi_mat[i, ] <- nume_Pi$a / demo_Pi
    Pi_mat[i, ][Pi_mat[i, ] < 1e-18] <- 1e-18

    # evalute the incomplete log-likelihood using updated Pi
    lkh <- mix_normal_lkh(z=z, Pi=Pi_mat[i, ], s=s, sigma=sigma_grid, lambda=lambda)
    log_lkh_list[i] <- log_lkh(lkh)
    penalty[i] <- sum((lambda-1) * log(Pi_mat[i, ]))

    # check the stop criterion
    tot <- log_lkh_list[i] + penalty[i]  - log_lkh_list[i-1] - penalty[i-1]
    if (tot < 1e-8) break

    if (verbose){
      cat(sprintf('Iteration %s, the penalized incomplete Log-Likelihood = %s', i, log_lkh_list[i]+ penalty[i]))
      cat('\n')
    }
  }

  # compute the posterior mean
  pos_mu <- apply(cbind(z, s^2), 1, post_mean_bind, Pi=Pi_mat[i, ], sigma=sigma_grid)

  # output
  if (return == 'list'){
    out <- list(pos_mu=pos_mu, log_lkh=log_lkh_list[1:i], Pi=Pi_mat[i, ], iters = i)
  }else if (return == 'mu'){
    out <- as.numeric(pos_mu)
  }
  
  return(out)
}
