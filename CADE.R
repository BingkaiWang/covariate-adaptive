CADE <- function(Y, A, Strata = NULL, W = NULL, tau = 0, pi = 0.5, family = "gaussian"){
  if(any(is.na(Y))){
    stop("No missing outcome is allowed.")
  }
  if(any(is.na(A))){
    stop("No missing treatment is allowed.")
  }
  if(is.null(Strata)){
    message("No strata variable is given. Treatment is assumed to be assigned randomly.")
    if(is.null(W)){
      message("No covarites is given. The unadjusted estimator is calculated.")
      est <- mean(Y[A==1]) - mean(Y[A==0])
      sd <- sqrt(var(Y[A==1])/pi + var(Y[A==0])/(1-pi))
    }else{
      message("The estimator adjusting for W is calculated.")
      d <- data.frame(Y, A, W)
      d1 <- data.frame(Y, A = 1, W)
      d0 <- data.frame(Y, A = 0, W)
      glm.fit <- glm(Y~ ., data = d, family = family)
      p1 <- predict(glm.fit, d1, type = "response")
      p0 <- predict(glm.fit, d0, type = "response")
      est <- mean(p1) - mean(p0)
      Y.1 <- Y[A==1]
      Y.0 <- Y[A==0]
      r.1 <- Y.1 - p0[A==1] * pi - p1[A==1] * (1 - pi)
      r.0 <- Y.0 - p0[A==0] * pi - p1[A==0] * (1 - pi)
      sd <- sqrt(var(r.1)/pi + var(r.0)/(1-pi))
    }
    return(c(est = est, sd = sd))
  }else{
    if(is.null(W)){
      message("No covarites is given. The unadjust estimator is calculated.")
      d <- data.frame(Y, A)
      d1 <- data.frame(Y, A = 1)
      d0 <- data.frame(Y, A = 0)
      glm.fit <- glm(Y~ ., data = d, family = family)
    } else {
      message("The estimator adjusting for Strata and W is calculated.")
      d <- data.frame(Y, A, Strata, W)
      d1 <- data.frame(Y, A = 1, Strata, W)
      d0 <- data.frame(Y, A = 0, Strata, W)
      glm.fit <- glm(Y~ 0 + ., data = d, family = family)
    }
  }
  p1 <- predict(glm.fit, d1, type = "response")
  p0 <- predict(glm.fit, d0, type = "response")
  est <- mean(p1) - mean(p0)
  # estimate standard error
  num.strata <- ncol(Strata)
  n <- length(Y)
  n.1 <- sum(A) 
  n.0 <- n - n.1
  Y.1 <- Y[A==1]
  Y.0 <- Y[A==0]
  r.1 <- Y.1 - p0[A==1] * pi - p1[A==1] * (1 - pi)
  r.0 <- Y.0 - p0[A==0] * pi - p1[A==0] * (1 - pi)
  tilde.r.1 <- NA
  tilde.r.0 <- NA
  Er1_S <- rep(NA, num.strata)
  Er0_S <- rep(NA, num.strata)
  p.s <- rep(NA, num.strata)
  for(s in 1:num.strata){
    index.1 <- which(I.S[A == 1,s] == 1)
    index.0 <- which(I.S[A == 0,s] == 1)
    Er1_S[s] <- mean(r.1[index.1])
    Er0_S[s] <- mean(r.0[index.0])
    p.s[s] <- sum(I.S[,s] == 1)/n
    tilde.r.1[index.1] <- r.1[index.1] - mean(r.1[index.1])
    tilde.r.0[index.0] <- r.0[index.0] - mean(r.0[index.0])
  }
  V.1 <- var(tilde.r.1)/pi + var(tilde.r.0)/(1 - pi)
  V.2 <- sum(p.s * ((Er1_S - mean(r.1)) - (Er0_S - mean(r.0)))^2)
  V.3 <- tau * sum(p.s * ((Er1_S - mean(r.1))/pi + (Er0_S - mean(r.0))/(1 - pi))^2)
  sd <- sqrt(V.1 + V.2 + V.3)
  return(c(est = est, sd = sd))
  
}