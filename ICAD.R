# Inference under Covariate-Adaptive Design (ICAD) with continuous or binary outcomes
ICAD <- function(Y, A, Strata, W = NULL, pi = 0.5, family = "gaussian"){
  if(any(is.na(A))){
    stop("No missing treatment is allowed.")
  }
  if(any(is.na(Strata))){
    stop("Strata must be provided for asymptotic efficiency.")
  }
  if(!is.null(W) & any(is.na(W))){
    stop("Please impute missing covariates first.") 
  }
  if(!any(is.na(Y))){
    message("No missing outcome presents. The unadjusted estimator and adjusted estimator are cacluated.")
    return(rbind(
      unadjsed = eCAD(Y, A, Strata, W, pi, family, method = "unadjusted"),
      adjusted = eCAD(Y, A, Strata, W, pi, family, method = "adjusted")
    ))
  }else{
    message("There are missing outcomes. The unadjusted estimator, adjusted estimator and DR-WLS estimator are calculated.")
    return(rbind(
      unadjsed = eCAD(Y[!is.na(Y)], A[!is.na(Y)], Strata[!is.na(Y)], W[!is.na(Y),], pi, family, method = "unadjusted"),
      adjusted = eCAD(Y[!is.na(Y)], A[!is.na(Y)], Strata[!is.na(Y)], W[!is.na(Y),], pi, family, method = "adjusted"),
      drwls = eCAD(Y, A, Strata, W, pi, family, method = "DR-WLS")
    )) 
  }
}

eCAD <- function(Y, A, Strata, W = NULL, pi = 0.5, family = "gaussian", method){
  if(method == "unadjusted"){
    d <- data.frame(Y, A, Strata)
    n <- nrow(d)
    d1 <- data.frame(Y, A = 1, Strata)
    d0 <- data.frame(Y, A = 0, Strata)
    glm.fit <- glm(Y~ ., data = d, family = family)
    p1 <- predict(glm.fit, d1, type = "response")
    p0 <- predict(glm.fit, d0, type = "response")
    est <- mean(p1) - mean(p0)
    Y.1 <- Y[A==1]
    Y.0 <- Y[A==0]
    r <- Y - predict(glm.fit, d, type = "response")
    r.1 <- Y.1 - p0[A==1] * pi - p1[A==1] * (1 - pi)
    r.0 <- Y.0 - p0[A==0] * pi - p1[A==0] * (1 - pi)
    tilde_var <- var(r.1)/pi + var(r.0)/(1-pi)
    E.S <- map_dbl(unique(Strata), ~mean(((2*A-1)*r)[Strata == .]))
    var <- tilde_var - (1-2*pi)^2/pi/(1-pi) * var(E.S)
    return(c(est = est, var = var/n, 
           CI.lower = qnorm(0.025, mean = est, sd =  sqrt(var/n)),
           CI.upper = qnorm(0.975, mean = est, sd =  sqrt(var/n))))
  }
  if(method == "adjusted"){
    d <- data.frame(Y, A, Strata, W)
    n <- nrow(d)
    d1 <- data.frame(Y, A = 1, Strata, W)
    d0 <- data.frame(Y, A = 0, Strata, W)
    glm.fit <- glm(Y~ ., data = d, family = family)
    p1 <- predict(glm.fit, d1, type = "response")
    p0 <- predict(glm.fit, d0, type = "response")
    est <- mean(p1) - mean(p0)
    Y.1 <- Y[A==1]
    Y.0 <- Y[A==0]
    r <- Y - predict(glm.fit, d, type = "response")
    r.1 <- Y.1 - p0[A==1] * pi - p1[A==1] * (1 - pi)
    r.0 <- Y.0 - p0[A==0] * pi - p1[A==0] * (1 - pi)
    tilde_var <- var(r.1)/pi + var(r.0)/(1-pi)
    E.S <- map_dbl(unique(Strata), ~mean(((2*A-1)*r)[Strata == .]))
    var <- tilde_var - (1-2*pi)^2/pi/(1-pi) * var(E.S)
    return(c(est = est, var = var/n, 
           CI.lower = qnorm(0.025, mean = est, sd =  sqrt(var/n)),
           CI.upper = qnorm(0.975, mean = est, sd =  sqrt(var/n))))
  }
  if(method == "DR-WLS"){
    d <- data.frame(Y, A, Strata, W)
    n <- nrow(d)
    d1 <- data.frame(Y, A = 1, Strata, W)
    d0 <- data.frame(Y, A = 0, Strata, W)
    M <- !is.na(Y)
    propensity.fit <- glm(M ~ ., data = data.frame(M, A, Strata, W), family = "binomial")
    propensity_score <- predict(propensity.fit, type = "response")
    glm.fit <- glm(Y~ ., data = d, family = family, weights = 1/propensity_score)
    X <- model.matrix(propensity.fit)
    p <- ncol(X)
    n <- nrow(d)
    pA <- predict(glm.fit, d, type = "response")
    p1 <- predict(glm.fit, d1, type = "response")
    p0 <- predict(glm.fit, d0, type = "response")
    est <- mean(p1) - mean(p0)
    if(family == "gaussian"){
      h_beta_A <- h_beta_0 <- h_beta_1 <- X
      h_beta_1[,"A"] <- 1
      h_beta_0[,"A"] <- 0
    }else if(family == "binomial"){
      h_beta_A <- h_beta_0 <- h_beta_1 <- (pA - pA^2) %*% t(rep(1, p)) * X
      h_beta_1[,"A"] <- pA - pA^2
      h_beta_0[,"A"] <- 0
    }
    e_beta_A <- (propensity_score - propensity_score^2) %*% t(rep(1, p)) * X
    MY <- ifelse(is.na(Y), 0, Y)
    c_1 <- colMeans(h_beta_1 - h_beta_0) %*% 
      solve(t(h_beta_A) %*% ((M/propensity_score) %*% t(rep(1, p)) * X)/n)
    c_2 <- colMeans(h_beta_1 - h_beta_0) %*% 
      solve(t(h_beta_A) %*% ((M/propensity_score) %*% t(rep(1, p)) * X)/n) %*%
      (t(e_beta_A) %*% (((M * (MY - pA)/propensity_score) %*% t(rep(1, p))) * X)/n) %*%
      solve(t(e_beta_A) %*% X/n)
    tilde_var <- var(p1 - p0 - 
                       X %*% t(c_1) * (M/propensity_score * (MY - pA)) - 
                       X %*% t(c_2) * (M - propensity_score) - est)
    E.S <- map_dbl(unique(Strata), ~mean(((2*A-1)*M*(MY-pA)/propensity_score)[Strata == .]))
    var <- tilde_var - (1-2*pi)^2/pi/(1-pi) * var(E.S)
    return(c(est = est, var = var/n, 
             CI.lower = qnorm(0.025, mean = est, sd =  sqrt(var/n)),
             CI.upper = qnorm(0.975, mean = est, sd =  sqrt(var/n))))
  }
}
