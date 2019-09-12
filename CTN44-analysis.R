library(tidyverse)
# import data and covariate imputation
CTN44 <- readRDS("~/Dropbox/research/clinical trial/covariate-adaptive/Covariate-adaptive/CTN44.rds")
CTN44$`0`[which(is.na(CTN44$`0`))] <- TRUE
pi = 0.5

# retention outcome -------------------------------
# specify outcome, treatment and covariates
Y <- CTN44$complete
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
W <- select(CTN44, age, gender, strata, `0`)

# n and % missing
n <- nrow(cbind(Y, A, W))
CTN44_complete <- CTN44[complete.cases(cbind(Y, A, W)),]
Y_complete <- Y[complete.cases(cbind(Y, A, W))]
A_complete <- A[complete.cases(cbind(Y, A, W))]
W_complete <- W[complete.cases(cbind(Y, A, W)),]
n_complete <- nrow(CTN44_complete)
missing_proportion <- 1 - n_complete/n

# number of strata
n_strata <- length(unique(CTN44$strata))

# Unadjusted estimator
unadj <- mean(Y_complete[A_complete==1]) - mean(Y_complete[A_complete==0])
var_unadj <- var(Y_complete[A_complete==1])/pi + var(Y_complete[A_complete==0])/(1-pi)
CI_unadj <- qnorm(c(0.025, 0.975), mean = unadj, sd =  sqrt(var_unadj))

# adjusted estimator
d <- data.frame(Y = Y_complete, A = A_complete, W_complete)
d1 <- data.frame(Y = Y_complete, A = 1, W_complete)
d0 <- data.frame(Y = Y_complete, A = 0, W_complete)
glm.fit <- glm(Y~ ., data = d, family = "binomial")
p1 <- predict(glm.fit, d1, type = "response")
p0 <- predict(glm.fit, d0, type = "response")
adj <- mean(p1) - mean(p0)
Y.1 <- Y_complete[A_complete==1]
Y.0 <- Y_complete[A_complete==0]
r.1 <- Y.1 - p0[A_complete==1] * pi - p1[A_complete==1] * (1 - pi)
r.0 <- Y.0 - p0[A_complete==0] * pi - p1[A_complete==0] * (1 - pi)
var_adj <- var(r.1)/pi + var(r.0)/(1-pi)
CI_adj <- qnorm(c(0.025, 0.975), mean = adj, sd =  sqrt(var_adj))

rbind(cbind(unadj, t(CI_unadj)),
      cbind(adj, t(CI_adj))) %>% round(2)



# average lab results -------------------
# specify outcome, treatment and covariates
Y <- apply(CTN44[,7:30], 1, function(x){mean(x, na.rm = T)})
Y <- ifelse(apply(CTN44[,7:30], 1, function(x){sum(is.na(x)) > 11}), NA, Y)
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
W <- select(CTN44, age, gender, strata, `0`)

# n and % missing
CTN44 <- cbind(Y, A, W)
n <- nrow(CTN44)
CTN44_complete <- CTN44[complete.cases(CTN44),]
Y_complete <- Y[complete.cases(CTN44)]
A_complete <- A[complete.cases(CTN44)]
W_complete <- W[complete.cases(CTN44),]
n_complete <- nrow(CTN44_complete)
missing_proportion <- 1 - n_complete/n

# number of strata
n_strata <- length(unique(CTN44$strata))

# Unadjusted estimator
unadj <- mean(Y_complete[A_complete==1]) - mean(Y_complete[A_complete==0])
var_unadj <- var(Y_complete[A_complete==1])/pi + var(Y_complete[A_complete==0])/(1-pi)
CI_unadj <- qnorm(c(0.025, 0.975), mean = unadj, sd =  sqrt(var_unadj))

# adjusted estimator
d <- data.frame(Y = Y_complete, A = A_complete, W_complete)
d1 <- data.frame(Y = Y_complete, A = 1, W_complete)
d0 <- data.frame(Y = Y_complete, A = 0, W_complete)
glm.fit <- glm(Y~ ., data = d, family = "gaussian")
p1 <- predict(glm.fit, d1, type = "response")
p0 <- predict(glm.fit, d0, type = "response")
adj <- mean(p1) - mean(p0)
Y.1 <- Y_complete[A_complete==1]
Y.0 <- Y_complete[A_complete==0]
r.1 <- Y.1 - p0[A_complete==1] * pi - p1[A_complete==1] * (1 - pi)
r.0 <- Y.0 - p0[A_complete==0] * pi - p1[A_complete==0] * (1 - pi)
var_adj <- var(r.1)/pi + var(r.0)/(1-pi)
CI_adj <- qnorm(c(0.025, 0.975), mean = adj, sd =  sqrt(var_adj))

# doubly-robust estimator
M <- !is.na(Y)
propensity.fit <- glm(M ~ ., data = data.frame(M, A, W), family = "gaussian")
propensity_score <- predict(propensity.fit, type = "response")
d <- data.frame(Y = Y_complete, A = A_complete, W_complete)
glm.fit <- glm(Y~ ., data = d, family = "gaussian", weights = 1/propensity_score[M==1])
p1 <- predict(glm.fit, data.frame(A = 1, W), type = "response")
p0 <- predict(glm.fit, data.frame(A = 0, W), type = "response")
e1 <- predict(propensity.fit, data.frame(A = 1, W), type = "response")
e0 <- predict(propensity.fit, data.frame(A = 0, W), type = "response")
dr <- mean(p1) - mean(p0)
Y.1 <- ifelse(is.na(Y[A==1]), 0, Y[A==1])
Y.0 <- ifelse(is.na(Y[A==0]), 0, Y[A==0])
M.1 <- M[A==1]
M.0 <- M[A==0]
r.1 <- M.1 * (Y.1 - p1[A==1]) / e1[A==1] + pi * p1[A==1] - pi * p0[A==1]
r.0 <- M.0 * (Y.0 - p0[A==0]) / e0[A==0] - (1-pi) * p1[A==0] + (1-pi) * p0[A==0]
var_dr <- var(r.1)/pi + var(r.0)/(1-pi)
CI_dr <- qnorm(c(0.025, 0.975), mean = dr, sd =  sqrt(var_dr))

rbind(cbind(unadj, t(CI_unadj)),
      cbind(adj, t(CI_adj)),
      cbind(dr, t(CI_dr))) %>% round(2)
