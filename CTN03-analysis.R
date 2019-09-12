library(tidyverse)
# import data and imputation for covariates
CTN03 <- readRDS("~/Dropbox/research/clinical trial/covariate-adaptive/Covariate-adaptive/CTN03.rds")
CTN03$baseline.is.negative[which(is.na(CTN03$baseline.is.negative))] <- median(CTN03$baseline.is.negative, na.rm = T)
CTN03$COWS[which(is.na(CTN03$COWS))] <- median(CTN03$COWS, na.rm = T)
CTN03$ARSW[which(is.na(CTN03$ARSW))] <- median(CTN03$ARSW, na.rm = T)
CTN03$VAS[which(is.na(CTN03$VAS))] <- median(CTN03$VAS, na.rm = T)
CTN03 <- CTN03[!is.na(CTN03$strata),]
CTN03$strata <- as.factor(CTN03$strata)
pi = 0.5

# lab outcome -------------------------------
# specify outcome, treatment and covariates
Y <- CTN03$outcome.is.negative
A <- as.numeric(CTN03$arm == "28-day taper")
W <- select(CTN03, sex, baseline.is.negative, COWS, ARSW, VAS, strata)

# n and % missing
n <- nrow(CTN03)
CTN03_complete <- CTN03[complete.cases(CTN03),]
Y_complete <- Y[complete.cases(CTN03)]
A_complete <- A[complete.cases(CTN03)]
W_complete <- W[complete.cases(CTN03),]
n_complete <- nrow(CTN03_complete)
missing_proportion <- 1 - n_complete/n

# number of strata
n_strata <- length(unique(CTN03$strata))

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

# doubly-robust estimator
M <- !is.na(Y)
propensity.fit <- glm(M ~ ., data = data.frame(M, A, W), family = "binomial")
propensity_score <- predict(propensity.fit, type = "response")
d <- data.frame(Y = Y_complete, A = A_complete, W_complete)
glm.fit <- glm(Y~ ., data = d, family = "binomial",, weights = 1/propensity_score[M==1])
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


# retention outcome -------------------------------
# specify outcome, treatment and covariates
Y <- CTN03$complete
A <- as.numeric(CTN03$arm == "28-day taper")
W <- select(CTN03, sex, baseline.is.negative, COWS, ARSW, VAS, strata)

# n and % missing
n <- nrow(CTN03)
CTN03_complete <- CTN03[complete.cases(CTN03),]
Y_complete <- Y[complete.cases(CTN03)]
A_complete <- A[complete.cases(CTN03)]
W_complete <- W[complete.cases(CTN03),]
n_complete <- nrow(CTN03_complete)
missing_proportion <- 1 - n_complete/n

# number of strata
n_strata <- length(unique(CTN03$strata))

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
