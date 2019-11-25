library(tidyverse)
# import data and covariate imputation
setwd("~/Dropbox/research/clinical trial/covariate-adaptive/Covariate-adaptive/")
source("ICAD-time-to-event.R")
source("ICAD.R")
CTN30 <- readRDS("CTN30.rds")
CTN30$age[which(is.na(CTN30$age))] <- median(CTN30$age, na.rm = T)
CTN30$BL[which(is.na(CTN30$BL))] <- FALSE
CTN30 <- mutate(CTN30, strata = interaction(CTN30$heroin_history, CTN30$chronic_pain))
pi <- 0.5

# average lab results -------------------
# handling missing data: regard all missing as positive lab results
CTN30_1 <- CTN30
CTN30_1[, 8:12][is.na(CTN30[,8:12])] <- 1
Y <- apply(CTN30_1[,8:12], 1, mean)
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
n <- nrow(CTN30)
missing_proportion <- 1 - sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))

ICAD(Y,A, Strata, W, pi = pi, family = "gaussian") %>% round(3)

# handling missing data: missing if one missed more than 2 weeks
Y <- apply(CTN30[,8:12], 1, function(x){mean(x, na.rm = T)})
missing <- (is.na(CTN30[,8]) & is.na(CTN30[,9])) + is.na(CTN30[,10]) + is.na(CTN30[,11]) + is.na(CTN30[,12])
Y <- ifelse(missing > 2, NA, Y)
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
n <- nrow(CTN30)
missing_proportion <- sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))

ICAD(Y,A, Strata, W, pi = pi, family = "gaussian") %>% round(2)

# retention outcome -------------------------------
Y <- CTN30$complete1
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
n <- nrow(CTN30)
missing_proportion <- sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))

# inference
ICAD(Y,A, Strata, W, pi = pi, family = "binomial") %>% round(2)




# survival outcome: time-to-first-negative-lab-result ------
# Y: time to first negative lab result
# M: time to first missing visit
Y <- apply(CTN30[,8:12], 1, function(x){which(x == 0)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN30[,8:12], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 5)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 5)
survplot(KM, conf = "none", xlim = c(0,7), time.inc = 1)


# survival analysis-------
# Y: time to first two consecutive negative lab result
# M: time to first two consecutive missing visits
event_table <- matrix(NA, nrow = nrow(CTN30), ncol = 4)
censor_table <- matrix(NA, nrow = nrow(CTN30), ncol = 4)
for(j in 1:4){
  temp_event <- CTN30[,7+j] + CTN30[8+j]
  censor_event <- is.na(CTN30[,7+j]) + is.na(CTN30[,8+j])
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
  censor_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 2, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(censor_table, 1, function(x){which(x == 1)[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 4)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 4)

# survival outcome: time-to-first-negative-lab-result ------
# Y: time to first two consecutive negative lab result
# M: time to first missing visit
event_table <- matrix(NA, nrow = nrow(CTN30), ncol = 4)
for(j in 1:4){
  temp_event <- CTN30[,7+j] + CTN30[8+j]
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN30[,8:12], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN30$arm1 == "Enhanced")
Strata <- as.factor(CTN30$strata)
W <- select(CTN30, sex, age, BL)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 5)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 5)