library(tidyverse)
# import data and covariate imputation
setwd("~/Dropbox/research/clinical trial/covariate-adaptive/Covariate-adaptive/")
source("ICAD-time-to-event.R")
source("ICAD.R")
CTN44 <- readRDS("CTN44.rds")
CTN44$`0`[which(is.na(CTN44$`0`))] <- TRUE
pi <- 0.5

# average lab results -------------------
# handling missing data: regard all missing as positive lab results
CTN44_1 <- CTN44
CTN44_1[, 7:30][is.na(CTN44[,7:30])] <- 1
Y <- apply(CTN44_1[,7:30], 1, mean)
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, gender, `0`)
n <- nrow(CTN44)
n_strata <- length(unique(Strata))

ICAD(Y,A, Strata, W, pi = pi, family = "gaussian") %>% round(2)

# handling missing data: missing if one missed more than 6 weeks
Y <- apply(CTN44[,7:30], 1, function(x){mean(x, na.rm = T)})
missing <- (is.na(CTN44[,seq(7, 30, by = 2)]) * is.na(CTN44[,seq(8, 30, by = 2)])) %>% apply(1,sum)
Y <- ifelse(missing > 6, NA, Y)
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, gender, `0`)

ICAD(Y,A, Strata, W, pi = pi, family = "gaussian") %>% round(2)

# retention outcome -------------------------------
Y <- CTN44$complete
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, gender, `0`)
n <- nrow(CTN44)
missing_proportion <- sum(is.na(Y))/length(Y)
n_strata <- length(unique(Strata))

ICAD(Y,A, Strata, W, pi = pi, family = "binomial") %>% round(2)


# survival outcome: time-to-first-negative-lab-result ------
# Y: time to first negative lab result
# M: time to first missing visit
Y <- apply(CTN44[,7:30], 1, function(x){which(x == 0)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN44[,7:30], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, `0`)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 6)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 6)

# survival analysis-------
# Y: time to first two consecutive negative lab result
# M: time to first two consecutive missing visits
event_table <- matrix(NA, nrow = nrow(CTN44), ncol = 23)
censor_table <- matrix(NA, nrow = nrow(CTN44), ncol = 23)
for(j in 1:23){
  temp_event <- CTN44[,6+j] + CTN44[7+j]
  censor_event <- is.na(CTN44[,6+j]) + is.na(CTN44[,7+j])
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
  censor_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 2, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(censor_table, 1, function(x){which(x == 1)[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, `0`)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 23)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 23)

# survival analysis-------
# Y: time to first two consecutive negative lab result
# M: time to first missing visits
event_table <- matrix(NA, nrow = nrow(CTN44), ncol = 23)
for(j in 1:23){
  temp_event <- CTN44[,6+j] + CTN44[7+j]
  event_table[,j] <- ifelse(!is.na(temp_event) & temp_event == 0, 1, 0)
}
Y <- apply(event_table, 1, function(x){which(x == 1)[1]})
Y[is.na(Y)] <- 999
M <- apply(CTN44[,7:30], 1, function(x){which(is.na(x))[1]})
M[is.na(M)] <- 999
E <- pmin(Y, M)
C <- Y <= M
A <- as.numeric(CTN44$arm == "Therapeutic Education System (TES)")
Strata <- as.factor(CTN44$strata)
W <- select(CTN44, age, `0`)
ICAD_tte(E, C, A, Strata, W, pi = 0.5, tau = 10)
ICAD_tte(E, C, A, Strata, W = NULL, pi = 0.5, tau = 10)
