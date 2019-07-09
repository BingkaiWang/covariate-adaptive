set.seed(123)
source("CADE.R")
summary_table <- matrix(NA, ncol = 2, nrow = 3)
colnames(summary_table) <- c("unadjusted", "adjusted")
rownames(summary_table) <- c("simulation", "theoretical", "ignored")
# estimate variance by simulation ---------------------
sim_size = 1000
results <- matrix(NA,  nrow = sim_size, ncol = 2)
# parameters of data generating mechanism
n <- 2000
num.strata <- 4
pi <- 0.5
for(iter in 1:sim_size){
  # generate baseline variable and potential outcomes --------
  W <- runif(n,-2,2)
  W2 <- runif(n,-2,2)
  epsi.1 <- abs(W2) * rnorm(n)
  epsi.0 <- abs(W2) * rnorm(n)
  m0 <- W2 * (abs(W2) >= 1) + W2^2 * (abs(W2) < 1) + sin(W)
  m1 <- W2^2 * (abs(W2) >= 1) + W2 * (abs(W2) < 1) + sin(W)
  Y.0 <- (m0 + epsi.0) > 0.5
  Y.1 <- (m1 + epsi.1) > 0.5
  
  # generate strata -----------
  bounds<-seq(-2,2,length.out = num.strata+1); # careful with bounds
  I.S   <- matrix(0,n,num.strata);
  for (s in 1:num.strata){
    I.S[,s]<- (W2>bounds[s])*(W2<=bounds[s+1]);
  } 
  
  # treatment assignment (stratified block randomziation) ------------
  A <- rep(NA, n)
  lambda <- 3/4
  D.s <- rep(0,num.strata)  # inbalance in each stratum s
  n.s <- rep(0,num.strata)  # sample size in each s
  for (j in 1:n){
    index   <- which(I.S[j,]==1);
    n.s[index]<-n.s[index]+1;
    phi  <- lambda*(D.s[index]<0)+pi*(D.s[index]==0)+(1-lambda)*(D.s[index]>0);
    A[j]    <-(runif(1)<=phi);
    D.s[index]<-(1/pi)*sum(A*(I.S[,index]==1), na.rm = T)-n.s[index]; 
  }
  A <- as.numeric(A)
  Y <- Y.1*A + Y.0*(1-A)
  
  # estimating ATE  ------------
  sim.data <- data.frame(Y, A, W, I.S)
  sim.data1 <- data.frame(Y, A = 1, W, I.S)
  sim.data0 <- data.frame(Y, A = 0, W, I.S)
  
  # unadjusted
  unadj <- mean(Y[A==1]) - mean(Y[A==0])
  
  # adjusting for both strata and covariates
  glm.fit3 <- glm(Y~ 0 + ., data = sim.data, family = "binomial")
  adj_both <- mean(predict(glm.fit3, sim.data1, type = "response") - predict(glm.fit3, sim.data0, type = "response"))
  
  results[iter,] <- c(unadj, adj_both)
}
summary_table[1,] <- apply(results, 2, sd) * sqrt(n)


# estimate variance by formula --------------
n <- 100000
# generate baseline variable and potential outcomes
W <- runif(n,-2,2)
W2 <- runif(n,-2,2)
epsi.1 <- abs(W2) * rnorm(n)
epsi.0 <- abs(W2) * rnorm(n)
m0 <- W2 * (abs(W2) >= 1) + W2^2 * (abs(W2) < 1) + sin(W)
m1 <- W2^2 * (abs(W2) >= 1) + W2 * (abs(W2) < 1) + sin(W)
Y.0 <- (m0 + epsi.0) > 0.5
Y.1 <- (m1 + epsi.1) > 0.5

# generate strata 
bounds<-seq(-2,2,length.out = num.strata+1); # careful with bounds
I.S   <- matrix(0,n,num.strata);
for (s in 1:num.strata){
  I.S[,s]<- (W2>bounds[s])*(W2<=bounds[s+1]);
}

# treatment assignment (stratified block randomziation)
A <- rep(NA, n)
lambda <- 3/4
D.s <- rep(0,num.strata)  # inbalance in each stratum s
n.s <- rep(0,num.strata)  # sample size in each s
for (j in 1:n){
  index   <- which(I.S[j,]==1);
  n.s[index]<-n.s[index]+1;
  phi  <- lambda*(D.s[index]<0)+pi*(D.s[index]==0)+(1-lambda)*(D.s[index]>0);
  A[j]    <-(runif(1)<=phi);
  D.s[index]<-(1/pi)*sum(A*(I.S[,index]==1), na.rm = T)-n.s[index]; 
}
A <- as.numeric(A)
Y <- Y.1*A + Y.0*(1-A)

summary_table[2,1] <- CADE(Y = Y, A = A, Strata = I.S, pi = pi, family = "binomial")[2]
summary_table[2,2] <- CADE(Y = Y, A = A, Strata = I.S, pi = pi, W = W, family = "binomial")[2]

# estimate variance if covariate-adaptive design is ignored
summary_table[3,1] <- CADE(Y = Y, A = A, pi = pi)[2]
summary_table[3,2] <- CADE(Y = Y, A = A, W = cbind(I.S[,-1], W), pi = pi, family = "binomial")[2]


round(summary_table^2, 2)
