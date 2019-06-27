set.seed(123)
source("CADE.R")
summary_table <- matrix(NA, ncol = 3, nrow = 3)
colnames(summary_table) <- c("unadj", "strata", "both")
rownames(summary_table) <- c("simulation", "theoretical", "ignored")
# estimate variance by simulation ---------------------
sim_size = 1000
results <- matrix(NA,  nrow = sim_size, ncol = 3)
# parameters of data generating mechanism
n <- 200
num.strata <- 4
pi <- 0.5
for(iter in 1:sim_size){
  # generate baseline variable and potential outcomes --------
  W <- runif(n,-2,2)
  W2 <- runif(n,-2,2)
  epsi.1 <- abs(W) * 2 * rt(n,3)/3
  epsi.0 <- abs(W) * rt(n,3)/3
  m0 <- W + abs(W2)
  m1 <- 2 * W + W2^2
  Y.0 <- (m0 + epsi.0) > 0.5
  Y.1 <- (m1 + epsi.1) > 0.5
  
  # generate strata -----------
  bounds<-seq(-2,2,length.out = num.strata+1); # careful with bounds
  I.S   <- matrix(0,n,num.strata);
  for (s in 1:num.strata){
    I.S[,s]<- (W2>bounds[s])*(W2<=bounds[s+1]);
  } 
  
  # treatment assignment (stratified block randomziation) ------------
  A <- rbinom(n, 1, pi)
  Y <- Y.1*A + Y.0*(1-A)
  
  # estimating ATE  ------------
  sim.data <- data.frame(Y, A, W, I.S)
  
  # unadjusted
  unadj <- mean(Y[A==1]) - mean(Y[A==0])
  
  # adjusting for strata only
  glm.fit1 <- glm(Y~ . + 0 - W, data = sim.data, family = "binomial")
  adj_strata <- glm.fit1$coefficients["A"]
  
  # adjusting for both strata and covariates
  glm.fit3 <- glm(Y~ 0 + ., data = sim.data, family = "binomial")
  adj_both <- glm.fit3$coefficients["A"]
  
  results[iter,] <- c(unadj, adj_strata, adj_both)
}
summary_table[1,] <- apply(results, 2, sd) * sqrt(n)


# estimate variance by formula --------------
n <- 2000
num.strata <- 4
pi <- 0.5
# generate baseline variable and potential outcomes
W <- runif(n,-2,2)
W2 <- runif(n,-2,2)
epsi.1 <- abs(W) * 2 * rt(n,3)/3
epsi.0 <- abs(W) * rt(n,3)/3
m0 <- W + abs(W2)
m1 <- 2 * W + W2^2
Y.0 <- (m0 + epsi.0) > 0.5
Y.1 <- (m1 + epsi.1) > 0.5

# generate strata 
bounds<-seq(-2,2,length.out = num.strata+1); # careful with bounds
I.S   <- matrix(0,n,num.strata);
for (s in 1:num.strata){
  I.S[,s]<- (W2>bounds[s])*(W2<=bounds[s+1]);
} 

# treatment assignment (stratified block randomziation)
A <- rbinom(n, 1, pi)
Y <- Y.1*A + Y.0*(1-A)

# variance calculating
n.1 <- sum(A) 
n.0 <- n-n.1
Y.0 <- Y*(1-A)
Y.1 <- Y*A
naI.S       <- I.S # here we replace zeros with NA
naI.S[naI.S==0]<-NA
n.0.s<-colSums(naI.S*(1-A),na.rm=TRUE)
n.1.s<-colSums(naI.S*A,na.rm=TRUE)
sum.0.s<-colSums(Y.0*naI.S,na.rm=TRUE)
sum.1.s<-colSums(Y.1*naI.S,na.rm=TRUE)
mu.0.s <- sum.0.s/n.0.s
mu.1.s <- sum.1.s/n.1.s
Y1.bar   <-sum(Y.1)/n.1 
Y0.bar   <-sum(Y.0)/n.0
# varsigma tilde Y
p.s <- colMeans(I.S)
term1<-(1/n.1)*sum(Y^2*A)-sum(p.s*mu.1.s^2,na.rm=TRUE)
term0<-(1/n.0)*sum(Y^2*(1-A))-sum(p.s*mu.0.s^2,na.rm=TRUE)
v.sg.tilde.Y <- term1/pi + term0/(1-pi)

# varsigma H  
t1    <- sum(p.s*(mu.1.s-Y1.bar)^2,na.rm=TRUE)
t2    <- sum(p.s*(mu.0.s-Y0.bar)^2,na.rm=TRUE)
t3    <- 2*sum(p.s*(mu.1.s-Y1.bar)*(mu.0.s-Y0.bar),na.rm=TRUE)
v.sg.H<- t1 + t2 - t3

summary_table[2,1] <- sqrt(v.sg.tilde.Y + v.sg.H)
summary_table[2,2] <- CADE(Y = Y, A = A, Strata = I.S, tau = pi * (1-pi), family = "binomial")[2]
summary_table[2,3] <- CADE(Y = Y, A = A, Strata = I.S, W = W, tau = pi * (1-pi), family = "binomial")[2]

# estimate variance if covariate-adaptive design is ignored
summary_table[3,1] <- CADE(Y = Y, A = A)[2]
summary_table[3,2] <- CADE(Y = Y, A = A, W = I.S[,-1], family = "binomial")[2]
summary_table[3,3] <- CADE(Y = Y, A = A, W = cbind(I.S[,-1], W), family = "binomial")[2]


round(summary_table,2)
