## Loading the data (replace this with online download)
elephant_data = read.csv(file="./data/elephant_data.csv")

# Defining time of day variable
elephant_data$timestamp = strptime(elephant_data$timestamp, 
                                   "%Y-%m-%d %H:%M:%S", tz = "GMT")
hours = as.numeric(format(elephant_data$timestamp, "%H"))
minutes = as.numeric(format(elephant_data$timestamp, "%M"))

elephant_data$timeOfDay = hours + (minutes/60)

hmm_data <- moveHMM::prepData(trackData = elephant_data, 
                              coordNames = c("location.long", "location.lat"), type = "LL")

## mllk

hmm_data$timeOfDay2 = round(hmm_data$timeOfDay)


mllk = function(theta.star, X, N=3, K=1, L=24){
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  mu.turn = theta.star[2*N+1:N]
  kappa = exp(theta.star[3*N+1:N])
  beta = matrix(theta.star[4*N+1:(N*(N-1)*(1+2*K))], ncol = 1+2*K)
  Gamma = Lcpp::tpm_p(seq(0.7, 23.7, length = 24), L, beta, degree=K)
  delta = Lcpp::stationary_p(Gamma, t = X$timeOfDay2[1])
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  -Lcpp::forward_p(delta, Gamma, allprobs, X$timeOfDay2)
}

theta.star = c(log(c(0.1, 0.4, 1, 0.1, 0.3, 0.6)), rep(0,3), log(c(0.1, 0.9, 2)), rep(-2,6), rep(0, 12))
N=3;K=1;L=24

hmm_data$step[which(hmm_data$step==0)] = 1e-10

mllk(theta.star, hmm_data)

mod = nlm(mllk, theta.star, X=hmm_data, iterlim = 1000, print.level = 2)
theta.star = mod$estimate
mu = exp(theta.star[1:N])
sigma = exp(theta.star[N+1:N])
mu.turn = theta.star[2*N+1:N]
kappa = exp(theta.star[3*N+1:N])
beta = matrix(theta.star[4*N+1:(N*(N-1)*(1+2*K))], ncol = 1+2*K)
rownames(beta) = c("2->1", "3->1", "1->2", "3->2", "1->3", "2->3")
colnames(beta) = c("intercept", "sin1", "cos1")

beta_moveHMM = t(mod_tod$mle$beta)

Gamma = Lcpp::tpm_p(1:24, 24, beta, degree=1)
Gamma_moveHMM = Lcpp::tpm_p(1:24, 24, beta_moveHMM, degree=1, byrow = TRUE)
Delta = Lcpp::stationary_p(Gamma)

plot(Gamma[2,3,], ylim = c(0,1), type = "l")
lines(Gamma_moveHMM[2,3,])

par(mfrow = c(1,2))
plot(Delta[,1], type = "l", lwd = 2, col = "orange")
lines(Delta[,2], lwd = 2, col = "deepskyblue")
lines(Delta[,3], lwd = 2, col = "seagreen2")

Delta_moveHMM = Lcpp::stationary_p(Gamma_moveHMM)
plot(Delta_moveHMM[,1], type = "l", lwd = 2, col = "orange")
lines(Delta_moveHMM[,2], lwd = 2, col = "deepskyblue")
lines(Delta_moveHMM[,3], lwd = 2, col = "seagreen2")

Rho = matrix(NA, 24, 3)
for(t in 1:24){
  Rho[t,] = Lcpp::stationary(Gamma_moveHMM[,,t])
}

plot(Rho[,1], type = "l", lwd = 2, col = "orange")
lines(Rho[,2], lwd = 2, col = "deepskyblue")
lines(Rho[,3], lwd = 2, col = "seagreen2")

par(mfrow = c(3,3))
for(i in 1:3){
  for(j in 1:3){
    plot(Gamma[i,j,], type = "l", ylim = c(0,1))
  }
}

Gamma2 = tpm_p(beta=beta_m2, degree = 2, byrow = TRUE)
Delta2 = stationary_p(Gamma2)

plot(Delta2[,1], type = "l", lwd = 2, col = "orange")
lines(Delta2[,2], lwd = 2, col = "deepskyblue")
lines(Delta2[,3], lwd = 2, col = "seagreen2")


# 2 state model

theta.star = c(log(c(0.1, 1, 0.1, 0.6)), rep(0,2), log(c(0.2, 2)), rep(-2,2), rep(0, 4))
mod = nlm(mllk, theta.star, X=hmm_data, N=2, iterlim = 1000, print.level = 2)
N=2
theta.star = mod$estimate
mu = exp(theta.star[1:N])
sigma = exp(theta.star[N+1:N])
mu.turn = theta.star[2*N+1:N]
kappa = exp(theta.star[3*N+1:N])
beta = matrix(theta.star[4*N+1:(N*(N-1)*(1+2*K))], ncol = 1+2*K)
Gamma = Lcpp::tpm_p(1:L, L, beta, degree=K)
Delta = Lcpp::stationary_p(Gamma)
plot(Delta[,1], ylim = c(0,1))

plot(Gamma[1,2,], ylim = c(0,0.5))
plot(Gamma[2,1,], ylim = c(0,0.5))

theta.star = c(log(c(0.1, 1, 0.1, 0.6)), rep(0,2), log(c(0.2, 2)), rep(-2,2), rep(0, 8))
mod2 = nlm(mllk, theta.star, X=hmm_data, N=2, K=2, iterlim = 1000, print.level = 2)
N=2; K=2
theta.star = mod2$estimate
mu = exp(theta.star[1:N])
sigma = exp(theta.star[N+1:N])
mu.turn = theta.star[2*N+1:N]
kappa = exp(theta.star[3*N+1:N])
beta = matrix(theta.star[4*N+1:(N*(N-1)*(1+2*K))], ncol = 1+2*K)
Gamma = Lcpp::tpm_p(1:L, L, beta, degree=K)
Delta = Lcpp::stationary_p(Gamma)
plot(Delta[,1], ylim = c(0,1))

plot(Gamma[1,2,], ylim = c(0,0.5))
plot(Gamma[2,1,], ylim = c(0,0.5))

Gamma2 = Lcpp::tpm_p(seq(0,24,length=200), 24, beta, degree=K)
plot(seq(0,24,length=200), Gamma2[1,2,], type = "l", ylim = c(0,0.5))
plot(seq(0,24,length=200), Gamma2[2,1,], type = "l", ylim = c(0,0.5))
