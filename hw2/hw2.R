library(truncnorm)
library(tidyverse)

# 1
m <- 10000
theta.hat <- se <- numeric(2)
g <- function(x) (x^2 / sqrt(2 * pi)) * exp(-x^2 / 2)

# f1: exp(1) truncated at 1
f1 <- function(x) exp(- (x - 1))
x1 <- rexp(m) + 1
f1g <- g(x1) / f1(x1)
theta.hat[1] <- mean(f1g)
se[1] <- sd(f1g)

# f2: N(0, 1) truncated at 1
f2 <- function(x) (1 / sqrt(2 * pi)) * exp(-0.5 * x^2) / (1 - pnorm(1))
x2 <- rtruncnorm(m, a = 1, b = Inf)
f2g <- g(x2) / f2(x2)
theta.hat[2] <- mean(f2g)
se[2] <- sd(f2g)

rbind(theta.hat, se)

# 2
m <- 10000 
x <- runif(m)
h <- function(x) exp(-x * cos(pi * x))

# control variate
g <- function(x) exp(-x)
mu <- 1 - exp(-1)
c_star <- - cov(h(x), g(x)) / var(g(x))
T1 <- h(x) + c_star * (g(x) - mu)
mean(T1)

# simple MC
T2 <- h(x)
mean(T2)

# variance reduction
(var(T2) - var(T1)) / var(T2)


# 3
n <- 10000
g <- function(x) exp(x^2)

# regular MC
regMC <- g(runif(n))
theta_regMC <- mean(regMC)
sd_regMC <- sd(regMC)

# antithetic variable
u <- runif(n / 2)
a <- 0.5 * g(u) + 0.5 * g(1 - u)
theta_a <- mean(a)
sd_a <- sd(a)


# control variate
x <- runif(n)
c <- - cov(g(x), x) / (1 / 12)
cv <- g(x) + c * (x - 0.5)
theta_cv <- mean(cv)
sd_cv <- sd(cv)

# combine antithetic and control variate
u <- runif(n / 2)
A <- g(u)
B <- g(1 - u)
V <- u - 0.5
fit <- lm(B ~ I(A - B) + V)
alpha_star <- -coef(fit)[2]
c_star2 <- -coef(fit)[3]
c(alpha_star, c_star2)     

mix <- alpha_star * A + (1 - alpha_star) * B + c_star2 * V
theta_mix <- mean(mix)
sd_mix <- sd(mix) 


name <- c("Regular MC", "Antithetic", "Control Variate", "Combine")
set_names(c(theta_regMC, theta_a, theta_cv, theta_mix), name)
set_names(c(sd_regMC, sd_a, sd_cv, sd_mix), name)


# 4
# p_tilde(x) = exp(-(0.5 - x)^2 / 2) * exp(-3 * x)
# q(x) = 3 * exp(-3 * x) ~ Exp(3)
n <- 10000
theta_hat <- numeric(n)
x <- rexp(n, 3)

for(i in 1:n) {
  xi <- x[1:i]
  w <- function(x) exp(-(0.5 - xi)^2 / 2) / 3
  w_tilde <- w(xi)
  w_adj <- w_tilde / sum(w_tilde)
  theta_hat[i] <- sum(xi * w_adj)
}
plot(theta_hat, type = "l", xlab = "sample size")


