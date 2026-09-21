# 1
# Candidate distribution: unif(0, 1)
beta.sim <- function(alpha, beta, c) {
    stopifnot(alpha >= 1 & beta >= 1)
    while(TRUE) {
        g <- runif(1, 0, 1)
        cu <- runif(1, 0, c)
        f <- g^(alpha - 1) * (1 - g)^(beta - 1) / (beta(alpha, beta))
        if(cu < f) return(g)
    }
}

# find the best c: we got the largest f(x)/g(x) at x = 2/3
# which means the best c we got is 
c <- (2/3)^2 * (1/3) / beta(3, 2)
n <- 1000
x <- rep(0, n)
for(i in 1:n) {
    x[i] <- beta.sim(3, 2, c)
}

hist(x, breaks = 100, freq = F, xlab = "x", ylab = "pdf f(x)", 
    main = "theoretical and simulated beta(3, 2) density")
y <- seq(0, 1, 0.001)
lines(y, dbeta(y, 3, 2))

# 2
# f(x) = (2 / (pi*R^2)) * sqrt(R^2 - x^2)
# g(x) = 1 / 2R
# f(x) / g(x) = (4 / (pi * R)) * sqrt(R^2 - x^2)
# f(x) / g(x) has the largest value at 4/pi when x = 0
# c = 4 / pi
semicircle.sim <- function(R) {
    while(TRUE) {
         y <- runif(1, -R, R)
         d <- (4 / (pi * R)) * sqrt(R^2 - y^2)
         u <- runif(1, 0, 4/pi)
         if(u < d) return(y)
    }
}

dsemicircle <- function(R, x) {
    (2 / (pi*R^2)) * sqrt(R^2 - x^2)
}

n <- 1000
x <- rep(0, n)
R <- 10
for(i in 1:n) {
    x[i] <- semicircle.sim(R)
}
hist(x, freq = F, xlab = "x", ylab = "pdf f(x)",
    main = "theoretical and simulated semicircle density")
y <- seq(min(x), max(x), 0.1)
lines(y, dsemicircle(R = 10, x = y))

# 3
# cdf: 1 - (a / (x + a))^b
# cdf inverse F^-1(u): a * ((1 - u)^(-1/b) - 1)
pareto.sim <- function(a, b) {
    u <- runif(1)
    x <- a * ((1 - u)^(-1/b) - 1)
    x
}
n <- 1000
y <- rep(0, n)
for(i in 1:n) {
    y[i] <- pareto.sim(2, 4)
}

hist(y, breaks = 50, freq = F, xlab = "x", ylab = "pdf f(x)", 
    main = "theoretical and simulated pareto(3, 2) density")
dpareto <- function(x, a, b) b * a^b / (x + a)^(b + 1)
z <- seq(0, 15, 0.01)
lines(z, dpareto(z, 2, 4))

# 4
library(MASS)
rwhis.general <- function(n, Sigma) {
    x <- mvrnorm(n, mu, sigma)
    M <- t(x) %*% x
    M
}
rwhis.bartlett <- function(d, sigma) {
    A <- matrix(rnorm(d * d), nrow = d, ncol = d)
    A[!lower.tri(A)] <- 0
    for(i in 1:d) {
        diag(A)[i] <- sqrt(rchisq(1, df = n - i + 1))
    }
    L <- t(chol(sigma))
    M <- L %*% A %*% t(A) %*% t(L)

}

n <- 100 # sample size
d <- 2 # dimension
N <- 2000 # iterations
sigma <- matrix(c(1, 0.5, 0.5, 2), 2,2, byrow = T)
mu <- numeric(d)

# general
set.seed(100)
system.time(for(i in 1:N){
    rwhis.general(n, Sigma)
})

# Bartlett
set.seed(100)
system.time(for(i in 1:N){
    rwhis.bartlett(d, Sigma)
})
