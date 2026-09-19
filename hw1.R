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
