#--------------------------------------------------------
# The Inverse Transform Method 
#--------------------------------------------------------

# Ex1: Sampling from finite populations
sample(0:1, size = 10, replace = TRUE)
sample(1:100, size = 6, replace = FALSE)

# permutation of letters a-z sample
sample(letters) 

# sample from a multinomial distribution
x <- sample(1:3, size = 10, replace = TRUE, prob = c(0.2, 0.5, 0.3))
table(x)

# Ex2: Continuous case
# Generate random sample for f_x(x) = 3x^2, 0 < x < 1
# F-1(x) = x^3, F^-1(u) = u^(1/3)
n <- 1000
for (i in 1:n){
    x[i] <- runif(1)^(1/3)
}
print(x)
hist(x, prob = TRUE, main = bquote(f(x)==3*x^2), xlab = "x", ylab = "Density")

# Compare with the theoretical density function
y <- seq(0, 1, .01)
lines(y, 3*y^2, col = "red", lwd = 2)

# Ex3: Exponential distribution
lambda <- 0.5
x <- -1/lambda * log(runif(1000))
x

# Ex4: Geometric distribution
u <- runif(1000)
p <- 0.25
q <- 1 - p
x1 <- ceiling(log10(1 - u) / log10(q)) - 1                                                                                                    
x1
x2 <- floor(log(1 - u) / log(q))
x2

# Ex5: Poisson distribution
simPoisson <- function(n.gen, lambda) {
    u <- runif(n.gen)
    sim.vector <- rep(0, n.gen)
    for(j in 1:n.gen){
        i <- 0
        p <- exp(-lambda)
        F <- p
        while(u[j] > F){
            p <- lambda * p / (i + 1)
            F <- F + p
            i <- i + 1
        }
        sim.vector[j] <- i
    }
    sim.vector
}
out1 <- simPoisson(1000, 5)
out1 
mean(out1)
xtabs(~ out1) / 1000
barplot(xtabs(~ out1) / 1000)

# Ex6: Logarithmic distribution
rlogarithmic <- function(n, theta) {
    # returns a random logarithmic(theta) sample size n
    u <- runif(n)

    # set the initial length of cdf vector
    N <- ceiling(-16 / log10(theta))
    k <- 1:N
    a <- -1 / log(1 - theta)
    fk <- exp(log(a) + k * log(theta) - log(k))
    Fk <- cumsum(fk)
    
    x <- integer(n)
    for(i in 1:n) {
        x[i] <- as.integer(sum(u[i] > Fk))
        while(x[i] == N) {
            logf <- log(a) + (N + 1) * log(theta) - log(N + 1)
            fk <- c(fk, exp(logf))
            Fk <- c(Fk, Fk[N] + fk[N + 1])
            N <- N + 1
            x[i] <- as.integer(sum(u[i] > Fk))
        }
    }
    x + 1
}

n <- 1000
theta <- 0.5
x <- rlogarithmic(n, theta)

# Compute density of logarithmic(theta) for comparison
k <- sort(unique(x))
p <- -1 / log(1 - theta) * theta^k / k
se <- sqrt(p * (1 - p) / n)
round(rbind(table(x) / n, p, se), 3)

#--------------------------------------------------------
# The Acceptance-Rejection Method
#--------------------------------------------------------

# Ex1: Triangular distribution
rejectionK <- function(fx, a, b, K) {
    while(TRUE) {
        x <- runif(1, a, b) # 從candidate g抽一個候選點y
        y <- runif(1, 0, K) # 抽一個高度cu or Ku
        if(y < fx(x)) return(x) 
    }
}

fx <- function(x) {
    if((0 < x) && (x < 1)) {
        return(x)
    } else if((1 < x) && (x < 2)) {
        return(2-x)
    } else {
        return(0)
    }
}

# generate a sample of size 3000
set.seed(21)
n <- 300000000
observation <- rep(0, n)
for(i in 1:n) {
    observation[i] <- rejectionK(fx, 0, 2, 1)
}
observation
hist(observation, breaks = seq(0, 2, by = 0.1), freq = FALSE,  
    ylim = c(0, 1.05), main = "")

# Ex2: Gamma distribution
gamma.sim <- function(lambda, m) {
    # sim a gamma(lambda, m) rv using rejection with an exp envelope
    # assumes m > 1 and lambda > 0
    f <- function(x) lambda^m * x^(m - 1) * exp(-lambda * x) / gamma(m)
    h <- function(x) lambda / m * exp(-lambda / m * x)
    k <- m^m * exp(1 - m) / gamma(m)
    
    while(TRUE) {
        X <- -log(runif(1)) * m / lambda
        Y <- runif(1, 0, k * h(X))
        if(Y < f(X)) return(X)
    }
}

set.seed(1999)
n <- 10000
g <- rep(0, n)
for(i in 1:n){
    g[i] <- gamma.sim(1, 2)
}
hist(g, breaks = 20, freq = F, xlab = "x", ylab = "pdf f(x)", 
    main = "theoretical and simulated gamma(1, 2) density")
x <- seq(0, max(g), 0.1)
lines(x, dgamma(x, 2, 1))

# Ex3: Poisson distribution
poisson.sim <- function(c) {
    while(TRUE) {
        j <- floor(log(runif(1)) / log(0.75))
        u <- runif(1)
        p <- exp(-3) * 3^j / factorial(j)
        q <- 0.25 * (1 - 0.25)^j
        if(u < p / (c * q)) return(j)
    }

}
n <- 3000
observation <- rep(0, n)
for(i in 1:n) {
    observation[i] <- poisson.sim(c = 2.12)
} 
hist(observation, breaks = 20, freq = F, xlab = "x", ylab = "pmf P(X = x)", 
    main = "simulated poisson density")

#--------------------------------------------------------
# Transformation Method
#--------------------------------------------------------
# Ex1: Beta distribution Beta(3, 2)
n <- 1000
a <- 3
b <- 2
u <- rgamma(n, shape = a, rate = 1)
v <- rgamma(n, shape = b, rate = 1)
x <- u / (u + v)

q <- qbeta(ppoints(n), a, b)
qqplot(q, x, cex = 0.25, xlab = "Beta(3, 2)", ylab = "Sample")
abline(0, 1)

# Ex2: Logarithmic distribution
rlogarithmic <- function(n, theta) {
    stopifnot(all(theta > 0 & theta < 1))
    th <- rep(theta, length = n)
    u <- runif(n)
    v <- runif(n)
    x <- floor(1 + log(v) / log(1 - (1 - th)^u))
    return(x)
}

n <- 1000
theta <- 0.5
u <- runif(n)
v <- runif(n)
x <- floor(1 + log(v)) / log(1 - (1 - theta)^u)
k <- 1:max(x)
p <- -1 / log(1 - theta) * theta^k / k 
se <- sqrt(p * (1 - p) / n)
p.hat <- tabulate(x) / n
print(round(rbind(p.hat, p, se), 3))

#--------------------------------------------------------
# Convolution 
#--------------------------------------------------------
# Chi-square distribution
# Z_1,...,Z_v are iid N(0, 1) r.v.
# V = Z_1 + ... + Z_v ~ chi_square(v)

n <- 1000
nu <- 2
X <- matrix(rnorm(n * nu), n, nu)^2 # matrix of sq. normals
# sum the squared normals across each row
# method 1
y <- rowSums(X)
# method 2
y <- apply(X, MARGIN = 1, FUN = sum)

mean(y)
mean(y^2)

# Convolution of two normal -> use mixture to handle
# X_1 ~ N(0, 1), X_2 ~ N(3, 1) and independent
# S = X_1 + X_2 is the convolution of X_1 and X_2 
normalConvolution.sim <- function() {
    w <- sample(c(1L, 2L), 1)
    if(w == 1) {
        x <- rnorm(1, mean = 0, sd = 1)
    } else {
        x <- rnorm(1, mean = 3, sd = 1)
    }
    return(x)
}
n <- 1000
s <- rep(0, n)
for(i in 1:n) {
    s[i] <- normalConvolution.sim()
}
s

# X_1 ~ Gamma(2, 2), X_2 ~ Gamma(2, 4)
# Generate the sample from the convolution S = X_1 + X_2 and the mixture 
n <- 1000
x1 <- rgamma(n, 2, 2)
x2 <- rgamma(n, 2, 4)
s <- x1 + x2
u <- runif(n)
k <- as.integer(u > 0.5)
x <- k * x1 + (1 - k) * x2

par(mfcol = c(1, 2))
hist(s, prob = TRUE)
hist(x, prob = TRUE)
par(mfcol = c(1, 1))

# Mixture of several gamma distribution
# my first try
k <- c(1L, 2L, 3L, 4L, 5L)
p <- (1:5) / 15
mixGamma.sim <- function(r) {
    x <- sample(k, 1, prob = p)
    if(x == 1) {
        return(rgamma(1, shape = r, rate = 1/x))
    } else if (x == 2) {
       return(rgamma(1, shape = r, rate = 1/x))
    } else if (x == 3) {
       return(rgamma(1, shape = r, rate = 1/x))
    } else if (x == 4) {
       return(rgamma(1, shape = r, rate = 1/x))
    } else {
       return(rgamma(1, shape = r, rate = 1/x))
    } 
}

n <- 5000
y <- rep(0, n)
for(i in 1:n) {
    y[i] <- mixGamma.sim(3)
}
y
plot(density(y), xlim = c(0, 40), ylim = c(0, .3), 
    lwd = 3, xlab = "x", main = "")

# lecture note
n <- 5000
k <- sample(1:5, size = n, replace = TRUE, prob = (1:5)/15)
rate = 1/k
x <- rgamma(n, shape = 3, rate = rate)
plot(density(x), xlim = c(0, 40), ylim = c(0, .3), 
    lwd = 3, xlab = "x", main = "")
for(i in 1:5) {
    lines(density(rgamma(n, 3, 1/i)))
}

# ex2 for mixture of several gamma distribution
n <- 5000
lambda <- seq(1, 3, 0.5)
theta <- c(.1, .2, .2, .3, .2) 
k <- sample(lambda, size = n,replace = T, prob = theta)
x <- rgamma(n, shape = 3, rate = k)

# A function to compute the density f(x)
f <- function(x, lambda, theta) {
    # density of the mixture at the point x
    sum(dgamma(x, 3, lambda) * theta)
}

p <- c(.1, .2, .2, .3, .2) 
lambda <- seq(1, 3, 0.5)
x <- seq(0, 8, length = 200)
dim(x) <- length(x) # need for apply
# compute density of the mixture f(x) along x
y <- apply(x, 1, f, lambda = lambda, theta = p)

plot(x, y, type = "l", ylim = c(0, .85), lwd = 3, ylab = "Density")
for(j in 1:5) {
    y <- apply(x, 1, dgamma, shape = 3, rate = lambda[j])
    lines(x, y)
}

# Poisson-Gamma mixture distribution
n <- 1000
r <- 4
beta <- 3
lambda <- rgamma(n, r, beta) # lambda is random
# Now supply the sample of lambda's as the Poisson mean
x <- rpois(n, lambda)
# Compare with negative binomial
mix <- tabulate(x + 1) / n
negbin <- round(dnbinom(0:max(x), r, beta / (1 + beta)), 3)
se <- sqrt(negbin * (1 - negbin) / n)
round(rbind(mix, negbin, se), 3)

#--------------------------------------------------------
# Multivariate Normal Distribution
#--------------------------------------------------------
# Ex15: Generate a bivariate normal sample 
# Spectral decomposition method
mu <- c(0, 0)
Sigma <- matrix(c(1, .9, .9, 1), nrow = 2, ncol = 2)
rmvn.eigen <- function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    ev <- eigen(Sigma, symmetric = TRUE)
    lambda <- ev$values
    V <- ev$vectors
    R <- V %*% diag(sqrt(lambda)) %*% t(V)
    Z <- matrix(rnorm(n * d), nrow = n, ncol = d)
    X <- Z %*% R + matrix(mu, n, d, byrow = T)
    X
}

# generate the sample
X <- rmvn.eigen(1000, mu, Sigma)
plot(X, xlab = "x", ylab = "y", pch = 20)
print(colMeans(X))
print(cor(X))

# SVD method
rmvn.svd <- function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    S <- svd(Sigma)
    R <- S$u %*% diag(sqrt(S$d)) %*% t(S$v) # sq. root Sigma
    Z <- matrix(rnorm(n * d), nrow = n, ncol =d)
    X <- Z %*% R + matrix(mu, n, d, byrow = TRUE)
    X
}

# Choleski Method
rmvn.Choleski <- function(n, mu, Sigma) {
    # generate n random vectors from MVN(mu, Sigma)
    # dimension is inferred from mu and Sigma
    d <- length(mu)
    Q <- chol(Sigma) # Choleski factorization of Sigma
    Z <- matrix(rnorm(n * d), nrow = n, ncol = d)
    X <- Z %*% Q + matrix(mu, n, d, byrow = TRUE)   
    X 
}

#--------------------------------------------------------
# Perfomance comparisons
#--------------------------------------------------------
library(MASS)
library(mvtnorm)
n <- 100 # sample size
d <- 30 # dimension
N <- 2000 # iterations
mu <- numeric(d)

# Spectral decomposition
set.seed(100)
system.time(for(i in 1:N){
    rmvn.eigen(n, mu, cov(matrix(rnorm(n * d), n, d)))
})

# SVD method
set.seed(100)
system.time(for(i in 1:N){
    rmvn.svd(n, mu, cov(matrix(rnorm(n * d), n, d)))
})

# Choleski method
set.seed(100)
system.time(for(i in 1:N){
    rmvn.Choleski(n, mu, cov(matrix(rnorm(n * d), n, d)))
})

# mvrnorm
set.seed(100)
system.time(for(i in 1:N){
    mvrnorm(n, mu, cov(matrix(rnorm(n * d), n, d)))
})

# rmvnorm
set.seed(100)
system.time(for(i in 1:N){
    rmvnorm(n, mu, cov(matrix(rnorm(n * d), n, d)))
})


