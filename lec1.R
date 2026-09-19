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
