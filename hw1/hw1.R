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
    main = "theoretical and simulated pareto(2, 4) density")
dpareto <- function(x, a, b) b * a^b / (x + a)^(b + 1)
z <- seq(0, 15, 0.01)
lines(z, dpareto(z, 2, 4))

# 4
library(MASS)
rwhis.general <- function(n, Sigma) {
    x <- mvrnorm(n, mu, Sigma)
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
    rwhis.general(n, sigma)
})

# Bartlett
set.seed(100)
system.time(for(i in 1:N){
    rwhis.bartlett(d, sigma)
})

# 5
library(survival)

# Function to simulate survival data with calibrated censoring
simulate_survival_AFT <- function(n, target_censor_rate, lambda, nu, beta) {
    # 1. Covariate: Treatment Assignment (1:1 allocated)
    X <- rbinom(n, size = 1, prob = 0.5)

    # 2. Uniform draws for inverse transform
    U_t <- runif(n)

    # 3. True event times via Weibull inversion 
    T_true <- ((-log(U_t) / lambda)^(1 / nu)) / exp(X * beta)

    # 4. Bisection search to calibrate tau for target censor rate
    U_c <- runif(n)
    low <- 0.1
    high <- as.numeric(quantile(T_true, 0.999))

    for(i in 1:1000) {
        mid_tau <- (low + high) / 2
        C_temp <- mid_tau * U_c
        current_rate <- mean(T_true > C_temp)
        
        if(abs(current_rate - target_censor_rate) < 1e-4) {
            tau <- mid_tau
            break
        } else if (current_rate > target_censor_rate) {
            low <- mid_tau
        } else {
            high <- mid_tau
        }
    }
    tau <- mid_tau

    # 5. Calculate observed variables
    C_time <- tau * U_c
    Y_obs <- pmin(T_true, C_time)
    d_event <- as.numeric(T_true <= C_time)

    data.frame(
        id  = 1:n,
        treatment = X,
        true_time = T_true,
        censor_time = C_time,
        time = Y_obs,
        status = d_event
    )
}

run_mc_power_r <- function(sample_sizes = c(50, 80, 100, 150, 200, 250, 300, 500),
                           n_sim = 1000, target_censor_rate,
                           seed = 42) {
    set.seed(seed)
    lam <- 0.04
    nu <- 1.4
    beta <- -0.7

    results <- list()

    for(N in sample_sizes) {
        p_values <- numeric(n_sim)

        for(s in 1:n_sim) {
            sim_df <- simulate_survival_AFT(n = N,
                                  target_censor_rate = target_censor_rate,
                                  lambda = lam,
                                  nu = nu,
                                  beta = beta)
            X <- sim_df$treatment
            Y_obs <- sim_df$time
            d_event <- sim_df$status

            # Cox Proportion Hazards Model fit
            fit <- coxph(Surv(Y_obs, d_event) ~ X)
            p_values[s] <- summary(fit)$coefficients["X", "Pr(>|z|)"]
        }
        
        power <- mean(p_values < 0.05)

        results[[as.character(N)]] <- data.frame(
            Sample_Size = N,
            Target_HR = round(exp(beta), 2),
            Censor_Rate = paste0(target_censor_rate * 100, "%"),
            Simulations = n_sim,
            Empirical_Power = round(power, 3)
        )
    }

    do.call(rbind, results)
}
# Run R simulation
# 30% censoring rate
power_table_30 <- run_mc_power_r(target_censor_rate = 0.30)
print(power_table_30)

# 50% censoring rate
power_table_50 <- run_mc_power_r(target_censor_rate = 0.50)
print(power_table_50)

# 70% censoring rate
power_table_70 <- run_mc_power_r(target_censor_rate = 0.70)
print(power_table_70)

plot(x = power_table_30$Sample_Size, y = power_table_30$Empirical_Power,
     type = "l", col = "darkgreen", ylim = c(0, 1),
     xlab = "Sample Size", ylab = "Empirical Power",
     main = "Comparative Power Curve")
lines(x = power_table_50$Sample_Size, y = power_table_50$Empirical_Power, col = "blue")
lines(x = power_table_70$Sample_Size, y = power_table_70$Empirical_Power, col = "red")

legend("bottomright", legend = c("censoring 30%", "censoring 50%", "censoring 70%"),
       col = c("darkgreen", "blue", "red"), lty = 1, lwd = 1)
