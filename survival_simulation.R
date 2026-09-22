#----------------------------------------------------------------
# Generate survival data (Y_i, d_i)
#----------------------------------------------------------------
library(survival)

# Function to simulate survival data with calibrated censoring
simulate_survival <- function(n = 1000,
                              target_censor_rate = 0.30,
                              lambda = 0.04,
                              nu = 1.4,
                              beta = -0.7,
                              seed = 42) {
    set.seed(seed)
    
    # 1. Covariate: Treatment Assignment (1:1 allocated)
    X <- rbinom(n, size = 1, prob = 0.5)

    # 2. Uniform draws for inverse transform
    U_t <- runif(n)

    # 3. True event times via Weibull inversion
    T_true <- (-log(U_t) / (lambda * exp(X * beta)))^(1 / nu)

    # 4. Bisection search to calibrate tau for target censor rate
    U_c <- runif(n)
    low <- 0.1
    high <- as.numeric(quantile(T_true, 0.999))

    for(i in 1:100) {
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

# Generate sample dataset
sim_df <- simulate_survival(n = 1000, target_censor_rate = 0.30)
head(sim_df)
cat("Achieved Censor Rate:", round(mean(sim_df$status == 0) * 100, 2), "%\n")

# Fitting Cox Proportional Hazards Model
cox_model <- coxph(Surv(time, status) ~ treatment, data = sim_df)
summary(cox_model)

#----------------------------------------------------------------
# Monte Carlo Power Sensitivity Analysis Across Censoring Rates
#----------------------------------------------------------------
library(survival)

run_mc_power_r <- function(sample_sizes = c(50, 80, 100, 150, 200, 250, 300, 500),
                           n_sim = 1000,
                           seed = 42) {
    set.seed(seed)
    lam <- 0.04
    nu <- 1.4
    beta <- -0.7
    tau <- 22.07

    results <- list()

    for(N in sample_sizes) {
        p_values <- numeric(n_sim)

        for(s in 1:n_sim) {
            X <- rbinom(N, 1, 0.5)
            U_t <- runif(N)
            T_true <- (-log(U_t) / (lam * exp(X * beta)))^(1 / nu)
            C_time <- tau * runif(N)
            Y_obs <- pmin(T_true, C_time)
            d_event <- as.numeric(T_true <= C_time)

            # Cox Proportion Hazards Model fit
            fit <- coxph(Surv(Y_obs, d_event) ~ X)
            p_values[s] <- summary(fit)$coefficients["X", "Pr(>|z|)"]
        }
        
        power <- mean(p_values < 0.05)

        results[[as.character(N)]] <- data.frame(
            Sample_Size = N,
            Target_HR = round(exp(beta), 2),
            Censor_Rate = "30.0%",
            Simulations = n_sim,
            Empirical_Power = paste0(round(power * 100, 1), "%")
        )
    }

    do.call(rbind, results)
}
# Run R simulation
power_table <- run_mc_power_r()
print(power_table)
