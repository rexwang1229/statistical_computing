#Ex1: Sampling from finite populations
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

