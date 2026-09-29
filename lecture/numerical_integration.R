# Trapezodial Rule
trapezoid <- function(ftn, a, b, n = 100) {
    # ftn is a function of a single variable
    # we assume a < b and n is a positive integer
    h <- (b - a) / n
    x.vec <- seq(a, b, by = h)
    f.vec <- sapply(x.vec, ftn)
    T <- h * (f.vec[1] / 2 + sum(f.vec[2:n]) + f.vec[n + 1] / 2)
    return(T)
}

ftn6 <- function(x) return(4 * x^3)
trapezoid(ftn6, 0, 1, n = 20)
trapezoid(ftn6, 0, 1, n = 40)
trapezoid(ftn6, 0, 1, n = 60)


# Simpson's Rule
simpson_n <- function(ftn, a, b, n = 100) {
    # ftn is a function of a single variable
    # we assume a < b and n is a positive even integer
    n <- max(c(2 * (n %/% 2), 4))
    h <- (b - a) / n
    x.vec1 <- seq(a + h, b - h, by = 2 * h) # x_1, x_3, ...x_n-1 odds
    x.vec2 <- seq(a + 2 * h, b - 2 * h, by = 2 * h) # x_2, x_4, x_n-2... even
    f.vec1 <- sapply(x.vec1, ftn)
    f.vec2 <- sapply(x.vec2, ftn)
    S <- h / 3 * (ftn(a) + ftn(b) + 4 * sum(f.vec1) + 2 * sum(f.vec2))
    return(S)
}

ftn6 <- function(x) return(4 * x^3)
simpson_n(ftn6, 0, 1, 20)

# Simpson's Rule with a set tolerance
simpson <- function(ftn, a, b, tol = 1e-8, verbose = FALSE) {
    # numerical integral of ftn from a to b
    # using Simpson's rule with tolerance tol
    n <- 4
    h <- (b - a) / 4
    fx <- sapply(seq(a, b, by = h), ftn)
    S <- sum(fx * c(1, 4, 2, 4, 1)) * h / 3
    S.diff <- tol + 1 # ensures we loop at least once

    # increase n until S changes by less than tol
    while(S.diff > tol) {
        cat('n = ', n, 'S = ', S, '\n')
        S.old <- S
        n <- 2 * n
        h <- h / 2
        fx[seq(1, n + 1, by = 2)] <- fx
        fx[seq(2, n, by = 2)] <- sapply(seq(a + h, b - h, by = 2 * h), ftn)
        S <- h / 3 * (fx[1] + fx[n + 1] + 4 * sum(fx[seq(2, n, by = 2)])+ 
        2 * sum(fx[seq(3, n - 1, by = 2)]))
        S.diff <- abs(S - S.old)
    }
    if(verbose) cat('partition size', n, '\n')
    return(S)
}
