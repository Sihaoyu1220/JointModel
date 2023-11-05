# Figure 1: censoring and measurement error 
par(mfrow = c(1, 2))
set.seed(2023)

# sample size n = 20
n = 30
# Generate random x
x <- runif(n, 0, 10)

# Generate random y = - 10 * x + 100, where e ~ N(0, 5)
e <- rnorm(n, 0, 3)
y <- - 3 * x + 60 + e

# censoring
# impute the response value to half the detection limit if the value is < 6.
dl <- 40
y_censored <- ifelse(y < dl, 0.5 * dl, y)

# plot the least square line for imputed y
plot(x, y_censored, pch=16, ylab = "y", ylim = c(20, 65))
mod3 <- lm(y_censored ~ x)
abline(coefficients(mod3), lwd=2, lty = 2)

# Plot the least square line for true y
points(x, y)
mod1 <- lm(y ~ x)
abline(coefficients(mod1), lwd=2, lty = 1)
abline(h = dl, lwd = 2, lty = 3)

legend("bottomleft",c("naive method","real line","observed values", "unobserved values"), lty=c(2,1,0,0),pch=c(NA,NA,16,1), cex=0.8, bty = "n")


# measurement error
# Assume x is measured with error, the observed x* = x + eps, where eps ~ N(0, 1.5)
eps <- rnorm(n, 0, 1.5)
x_star <- x + eps

# Plot the least square line for x without measurement error
plot(x, y, pch=16, ylim = c(20, 65))
abline(coefficients(mod1), lwd=2, lty = 1)

# Plot the least square line for x with measurement error
points(x_star, y)
mod2 <- lm(y ~ x_star)
abline(coefficients(mod2), lwd=2, lty = 2)

legend("bottomleft",c("with measurement error","without measurement error"),lty=c(2,1),pch=c(1, 16), cex=.8, bty = "n")





