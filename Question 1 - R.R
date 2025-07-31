#PART 01 A

# Load necessary libraries
library(ggplot2)

# Define the probability density function f(x)
f <- function(x) {
  return(1/2 * exp(-abs(x)))
}

# Metropolis-Hastings algorithm
metropolis_hastings <- function(N, s, x0) {
  samples <- c(x0)  # Initial value
  x_prev <- x0  # Set the previous sample to the initial value
  
  for (i in 2:N) {
    x <- rnorm(1, mean = x_prev, sd = s)
    
    # Compute the ratio
    log_r <- log(f(x)) - log(f(x_prev))
    log_u <- log(runif(1))
    
    if (log_u < log_r) {
      samples <- c(samples, x)
      x_prev <- x
    } else {
      samples <- c(samples, x_prev)
    }
  }
  
  return(samples)
}

# Parameters
N <- 10000
s <- 1
x0 <- 1  # Initial value

# Generate samples
samples <- metropolis_hastings(N, s, x0)

# Create a data frame for plotting
df <- data.frame(samples)

# Plot histogram and kernel density plot
ggplot(df, aes(x = samples)) +
  geom_histogram(bins = 50, fill = 'green', alpha = 0.5, aes(y = ..density..)) +
  geom_density(color = 'darkblue') +
  stat_function(fun = f, color = 'red', size = 1) +
  xlim(-10, 10) +
  theme_minimal()


# Calculate sample mean and standard deviation
sample_mean <- mean(samples)
sample_std <- sd(samples)

print(paste("Sample Mean:", sample_mean))
print(paste("Sample Standard Deviation:", sample_std))


#PART 01 B

library(coda)

N <- 2000
s_values <- seq(0.001, 1, length.out = 250)
J <- 4

compute_R <- function(N, s, J) {
  chains <- replicate(J, metropolis_hastings(N, s, runif(1, 0, 1)))
  Mj <- apply(chains, 2, mean)
  Vj <- apply(chains, 2, var)
  W <- mean(Vj)
  M <- mean(Mj)
  B <- N / (J - 1) * sum((Mj - M)^2)
  R <- sqrt((W + B) / W)
  return(R)
}

R_values <- sapply(s_values, function(s) compute_R(N, s, J))

print(R_values)

plot(s_values, R_values, type = "l", xlab = "s", ylab = "R value", main = "Convergence (R_value) vs. s", col = "blue")
grid()


