# Required packages
library(MASS)
library(lme4)
library(ggplot2)

RNGkind("default")

# Initial inputs
set.seed(1)
siglevel <- 0.025
z1score <- abs(qnorm(siglevel))
simus <- 1000
n1 <- 5
n2 <- 4
n3low <- 10
n3high <- 60
n3step <- 5
beta <- 2.5
effectbeta <- abs(beta)
rand3column <- 1
rand2column <- 1
sigma2u <- matrix(16, 1, 1)
sigma2v <- matrix(16, 1, 1)
sigmae <- sqrt(64)
n3range <- seq(n3low, n3high, n3step)
n3size <- length(n3range)

output <-
  data.frame(
    n3 = n3range,
    n2 = rep(n2, n3size),
    n1 = rep(n1, n3size),
    spb0 = NA
  )

rowcount <- 1

# Creating grouped data
modelformula <- y ~ 1 + (1 | l2id:l3id) + (1 | l3id)
data <- list(l2id = NULL,
             l3id = NULL,
             y = NULL,
             x0 = NULL)

cat("     The programme was executed at", date(), "\n")
cat("--------------------------------------------------------------------\n")

# Sample size combination
for (n3 in seq(n3low, n3high, n3step)) {
  cat("Simulation for ",
      n1,
      ", 1-level ",
      n2,
      " 2-level and ",
      n3,
      " 3-level units\n")
  nlen <- n1 * n2 * n3
  y <- rep(0, nlen)
  z2 <- rep(1, nlen)
  z3 <- rep(1, nlen)
  data$x0 <- rep(1, nlen)
  data$l2id <- rep(1:(n2 * n3), each = n1)
  data$l3id <- rep(1:n3, each = n1 * n2)
  sdepower <- rep(0, simus)
  
  # Simulation step
  for (iter in 1:simus) {
    if (iter / 10 == floor(iter / 10)) {
      cat(".")
    }

    # Inputs for model fitting
    e <- rnorm(nlen, 0, sigmae)
    u <- mvrnorm(n2 * n3, 0, sigma2u)
    v <- mvrnorm(n3, 0, sigma2v)
    fixpart <- data$x0 * beta
    rand2part <- z2 * u[data$l2id]
    rand3part <- z3 * v[data$l3id]
    data$y <- fixpart + rand2part + rand3part + e
    
    # Fitting the model using lmer function
    fitmodel <- lmer(modelformula, data, REML = TRUE)
    
    # To obtain the power of parameter(s)
    sdepower[iter] <- sqrt(diag(vcov(fitmodel)))
  }
  
  estsde <- sqrt(mean(sdepower^2))  
  powsde <- pnorm(effectbeta / estsde - z1score)
  
  output[rowcount, "spb0"] <- powsde
  
  rowcount <- rowcount + 1
  cat("\n")
}

output$power <- output$spb0
output$zscore <- qnorm(output$power)
output$sqrtn3 <- sqrt(output$n3)

# Figure 8a
cutoff <- 0.8
figure8a <- ggplot(data = output, aes(x = n3, y = power)) +
  geom_point(colour = "black") +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 0.8) +
  geom_line(y = cutoff, colour = "grey25") +
  xlab(expression(n[3])) +
  theme_bw()

# Figure 8b
cutoff <- qnorm(0.8)
figure8b <- ggplot(data = output, aes(x = sqrtn3)) +
  geom_point(aes(y = zscore), colour = "black") +
  geom_smooth(aes(y = zscore), method = "loess", se = FALSE, colour = "black", linewidth = 0.8) +
  geom_line(y = cutoff, colour = "grey25") +
  xlab(expression(sqrt(n[3]))) + 
  ylab("z-score") +
  theme_bw()

figure8 <- gridExtra::arrangeGrob(figure8a, figure8b, ncol = 2)
ggsave("figure8.svg", plot=figure8, width=16, height=8, units="cm")