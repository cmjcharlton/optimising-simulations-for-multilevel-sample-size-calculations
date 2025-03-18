# Required packages
library(MASS)
library(lme4)
library(ggplot2)

RNGkind("default")

# Initial inputs
set.seed(1)
siglevel <- 0.025
z1score <- abs(qnorm(siglevel))
simus <- 5000
n1low <- 2
n1high <- 15
n1step <- 1
n2 <- 4
n3 <- 30
beta <- 2.5
effectbeta <- abs(beta)
sigma2u <- matrix(16, 1, 1)
sigma2v <- matrix(16, 1, 1)
sigmae <- sqrt(64)
n1range <- seq(n1low, n1high, n1step)
n1size <- length(n1range)
output <-
  data.frame(
    n3 = rep(n3, n1size),
    n2 = rep(n2, n1size),
    n1 = n1range,
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
for (n1 in seq(n1low, n1high, n1step)) {
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
    if (iter / 1000 == floor(iter / 1000)) {
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
output$sqrtn1 <- sqrt(output$n1)
output$function_of_n1 <- sqrt(output$n1 / (64 + 80 * output$n1))

# Figure 9b
cutoff <- qnorm(0.8)
ggplot(data = output, aes(x = function_of_n1, y = zscore)) +
  geom_point(colour = "black") +
  geom_smooth(method = "loess", se = FALSE, colour = "grey50") +
  geom_smooth(method = 'lm', se = FALSE, colour = "grey75") +
  geom_line(y = cutoff, colour = "grey25") +
  xlab(expression(f(n[1]))) + 
  ylab("z-score") + 
  theme_bw()

ggsave("figure9b.svg", width=8, height=8, units="cm")
