# Required packages
library(MASS)
library(ggplot2)

RNGkind("default")

# Initial inputs
set.seed(1)
siglevel <- 0.025
z1score <- abs(qnorm(siglevel))
simus <- 1000
n1low <- 10
n1high <- 100
n1step <- 5
beta <- 3
effectbeta <- abs(beta)
sigmae <- sqrt(81)
n1range <- seq(n1low, n1high, n1step)
n1size <- length(n1range)
rowcount <- 1

# Inputs for model fitting
data <- list(y = NULL, x0 = NULL)
modelformula <- y ~ 1

cat("     The programme was executed at", date(), "\n")
cat("--------------------------------------------------------------------\n")

output <- data.frame(
  n = n1range,
  spb0 = NA,
  spb200 = NA,
  spb50 = NA
)

for (n1 in seq(n1low, n1high, n1step)) {
  cat(" Start of simulation for sample sizes of ", n1, " units\n")
  data$x0 <- rep(1, n1)
  sdepower <- rep(0, simus)
  
  # Simulation step
  for (iter in seq_len(simus)) {
    if (iter / 10 == floor(iter / 10)) {
      cat(".")
    }
    fixpart <- data$x0 * beta
    e <- rnorm(n1, 0, sigmae)
    data$y <- fixpart + e
    
    # Fitting the model using lm function
    fitmodel <- lm(modelformula, data)
    
    # To obtain the power of parameter(s)
    sdepower[iter] <- sqrt(diag(vcov(fitmodel)))
  }
  
  estsde <- sqrt(mean(sdepower^2))  
  estsde50 <- sqrt(mean(sdepower[1:50]^2))  
  estsde200 <- sqrt(mean(sdepower[1:200]^2))  
  
  powsde <- pnorm(effectbeta / estsde - z1score)
  powsde50 <- pnorm(effectbeta / estsde50 - z1score)
  powsde200 <- pnorm(effectbeta / estsde200 - z1score)
  
  output[rowcount, c("spb0", "spb200", "spb50")] <- c(powsde, powsde200, powsde50)
  
  rowcount <- rowcount + 1
  cat("\n")
}

output$power <- output$spb0

zeropointeight <- predict(loess(n ~ power, data = output), newdata = list(power = 0.8))

# Figure 2
ggplot(data = output, aes(x = n, y = power)) +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linewidth = 0.8) +
  geom_point(colour = "black") +
  geom_hline(yintercept = 0.8, linetype = "solid", colour = "grey25") +
  geom_vline(xintercept = zeropointeight, linetype = "solid", colour = "grey25") +
  theme_bw()

ggsave("figure2.svg", width=10, height=10, units="cm")

output$zscore <- qnorm(output$power)
output$zscore200 <- qnorm(output$spb200)
output$zscore50 <- qnorm(output$spb50)
output$sqrtn <- sqrt(output$n)

# Figure 3
ggplot(data = output) +
  geom_smooth(aes(x = sqrtn, y = zscore), method = "loess", se = FALSE, colour = "black") +
  geom_point(aes(x = sqrtn, y = zscore, colour = "1000")) +
  geom_point(aes(x = sqrtn, y = zscore200, colour = "200")) +
  geom_point(aes(x = sqrtn, y = zscore50, colour = "50")) +
  scale_colour_grey(name = "simulations", start = 0, end = .9) +
  xlab(expression(sqrt(n))) +
  ylab("z-score") + 
  theme_bw()

ggsave("figure3.svg", width=13, height=10, units="cm")
