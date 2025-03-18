# Required packages
library(MASS)
library(boot)
library(lme4)
library(ggplot2)
library(patchwork)

RNGkind("default")

# Initial inputs
set.seed(1)
siglevel <- 0.025
z1score <- abs(qnorm(siglevel))
simus <- 5000
n1 <- 20
n2low <- 10
n2high <- 70
n2step <- 5
beta <- 0.40547
effectbeta <- abs(beta)
sigma2u <- matrix(0.5, 1, 1)
n2range <- seq(n2low, n2high, n2step)
n2size <- length(n2range)
output <-
  data.frame(
    N = n2range,
    n = rep(n1, n2size),
    zpb0 = NA,
    spb0 = NA
  )
rowcount <- 1

# Inputs for model fitting
modelformula <- y ~ 1 + (1 | l2id)
data <- list(l2id = NULL, y = NULL, x0 = NULL)

cat("     The programme was executed at", date(), "\n")
cat("--------------------------------------------------------------------\n")

for (n2 in seq(n2low, n2high, n2step)) {
  nlen <- n1 * n2
  z <- rep(1, nlen)
  data$x0 <- rep(1, nlen)
  data$l2id <- rep(1:n2, each = n1)
  sdepower <- rep(0, simus)
  powaprox <- 0
  
  cat(" Start of simulation for sample sizes of ",
      n1,
      " micro and ",
      n2,
      "macro units\n")

  for (iter in 1:simus) {
    if (iter / 10 == floor(iter / 10)) {
      cat(".")
    }

    # Inputs for model fitting
    u <- mvrnorm(n2, 0, sigma2u)
    fixpart <- data$x0 * beta
    randpart <- z * u[data$l2id]
    binomprob <- inv.logit(fixpart + randpart)
    data$y <- rbinom(nlen, 1, binomprob)

    # Fitting the model using lmer function
    fitmodel <-
        glmer(
          modelformula,
          data,
          family = binomial(link = "logit"),
          nAGQ = 3
        )
    
    # To obtain the power of parameter(s)
    estbeta <- fixef(fitmodel)
    sdebeta <- sqrt(diag(vcov(fitmodel)))
    sgnbeta <- sign(estbeta)
    cibeta <- estbeta - sgnbeta * z1score * sdebeta
    powaprox <- powaprox + as.integer(estbeta * cibeta > 0)
    sdepower[iter] <- sdebeta
  }
  
  # Powers and their CIs
  meanaprox <- powaprox / simus
  estsde <- sqrt(mean(sdepower^2))
  powsde <- pnorm(effectbeta / estsde - z1score)
  
  output[rowcount, c("zpb0", "spb0")] <- c(meanaprox, powsde)
  
  rowcount <- rowcount + 1
  cat("\n")
}

output$power <- output$spb0
output$zscore <- qnorm(output$power)
output$sqrtN <- sqrt(output$N)
output$power2 <- output$zpb0
output$zscore2 <- qnorm(output$power2)

# Figure 13a
cutoff <- 0.8
figure13a <- ggplot(data = output, aes(x = N)) +
  geom_point(aes(y = power, colour = "SE")) +
  geom_smooth(aes(y = power, colour = "SE"),
              method = "loess",
              se = FALSE) +
  geom_point(aes(y = power2, colour = "0/1")) +
  geom_smooth(aes(y = power2, colour = "0/1"),
              method = "loess",
              se = FALSE) +
  geom_line(y = cutoff, colour = "grey25") +
  scale_colour_manual(name = "method", values = c("black", "grey75")) +
  theme_bw() +
  theme(legend.position = "bottom")

# Figure 13b
cutoff <- qnorm(0.8)
figure13b <- ggplot(data = output, aes(x = sqrtN)) +
  geom_point(aes(y = zscore, colour = "SE")) +
  geom_smooth(aes(y = zscore, colour = "SE"),
              method = "loess",
              se = FALSE) +
  geom_point(aes(y = zscore2, colour = "0/1")) +
  geom_smooth(aes(y = zscore2, colour = "0/1"),
              method = "loess",
              se = FALSE) +
  geom_line(y = cutoff, colour = "grey25") +
  scale_colour_manual(name = "method", values = c("black", "grey75")) +
  xlab(expression(sqrt(N))) +
  ylab("z-score") +
  theme_bw() + 
  theme(legend.position = "bottom")

figure13 <- figure13a + figure13b + plot_layout(guide = "collect") & theme(legend.position = "bottom")
ggsave("figure13.svg", plot=figure13, width=16, height=10, units="cm")