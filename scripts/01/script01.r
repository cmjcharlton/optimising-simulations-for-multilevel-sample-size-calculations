# Required packages
library(MASS)
library(ggplot2)

RNGkind("default")

# Initial inputs
set.seed(1)
siglevel <- 0.025
n1 <- 70
z1score <- abs(qnorm(siglevel))
beta <- 3
effectbeta <- abs(beta)
sigmae <- sqrt(81)

# Inputs for model fitting
data <- list(y = NULL, x0 = NULL)
modelformula <- y ~ 1

# Initial input for power in two approaches
nreps <- 1001

cat("     The programme was executed at", date(), "\n")
cat("--------------------------------------------------------------------\n")

output <- data.frame(
  n = NA,
  zpb0 = NA,
  spb0 = NA
)

# This part is constant through simulations
data$x0 <- rep(1, n1)
fixpart <- data$x0 * beta

rowcount <- 1
for (simus in c(seq(50, 950, 50), seq(1000, 5000, 500))) {
  cat("Start of simulation for simulation sizes of ",
      simus,
      " units\n")
  for (reps in seq_len(nreps)) {

    estpower <- rep(0, simus)
    sdepower <- rep(0, simus)
    powaprox <- 0
    
    # Simulation step
    for (iter in seq_len(simus)) {
      e <- rnorm(n1, 0, sigmae)
      data$y <- fixpart + e

      # Fitting the model using lm function
      fitmodel <- lm(modelformula, data)
      
      # To obtain the power of parameter(s)
      estbeta <- coef(fitmodel)
      sdebeta <- sqrt(diag(vcov(fitmodel)))
      sgnbeta <- sign(estbeta)
      cibeta <- estbeta - sgnbeta * z1score * sdebeta
      powaprox <- powaprox + as.integer(estbeta * cibeta > 0)
      sdepower[iter] <- sdebeta
      estpower[iter] <- estbeta / sdebeta
    }
    
    # Powers and their CIs
    meanaprox <- powaprox / simus
    estsde <- sqrt(mean(sdepower^2))
    powsde <- pnorm(effectbeta / estsde - z1score)
    
    output[rowcount, ] <- c(simus, meanaprox, powsde)
    rowcount <- rowcount + 1
    
  }
}

quants <- c(0.05, 0.5, 0.95)
zresults <- data.frame(do.call("rbind", by(output, output$n, function(x) {
  zquant <- quantile(x$zpb0, probs=quants)
  return(c(simulations = x$n[1], lower = zquant[[1]], power = mean(x$zpb0), upper = zquant[[3]]))
})))
zresults$power_method <- "0/1"

sresults <- data.frame(do.call("rbind", by(output, output$n, function(x) {
  squant <- quantile(x$spb0, probs=quants)
  return(c(simulations = x$n[1], lower = squant[[1]], power = mean(x$spb0), upper = squant[[3]]))
})))
sresults$power_method <- "SE"

outputlines <- rbind(zresults, sresults)
outputlines$power_method <- factor(outputlines$power_method, levels=c("0/1", "SE"))

# Figure 1
ggplot(
  data = outputlines,
  aes(
    x = simulations,
    y = power,
    ymin = lower,
    ymax = upper,
    fill = power_method
  )
) +
  geom_line(aes(colour = power_method)) +
  geom_ribbon(alpha = 1) +
  scale_fill_grey() +
  labs(colour = "method", fill = "method") +
  theme_bw()

ggsave("figure1.svg", width=13, height=10, units="cm")
