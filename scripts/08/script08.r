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
n1 <- 19
n2 <- 148
randsize <- 2
totalsamplelow <- 4
totalsamplehigh <- 20
totalsamplestep <- 2
xctable <- matrix(scan("fife2.txt"), n1, n2, byrow=TRUE)
beta <- 3.0
effectbeta <- abs(beta)
sgnbeta <- sign(beta)
sigma2u <- matrix(9.000, 1, 1)
sigma2v <- matrix(36.000, 1, 1)
sigmae <- sqrt(64.000)
Tsamplerange <- seq(totalsamplelow, totalsamplehigh, totalsamplestep)
totalsize <- length(Tsamplerange)

output <-
  data.frame(
    XC2 = rep(n2, totalsize),
    XC1 = rep(n1, totalsize),
    Tninrow = Tsamplerange,
    spb0 = NA
  )

rowcount <- 1

# Creating grouped data
modelformula <- y ~ 1 + (1|l1id) + (1|l2id)
data <- vector("list", 4)
names(data) <- c("l1id", "l2id", "y", "x0")

cat("     The programme was executed at", date(),"\n")
cat("--------------------------------------------------------------------\n")

cat("Start of simulation for ", n2, " XC2, ", n1, " XC1 and 1 1-level units\n")
for (totalsample in Tsamplerange) {
  cat("Start of simulation for ", totalsample, " observations per XC1 unit\n")
  sdepower <- matrix(0, 1, simus)

  # Simulation step
  for (iter in 1:simus) {
    if (iter/10 == floor(iter/10)) {
      cat(" Iteration remain=", simus-iter, "\n")
    }

    # conditional sampling from XC2 given XC1
    xcunbal <- rmultinom(1, totalsample, xctable[1, ])
    for (i in 2:n1) {
      xcunbal <- c(xcunbal, rmultinom(1, totalsample, xctable[i, ]))
    }
    # Inputs for model fitting
    data$l1id <- rep(rep(1:n1, each=n2), xcunbal)
    data$l2id <- rep(rep(1:n2, n1), xcunbal) 
    cumxc <- c(0, cumsum(xcunbal))
    nlen <- cumxc[length(cumxc)]
    data$x0 <- matrix(1, nlen, 1)
    z1 <- matrix(1, nlen, 1)
    z2 <- matrix(1, nlen, 1)
    u <- mvrnorm(n1, 0, sigma2u)
    v <- mvrnorm(n2, 0, sigma2v)
    fixpart <- data$x0*beta
    rand1part <- rowSums(z1*u[data$l1id, ])
    rand2part <- rowSums(z2*v[data$l2id, ])
    e <- rnorm(nlen, 0, sigmae)
    data$y <- fixpart+rand1part+rand2part+e

    # Fitting the model using lmer funtion
    fitmodel <- lmer(modelformula, data, REML = TRUE)

    # To obtain the power of parameter(s)
    sdepower[iter] <- sqrt(diag(vcov(fitmodel)))
  }

  estsde <- sqrt(mean(sdepower^2))
  powsde <- pnorm(effectbeta / estsde - z1score)
  
  output[rowcount, "spb0"] <- powsde

  rowcount <- rowcount+1
  cat("--------------------------------------------------------------------\n")
}

calc_c <- function(n, p, r, step) {
  fife <- data.frame(n, p)
  fife$sqrtn <- sqrt(n)
  fife$zscore <- qnorm(p)

  d <- seq(from = 0, to = r, by = step)
  nvals <- length(d)
  ss <- rep(0, nvals)
  for (dval in seq(1, nvals)) {
    fife$x <- sqrt(fife$n / (1 + d[dval] * fife$n))
    ss[dval] <- summary(lm(fife, formula = "zscore ~ x"))$sigma
  }
  out <- data.frame(d, ss)
  fife$x <- sqrt(fife$n / (1 + out[which.min(out$ss), "d"] * fife$n))
  return(list(fife = fife, out = out))
}

result <- calc_c(output$Tninrow, output$spb0, 0.5, 0.005)
out <- result$out
fife <- result$fife

# Figure 10a
figure10a <- ggplot(data = out, aes(d, y = ss)) +
  geom_line(colour = "black") +
  xlab("c") +
  ylab("SSD") +
  theme_bw()

# Figure 10b
cutoff <- qnorm(0.8)
figure10b <- ggplot(data = fife, aes(x = x, y = zscore)) +
  geom_point(colour = "black") +
  geom_smooth(method = "loess", se = FALSE, colour = "grey50") +
  geom_smooth(method = 'lm', se = FALSE, colour = "grey75") +
  geom_line(y = cutoff, colour = "grey25") +
  xlab(expression(f(n))) +
  ylab("z-score") +
  theme_bw()

figure10 <- gridExtra::arrangeGrob(figure10a, figure10b, ncol = 2)
ggsave("figure10.svg", plot=figure10, width=16, height=8, units="cm")