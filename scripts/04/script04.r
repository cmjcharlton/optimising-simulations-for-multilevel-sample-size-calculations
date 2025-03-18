# Required packages
library(MASS)
library(lme4)
library(ggplot2)
library(doParallel)

RNGkind("L'Ecuyer-CMRG")

# Initial inputs
siglevel <- 0.025
z1score <- abs(qnorm(siglevel))
simus <- 1000
saveoutput <- c(50, 200, simus)
seeds <- seq(1, 100)
n1 <- 20
n2low <- 10
n2high <- 50
n2step <- 5
beta <- 2.5
effectbeta <- abs(beta)
sigma2u <- matrix(16, 1, 1)
sigmae <- sqrt(81)
n2range <- seq(n2low, n2high, n2step)
n2size <- length(n2range)

# Inputs for model fitting
data <- list(l2id = NULL, y = NULL, x0 = NULL)
modelformula <- y ~ 1 + (1 | l2id)

cat("     The programme was executed at", date(), "\n")
cat("--------------------------------------------------------------------\n")

cl <- makeCluster(detectCores(logical = FALSE))
registerDoParallel(cl)

r <- foreach(seed = seeds, .packages = c("MASS", "lme4")) %dopar% {
  outputs <- list()
  set.seed(seed)
  for (iter in saveoutput) {
    # Set up output data frame
    dfname <- paste0(iter, "_", seed)
    outputs[[dfname]] <-
      data.frame(N = n2range,
                 n = rep(n1, n2size))
    outputs[[dfname]]$spb0 <- NA
  }
  rowcount <- 1
  for (n2 in seq(n2low, n2high, n2step)) {
    nlen <- n1 * n2
    data$x0 <- rep(1, nlen)
    z <- matrix(1, nlen, 1)
    data$l2id <- rep(1:n2, each = n1)
    sdepower <- rep(0, simus)
    for (iter in 1:simus) {
      # Inputs for model fitting
      e <- rnorm(nlen, 0, sigmae)
      u <- mvrnorm(n2, 0, sigma2u)
      fixpart <- data$x0 * beta
      randpart <- rowSums(z * u[data$l2id, ])
      data$y <- fixpart + randpart + e

      # Fitting the model using lmer function (ML)
      fitmodel <- lmer(modelformula, data, REML = TRUE)

      # To obtain the power of parameter(s) (se method)
      sdepower[iter] <- sqrt(diag(vcov(fitmodel)))

      if (iter %in% saveoutput) {
        dfname <- paste0(iter, "_", seed)
        estsde <- sqrt(mean(sdepower[1:iter]^2))
        powsde <- pnorm(effectbeta / estsde - z1score)
        outputs[[dfname]][rowcount, "spb0"] <- powsde
      }
    }
    rowcount <- rowcount + 1
  }
  return(outputs)
}
stopCluster(cl)

outputs <- unlist(r, recursive = FALSE)

powset <- 0.8
cutoff <- qnorm(powset)

results <- NULL
for (seed in seeds) {
  for (iter in saveoutput) {
    outname <- paste0(iter, "_", seed)
    answers <- matrix(NA, 1, 6)
    colnames(answers) <-
      c("iter", "seed", "clusters", paste0("scenario", 1:3))
    rownames(answers) <- n1
    
    answers[, 1] <- iter
    answers[, 2] <- seed
    answers[, 3] <- n1
        
    # 1. Interpolate from SE
    
    row1 <-
      tail(subset(outputs[[outname]], spb0 < powset), 1)
    row2 <-
      head(subset(outputs[[outname]], spb0 >= powset), 1)
    if (nrow(row1) > 0 && nrow(row2) > 0) {
      step <- (row2$spb0 - row1$spb0) / n2step
      diff <- powset - row1$spb0
      inc <- ceiling(diff / step)
      answers[1, 4] <- row1$N + inc
    }
    
    # 2. Predict from 10 + 50 clusters (SE method)
    try({
      subdata <- subset(outputs[[outname]], N == 10 | N == 50)
      msub <- lm(qnorm(spb0) ~ sqrt(N), data = subdata)
      preddata <-
        data.frame(N = seq(1, 100),
                   n = rep(n1, 100))
      pred <- predict(msub, newdata = preddata)
      preddata[["pred"]] <- pred
      answers[1, 5] <-
        preddata$N[which.max(preddata$pred > cutoff)]
    })
    
    # 3. Predict from all clusters (SE method)
    try({
      msub <- lm(qnorm(spb0) ~ sqrt(N), data = outputs[[outname]])
      preddata <-
        data.frame(N = seq(1, 100),
                   n = rep(n1, 100))
      pred <- predict(msub, newdata = preddata)
      preddata[["pred"]] <- pred
      answers[1, 6] <-
        preddata$N[which.max(preddata$pred > cutoff)]
    })
    print(answers)
    results <- rbind(results, data.frame(answers))
  }
}

resultslong <-
  reshape(
    results,
    direction = "long",
    timevar = "scenario",
    varying = paste0("scenario", 1:3),
    v.names = "Ncluster"
  )
# values of one are unlikely to be valid, so remove them
resultslong$Ncluster[resultslong$Ncluster == 1] <- NA

resultslong$scenario <-
  factor(
    resultslong$scenario,
    labels = c(
      "interp",
      "10+50 pred",
      "all pred"
    )
  )

# Figure 7
ggplot(resultslong, aes(x = Ncluster)) + 
  geom_histogram(binwidth = 1) + 
  facet_grid(
    rows = vars(scenario),
    cols = vars(iter),
    labeller = label_bquote(
      cols = sims==.(formatC(iter, format = "d", big.mark = ","))
    )
  ) +
  scale_x_continuous(limits = c(10, 90), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  xlab("no. of clusters") +
  ylab("count") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    strip.background  = element_rect(
      fill = "grey95", colour = "grey50",
      linetype = 1, linewidth = rel(0.8)
    ),
    panel.border = element_rect(
      colour = "grey50",
      fill = NA, linetype = 1,
      linewidth = rel(0.8)
    ),
    panel.spacing = unit(0.7, "lines", data = NULL)
  )

ggsave("figure7.svg", width=16, height=16, units="cm")

# Way of getting numbers from graphs
table(resultslong[resultslong$iter == 1000 &
                    resultslong$scenario == "interp",]$Ncluster)
table(resultslong[resultslong$iter == 1000 &
                    resultslong$scenario == "10+50 pred",]$Ncluster)
table(resultslong[resultslong$iter == 1000 &
                    resultslong$scenario == "all pred",]$Ncluster)

table(resultslong[resultslong$iter == 200 &
                    resultslong$scenario == "interp",]$Ncluster)
table(resultslong[resultslong$iter == 200 &
                    resultslong$scenario == "10+50 pred",]$Ncluster)
table(resultslong[resultslong$iter == 200 &
                    resultslong$scenario == "all pred",]$Ncluster)

table(resultslong[resultslong$iter == 50 &
                    resultslong$scenario == "interp",]$Ncluster)
table(resultslong[resultslong$iter == 50 &
                    resultslong$scenario == "10+50 pred",]$Ncluster)
table(resultslong[resultslong$iter == 50 &
                    resultslong$scenario == "all pred",]$Ncluster)
