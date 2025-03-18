library(ggplot2)
library(patchwork)

RNGkind("default")

calcpint <- function(start, step, end, n2min, n2step, n2max, RESID, TAU2) {
  N = rep(seq(n2min, n2max, n2step), times = ((end - start) / step) + 1)
  n = rep(seq(start, end, step), each = ((n2max - n2min) / n2step) + 1)
  output <- data.frame(
    total = N * n,
    N = N,
    n = n,
    se = sqrt((TAU2 + RESID/n) /N)
  )
  return(output)
}

pint <- calcpint(10, 5, 60, 10, 2, 50, 81, 16)
pint$zs <- 2.5 / pint$se
siglevel <- 0.025
z1score <- abs(qnorm(siglevel))
pint$power <- pnorm(pint$zs - z1score)
pint$zscore <- qnorm(pint$power)
pint$sqrtn <- sqrt(pint$n)
pint$sqrtN <- sqrt(pint$N)

cutoff <- 0.8

# Figure 4a
figure4a <- ggplot(data = pint, aes(x = N, y = power, group = n)) +
  geom_point(aes(colour = n)) +
  geom_line(aes(colour = n)) +
  geom_line(y = cutoff, colour = "grey25") +
  scale_color_distiller(type = "seq", direction = -1, palette = "Greys") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("figure4a.svg", width=8, height=10, units="cm")

# Figure 5a
figure5a <- ggplot(data = subset(pint, N %in% seq(from=10, to=50, by=4)), aes(x = n, y = power, group = N)) +
  geom_point(aes(colour = N)) +
  geom_line(aes(colour = N)) +
  geom_line(y = cutoff, colour = "grey25") +
  scale_color_distiller(type = "seq", direction = -1, palette = "Greys") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("figure5a.svg", width=8, height=10, units="cm")

cutoff <- qnorm(0.8)

# Figure 4b
figure4b <- ggplot(data = pint , aes(x = sqrtN, y = zscore, group = n)) +
  geom_point(aes(colour = n)) +
  geom_line(aes(colour = n)) +
  geom_line(y = cutoff, colour = "grey25") +
  xlab(expression(sqrt(N))) +
  ylab("z-score") +
  scale_color_distiller(type = "seq", direction = -1, palette = "Greys") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("figure4b.svg", width=8, height=10, units="cm")

# Figure 5b
figure5b <- ggplot(data = subset(pint, N %in% seq(from=10, to=50, by=4)), aes(x = sqrtn, y = zscore, group = N)) +
  geom_point(aes(colour = N)) +
  geom_line(aes(colour = N)) +
  geom_line(y = cutoff, colour = "grey25") +
  xlab(expression(sqrt(n))) +
  ylab("z-score") +
  scale_color_distiller(type = "seq", direction = -1, palette = "Greys") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("figure5b.svg", width=8, height=10, units="cm")

figure4 <- figure4a + figure4b + plot_layout(guide = "collect") & theme(legend.position = "bottom")
ggsave("figure4.svg", plot=figure4, width=16, height=10, units="cm")

figure5 <- figure5a + figure5b + plot_layout(guide = "collect") & theme(legend.position = "bottom")
ggsave("figure5.svg", plot=figure5, width=16, height=10, units="cm")

pint$fofn <- sqrt(pint$n / (1 + 16 / 81 * pint$n))

# Figure 6
ggplot(data = subset(pint, N %in% seq(from=10, to=50, by=4)), aes(x = fofn, y = zscore, group = N)) +
  geom_point(aes(colour = N)) +
  geom_line(aes(colour = N)) +
  geom_line(y = cutoff, colour = "grey25") +
  xlab(expression(f(n))) +
  ylab("z-score") +
  scale_color_distiller(type = "seq", direction = -1, palette = "Greys") +
  theme_bw()

ggsave("figure6.svg", width=13, height=10, units="cm")
