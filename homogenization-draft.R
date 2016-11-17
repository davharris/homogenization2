knitr::opts_chunk$set(echo = FALSE)
library(tidyverse)
library(blender)
library(viridis)
library(ggrepel)
data(PLANTS)
# Shuffling error -------------------------------

# Quantify the errors introduced by shuffling species' locations
# in the 47 landscapes described by the USDA PLANTS data set.
abbreviations = strsplit(names(PLANTS), " ") %>% 
  map(1) %>% 
  unlist() %>% 
  unique()

# Save a vector of mean absolute errors, one for each landscape.
# The value reported in the manuscript is the mean of this vector.
errors = sapply(
  abbreviations, 
  function(ab){
    x = bind_rows(PLANTS[grepl(ab, names(PLANTS))])
    shuffled_values = replicate(10, jbar(t(apply(x, 1, sample))))
    mean(abs(shuffled_values - jbar(x)))
  }
)

# Figure captions -------------------------------
datacap = "\\label{fig:datascoop} 
**A.** 
Across data sets from 47 US states, mean similarity is almost entirely determined by species' occupancy rates (as calculated with the $J^*$ approximation in Equation \ref{jstar}).
**B.**
The proportion of sites occupied by different exotic species in Louisiana explains almost all of the variance in their effects on mean similarity [@harris_occupancy_2011].
Louisiana was chosen because it had the median $R^2$ value of the 47 available landscapes. 
As a species invades and spreads across the landscape, it passes through three phases, indicated by colors.
The boundaries between these phases are given by $p^*/2$ and $p^*$, described in Equation \\ref{pstar}.
*Phase I*: the invader is rare, and increasing its occupancy magnifies its differentiating influence.
*Phase II*: spreading the invader makes it more similar to the background of native species, so its net effect approaches zero as the invader spreads.
*Phase III*: the invader is more common than the background of native species, and incresaes mean similarity as it spreads.
"

covcap = "\\label{fig:covariance} 
The effect of covariance on $\\bar{J}$ (Equation \\ref{decomposition}), holding each species' occurrence rate constant at 50% and the occupancy component constant at 0.33. 
In **A**, there are only two kinds of communities: light blue communities that contain species 1-50 and dark blue communities that contain species 51-100. 
Here, the covariance between $T_{ij}$ and $S_{ij} / T_{ij}$ is strongly negative and covariance contributes substantially to $\\bar{J}$. 
In **B**, the two community types are not as distinct, and the covariance effect is much smaller.
In **C**, there are no distinct community types and the covariance is not distinguishable from zero. Counter-intuitively, randomly smearing species across the landscape causes anti-homogenization, in terms of mean similarity."

# Figure 1---------------------------------
all = blend(PLANTS)
par(mfrow = c(2, 1))
par(mar = c(5, 4 + 2, 4, 2) + 0.1)
par(bty = "l")
state = all$LA

brewer_cols = RColorBrewer::brewer.pal(3, "PuOr")
colors = as.character(
  cut(state$species.delta.table$occupancy, 
      c(-Inf, state$p.Star/2, state$p.Star, Inf),
      labels = c(brewer_cols[1], "black", brewer_cols[3]))
)

states = ggplot(all$summary, 
                aes(x = J.Star, y = J.Bar, label = rownames(all$summary))) + 
  geom_point() + 
  geom_text_repel(segment.alpha = 0.25) + 
  cowplot::theme_cowplot() + 
  geom_abline(intercept = 0, slope = 1, col = alpha("black", 1/3)) + 
  xlab("Predicted mean similarity (J*)") + 
  ylab(expression(paste("Observed mean similarity (", bar(J), ")"))) +
  theme(plot.margin = unit(c(1, 1, 1.5, 1), "lines"))


LA = ggplot() + 
  geom_point(data = state$species.delta.table, aes(x = occupancy, y = delta.J.Bars), 
             color = colors) + 
  geom_path(data = state$scoop, aes(x = x, y = y), size = 1, alpha = 0.75) + 
  cowplot::theme_cowplot() +
  geom_vline(xintercept = state$p.Star / c(1,2), alpha = 0.5) + 
  geom_hline(yintercept = 0, alpha = 0.5) + 
  xlim(0, 1.01) + 
  coord_cartesian(expand = FALSE) + 
  geom_text(aes(x = state$p.Star / c(1,2) + 0.01, 
                y = max(state$species.delta.table$delta.J.Bars - 5E-5), 
                label = c("p*", "p*/2")),
            hjust = 0, fontface = "bold", size = 5) +
  expand_limits(
    y = c(
      min(state$species.delta.table$delta.J.Bars - 1E-4),
      max(state$species.delta.table$delta.J.Bars + 1E-4))
  ) +
  xlab("Proportion of sites occupied by invader") +
  ylab("Effect on mean similarity") +
  theme(plot.margin = unit(c(1.5, 1, 1, 1), "lines"))

cowplot::plot_grid(states, LA, nrow = 2, align = "v",
                   labels = c("A", "B"), vjust = 1.25)
# Figure 2 -------------------------------

invisible(covcap)
N = 80
K = 100

pa1 = matrix(NA, K, N)
for (i in 1:K) {
  if (i <= K/2) {
    pa1[i, ] = 1:N <= N/2
  } else {
    pa1[i, ] = 1:N > N/2
  }
}
pa2 = t(apply(pa1, 1, sample, prob = 1/seq_len(N)^2))
pa3 = t(apply(pa1, 1, sample))

pa = list(pa1, pa2, pa3)
labels = sapply(
  1:3,
  function(i){
    S = as.numeric(as.dist(crossprod(pa[[i]])))
    T = as.numeric(as.dist(K - crossprod(1 - pa[[i]])))
    
    paste0(LETTERS[i], ": mean similarity = ",
           format(jbar(pa[[i]]), nsmall = 2, digits = 2),
           "\noccupancy component: ",
           format(jstar(pa[[i]]), nsmall = 2, digits = 2),
           "\ncovariance component: ", 
           format(round(-cov(T, S/T) / mean(T), 2),
                  scientific = FALSE, nsmall = 2))
  })

plots = lapply(
  1:3,
  function(i){
    x = pa[[i]]
    
    intensity = colSums(1E-6 + x[seq_len(K/2), ]) / K * 2
    x = t(t(x) * intensity)
    x = x[ , order(intensity)]
    
    out = reshape2::melt(x * 1)
    colnames(out) = c("species", "site", "value")
    out[out == 0] = NA
    cbind(out, panel = labels[[i]], stringsAsFactors = FALSE)
  }
) %>% bind_rows()

ggplot(plots, aes(x = species, y = site, fill = value)) + 
  geom_tile() + 
  facet_grid(~panel) +
  scale_fill_gradient(na.value = "gray97", limits = c(0, 1 + 1E-6)) +
  guides(fill = FALSE) + 
  coord_cartesian(c(0, K), c(0, N), FALSE) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"), panel.spacing.x = unit(1, "lines"),
        strip.background = element_rect(fill = "white"), text = element_text(size = 11), axis.text = element_text(size = 9))
