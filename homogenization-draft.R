knitr::opts_chunk$set(echo = FALSE)
library(tidyverse)
library(blender)
library(viridis)
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

cap1 = "Effect of an exotic species on $J^*$, given a landscape with 120 sites and 20 native species whose occupancy rates produce $p^{*} = 2/3$. The exotic species' effects on $J^*$ occur in three phases as it spreads across the landscape. When an exotic species begins spreading, it will decrease mean similarity until its occupancy rate reaches $p^{*}/2$ (blue arrow). As exotic occupancy continues to increase toward $p^*$, and $J^*$ will return to its initial value (purple arrow). Only once the species becomes more widespread than its native counterparts will it be able to raise $J^*$ above its initial value (red arrow)."

cap2 = "The effect of covariance on $\\bar{J}$ (Equation \\ref{decomposition}), holding each species' occurrence rate constant at 50% and the occupancy component constant at 0.33. In **A**, there are only two kinds of communities: light blue communities that contain species 1-50 and dark blue communities that contain species 51-100. Here, the covariance between $T_{ij}$ and $S_{ij} / T_{ij}$ is strongly negative and covariance contributes substantially to $\\bar{J}$. In **B**, the two community types are not as distinct, and the covariance effect is much smaller.  In **C**, there are no distinct community types and the covariance is not distinguishable from zero. Counter-intuitively, randomly smearing species across the landscape causes anti-homogenization, in terms of mean similarity."

# Figure 1 -------------------------------

invisible(cap1)
baseline = 2/3
K = 20
N = 120

# x and y values for the graph
values = seq(0, 1, length = N + 1)
scoop = sapply(values, function(x){jstar(c(rep(baseline, K), x), N)})

{
  plot(values, scoop, type = "n", xaxs = "i", bty = "l",
       ylab = expression(paste(J, "*")), axes = FALSE, 
       xlab = "Focal species occupancy", xlim = c(0, 1.01),
       ylim = c(min(scoop) - 0.005, max(scoop) + 0.005))
  
  axis(1, c(0, .333, .667, 1))
  axis(2, seq(0, 1, .01), las = 1)
  
  segments(baseline, 0, baseline, 1, col = "#00000050")
  segments(baseline/2, 0, baseline/2, 1, col = "#00000050")
  abline(h = scoop[1], col = "#00000050")
 
  text(baseline, max(scoop) - .005, expression(paste(p,"*")), pos = 4)
  text(baseline/2, max(scoop) - .005, expression(paste(p,"*/2")),
       pos = 4)
  
  is_blue = values < baseline/2
  is_red = values > baseline
  is_purple = values > baseline/2 & values < baseline
  
  lines(scoop ~ values, subset = is_blue, col = "blue", lwd = 2)
  lines(scoop ~ values, subset = is_purple, col = "purple", lwd = 2)
  lines(scoop ~ values, subset = is_red, col = "red", lwd = 2)
  
  # Differentiating invader
  i = sum(is_blue)
  arrows(values[i], scoop[i], values[i + 1], scoop[i + 1], col = "blue", length = .125, lwd = 2)
  
  # Un-differentiating invader
  i = match(TRUE, is_red) - 2
  arrows(values[i], scoop[i], values[i + 1], scoop[i + 1], col = "purple", length = .125, lwd = 2)
  
  # Homogenizing invader
  i = N
  arrows(values[i], scoop[i], values[i + 1], scoop[i + 1], col = "red", length = .125, lwd = 2)
}
# Figure 2 -------------------------------

invisible(cap2)
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
        strip.background = element_rect(fill = "white"), text = element_text(size = 11))
