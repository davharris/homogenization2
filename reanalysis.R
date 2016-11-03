library(tidyverse)
library(readxl)

data = read_excel("Appendix A3--Initial and current species lists.xlsx", sheet = 2)
data[is.na(data)] = 0

names = colnames(data)

data = data %>% 
  group_by(`Scientific name`) %>% 
  by_slice(colSums, .collate = "cols", .labels = TRUE)
names(data) = names
data[-1][data[-1] > 1] = 1


initial = select(data, matches("ini?tially"))
final = select(data, matches("currently"))

library(blender)


shuffle = function(x) {
  t(apply(x, 1, sample))
}


initial_jbar = jbar(initial)
shuffled_diffs = replicate(
  100,
  jbar(shuffle(final)) - jbar(shuffle(initial))
)
round(quantile(shuffled_diffs, c(.025, .975)), 3)
mean(jbar(final) - jbar(initial) - shuffled_diffs)


w = numeric(nrow(final))
x = numeric(nrow(final))
y = numeric(nrow(final))
z = numeric(nrow(final))
for (i in 1:nrow(final)) {
  temp = initial
  temp[i, ] = final[i, ]
  w[i] = sum(final[i, ] - initial[i, ])
  x[i] = jstar(temp) - jstar(initial)
  y[i] = jbar(temp) - initial_jbar
  z[i] = pstar(temp)
}


plot(y ~ I(w/20), subset = rowSums(initial) == 0 | rowSums(final) == 0)
lines(
  x[order(w)] ~ I(sort(w)/20), 
  subset = (rowSums(initial) == 0 | rowSums(final) == 0)[order(w)])

# % of variance explained
1 - Metrics::mse(y, x) / Metrics::mse(y, mean(y))

plot(y ~ z)
lines(x[order(z)] ~ sort(z))
abline(v = pstar(initial))


S = crossprod(as.matrix(initial))
S = S[row(S) > col(S)]

T = nrow(initial) - crossprod(1 - as.matrix(initial))
T = T[row(T) > col(T)]

N = choose(ncol(initial), 2)
correction = cov(T, S/T) / mean(T) * (N - 1) / N

jstar(initial) - jbar(initial)

jstar(initial) - correction
jbar(initial)


jbar(initial) + correction
jstar(initial)



a = c(0, 1)
b = c(1, 0)


x = cbind(replicate(1E2, a), replicate(1E2, b))
jbar(x)
jstar(x)


S = crossprod(as.matrix(x))
S = S[row(S) > col(S)]

T = nrow(x) - crossprod(1 - as.matrix(x))
T = T[row(T) > col(T)]

N = choose(ncol(x), 2)
correction = cov(T, S/T) / mean(T) * (N - 1) / N



jbar(x) + correction
jstar(x)

