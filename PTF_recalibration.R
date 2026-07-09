# Bayesian recalibration of Alexander's (1980) PTF using national soil data
#
# Author: Gabor Szatmari
# E-mail: szatmari.gabor@atk.hun-ren.hu


# 1. Packages ----
library(BayesianTools)
library(dplyr)
library(ggplot2)
library(hexbin)


# 2. Read national soil data ----
load("data.RData")


# 3. Bayesian recalibration ----

## 3.1. Likelihood ----
likelihood <- function(params){
  
  b0 <- params[1]
  b1 <- params[2]
  s0 <- params[3]
  s1 <- params[4]
  
  sigma <- s0 + s1 * oc # heteroscedastic error model
  
  if(any(sigma <= 0)) return(-Inf)
  
  pred <- b0 + b1 * sqrt(pmax(oc, 0)) # Alexander's (1980) PTF
  
  if(any(!is.finite(pred))) return(-Inf)
  
  ll <- sum(dnorm(bd_obs, mean = pred, sd = sigma, log = TRUE))
  
  if(!is.finite(ll)) return(-Inf)
  
  return(ll)
}


## 3.2. Priors ----
prior <- createUniformPrior(lower = c(1.0, -1.0, 0.01, -0.1),
                            upper = c(2.5,  0.0,  0.5,  0.1))


## 3.3. Setup ----
setup <- createBayesianSetup(likelihood = likelihood,
                             prior = prior,
                             names = c("b0","b1","s0", "s1"))


## 3.4. Run MCMC ----
set.seed(1234)
out <- runMCMC(setup,
               sampler = "DEzs",
               settings = list(iterations = 50000, nrChains = 3))

saveRDS(out, "New_PTF.rds")


## 3.5. Check diagnostics ----
plot(out)
summary(out)


# 4. Graphical presentation of the results ----
samples <- getSample(out, start = 5000)

oc_seq <- seq(min(oc), max(oc), length.out = 200)

pred <- lapply(oc_seq, function(x){ # propagation of posterior samples
  
  mu <- samples[, "b0"] + samples[, "b1"] * sqrt(x) # mean
  
  sigma <- samples[, "s0"] + samples[, "s1"] * sqrt(x) # standard deviation
  
  yrep <- rnorm(n = length(mu), mean = mu, sd = sigma) # BD posterior distribution
  
  data.frame(oc = x,
             mean = mean(mu),
             lower95 = quantile(yrep, 0.025),
             upper95 = quantile(yrep, 0.975),
             lower90 = quantile(yrep, 0.05),
             upper90 = quantile(yrep, 0.95))
})

pred <- bind_rows(pred)

pred$bd_original <- 1.72 - 0.294 * sqrt(pred$oc) # Applying Alexander's (1980) original PTF

p1 <- ggplot() +
  geom_hex(data = cbind(oc, bd_obs), aes(x = oc, y = bd_obs), bins = 45) +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  geom_ribbon(data = pred, aes(x = oc, ymin = lower90, ymax = upper90), alpha = 0.25) +
  geom_line(data = pred, aes(oc, mean, colour = "Recalibrated PTF"), linewidth = 1.3) +
  geom_line(data = pred, aes(oc, bd_original, colour = "Original PTF"), linewidth = 1.3) +
  scale_colour_manual(values = c("Original PTF" = "red", "Recalibrated PTF" = "blue"), name = NULL) +
  labs(x = "Organic carbon (%)", y = expression(paste("Bulk density (g ", cm^{-3}, ")"))) +
  theme_bw(base_size = 14)

p1 # plot the figure

ggsave("Fig01.jpg", # save the figure
       p1,
       width = 10,
       height = 6.5,
       dpi = 600)
