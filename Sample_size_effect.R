# Evaluation of the effect of sample size on recalibration
#
# Author: Gabor Szatmari
# E-mail: szatmari.gabor@atk.hun-ren.hu


# 1. Packages ----
library(BayesianTools)
library(caret)
library(dplyr)
library(ggplot2)
library(patchwork)


# 2. Read and pre-process data ----
load("data.RData")
data <- data.frame(oc = oc, bd = bd_obs)


# 3. Evaluation of the effect of sample size ----
set.seed(123)
idx <- createDataPartition(y = data$bd,
                           p = 0.8,
                           times = 10) # 80% recalibration and 20% validation

n_vals <- c(20, 50, 100, 200, 500, 1000, 2000, 5000) # recalibration sample sizes
n_rep1  <- 10 # repetition of data splitting
n_rep2 <- 10 # repetition of drawing recalibration samples

results <- data.frame() # container

for(rep1 in 1:n_rep1){
  
  train_full <- data[idx[[rep1]], ] # (full) recalibration set
  test <- data[-idx[[rep1]], ] # validation set
  
  for(n in n_vals){
    for(rep2 in 1:n_rep2){
      
      cat("Running n =", n, "rep_1 =", rep1, "rep_2 =", rep2, "\n")
      
      ## Recalibration sample
      idx_sub <- sample(1:nrow(train_full), n)
      train <- train_full[idx_sub, ] # recalibration set
      
      ## Likelihood
      likelihood <- function(params){
        
        b0 <- params[1]
        b1 <- params[2]
        s0 <- params[3]
        s1 <- params[4]
        
        sigma <- s0 + s1 * train$oc
        
        if(any(sigma <= 0)) return(-Inf)
        
        pred <- b0 + b1 * sqrt(pmax(train$oc, 0))
        
        if(any(!is.finite(pred))) return(-Inf)
        
        ll <- sum(dnorm(train$bd, mean = pred, sd = sigma, log = TRUE))
        
        if(!is.finite(ll)) return(-Inf)
        
        return(ll)
      }
      
      ## Priors
      prior <- createUniformPrior(lower = c(1.0, -1.0, 0.01, -0.1),
                                  upper = c(2.5,  0.0, 0.5,   0.1))
      
      ## Setup
      setup <- createBayesianSetup(likelihood = likelihood,
                                   prior = prior,
                                   names = c("b0","b1","s0", "s1"))
      
      ## Run MCMC
      out <- runMCMC(setup,
                     sampler = "DEzs",
                     settings = list(iterations = 8000, nrChains = 3))
      
      samples <- getSample(out, coda = FALSE)
      
      b0_ci <- quantile(samples[,"b0"], 0.975) - quantile(samples[,"b0"], 0.025)
      b1_ci <- quantile(samples[,"b1"], 0.975) - quantile(samples[,"b1"], 0.025)
      s0_ci <- quantile(samples[, "s0"], 0.975) - quantile(samples[, "s0"], 0.025)
      s1_ci <- quantile(samples[, "s1"], 0.975) - quantile(samples[, "s1"], 0.025)
      
      ## Prediction
      pred_matrix <- apply(samples, 1, function(p){
        
        mu <- p["b0"] + p["b1"] * sqrt(pmax(test$oc, 0))
        
        sigma <- p["s0"] + p["s1"] * test$oc
        
        rnorm(length(mu), mean = mu, sd = sigma)
      })
      
      pred_mean <- rowMeans(pred_matrix)
      
      ## Compute RMSE and coverage of 95% CI
      rmse <- sqrt(mean((pred_mean - test$bd)^2, na.rm=TRUE))
      lower <- apply(pred_matrix, 1, quantile, 0.025, na.rm=TRUE)
      upper <- apply(pred_matrix, 1, quantile, 0.975, na.rm=TRUE)
      coverage <- mean(test$bd >= lower & test$bd <= upper, na.rm=TRUE)
      
      ## Save the results
      results <- rbind(results, data.frame(n = n,
                                           rep1 = rep1,
                                           rep2 = rep2,
                                           rmse = rmse,
                                           coverage = coverage,
                                           b0_ci = b0_ci,
                                           b1_ci = b1_ci,
                                           s0_ci = s0_ci,
                                           s1_ci = s1_ci))
    }
  }
}

saveRDS(results, "Results_df.rds")

results_mean <- aggregate(.~n, results, mean) # compute mean
results_sd <- aggregate(.~n, results, sd) # compute standard deviation


# 4. Graphical presentation of the results ----
p1 <- ggplot(results_mean, aes(x = n, y = rmse)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = rmse - results_sd$rmse, ymax = rmse + results_sd$rmse), width = 0.05) +
  scale_x_log10() +
  labs(x = "Sample size", y = "RMSE", title = "(A)") +
  theme_bw()

p2 <- ggplot(results_mean, aes(x = n, y = coverage)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = coverage - results_sd$coverage, ymax = coverage + results_sd$coverage), width = 0.05) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  scale_x_log10() +
  labs(x = "Sample size", y = "Coverage", title = "(B)" ) +
  theme_bw()

p3 <- ggplot(results_mean,  aes(x = n, y = b1_ci)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = b1_ci - results_sd$b1_ci, ymax = b1_ci + results_sd$b1_ci), width = 0.05) +
  scale_x_log10() +
  labs(x = "Sample size", y = "95% CI width of b1", title = "(C)") +
  theme_bw()

p4 <- ggplot(results_mean, aes(x = n, y = s1_ci)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = s1_ci - results_sd$s1_ci, ymax = s1_ci + results_sd$s1_ci), width = 0.05) +
  scale_x_log10() +
  labs(x = "Sample size", y = "95% CI width of s1", title = "(D)") +
  theme_bw()

final_plot <- (p1 | p2) / (p3 | p4)

final_plot # plot the figure

ggsave("Fig02.jpg", # save the figure
       final_plot,
       width = 10,
       height = 8,
       dpi = 600)
