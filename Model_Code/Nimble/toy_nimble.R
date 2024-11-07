#==============================
# NIMBLE toy example 
# Date updated: 2024/10/23
#=============================

# Load necessary packages
library(nimble)
library(coda)
library(mcmcplots)
library(sfsmisc) #mult.fig

wd <- "~/Documents/Github/WRPC/"
res_dir <- "Code/NIMBLE/Output/"
source(paste0(wd, "Code/simulate_data_wrpc.R"))



#================ (1) Simulate data ============================================

# Population size and strata dimensions
N <- 1000; H <- 2; S <- 2; J <- 10; R <- 4 
N_h <- c(500, 500)
N_s <- c(650, 350)

# Generate C ~ S
K <- 3  
formula_c <- "~ s_all"
V_unique <- data.frame(s_all = as.factor(1:S))
# Corresponds to an overall true_pi ~= (0.253, 0.522, 0.225)
# Matrix of global class assignment probabilities for each level of s_i
pi_global_mat <- matrix(c(0.3, 0.5, 0.2,   # global class probs for s_i=1
                          0.1, 0.6, 0.3),  # global class probs for s_i=2
                        byrow = TRUE, nrow = S, ncol = K)
beta_mat_c <- get_betas_c(pi_mat = pi_global_mat, formula_c = formula_c, 
                          V_unique = V_unique)

# Generate H 
h_all <- c(rep(1, times = N_h[1]), rep(2, times = N_h[2]))

# Global-local allocation probability of global assignment
nu <- matrix(c(0.3, 0.7), byrow = FALSE, nrow = H, ncol = J)

# Generate X ~ C
formula_x_global <- "~ c_all"
formula_x_local <- "~ h_all"
modal_theta_prob <- 0.85
V_unique <- expand.grid(c_all = as.factor(1:K), h_all = as.factor(1:H))
# Matrix of global profiles
profiles_global <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J), 
                                               rep(3, times = 0.5 * J)),
                                        C2 = c(rep(4, times = 0.2 * J), 
                                               rep(2, times = 0.8 * J)),
                                        C3 = c(rep(3, times = 0.3 * J), 
                                               rep(4, times = 0.4 * J),
                                               rep(1, times = 0.3 * J))))
# List of beta matrices for each item j
beta_list_x_global <- get_betas_x(profiles = profiles_global, R = R, 
                                  modal_theta_prob = modal_theta_prob, 
                                  formula_x = formula_x_global, 
                                  V_unique = V_unique)
# List of matrices of local profiles for each subpopulation
profiles_local <- as.matrix(data.frame(H1 = rep(1, times = J), 
                                       H2 = rep(4, times = J)))
# List over subpopulation of list of matrices of beta matrices for each item j
beta_list_x_local <- get_betas_x(profiles = profiles_local, R = R, 
                                 modal_theta_prob = modal_theta_prob, 
                                 formula_x = formula_x_local, 
                                 V_unique = V_unique)

# Simulate population
pop_seed <- 1  # Set seed
sim_pop <- simulate_pop_wrpc_v2(N = N, H = H, S = S, J = J, R = R, K = K, 
                             N_h = N_h, N_s = N_s,
                             modal_theta_prob = modal_theta_prob,
                             formula_c = formula_c, 
                             formula_x_global = formula_x_global, 
                             formula_x_local = formula_x_local, 
                             beta_mat_c = beta_mat_c, 
                             beta_list_x_global = beta_list_x_global,
                             beta_list_x_local = beta_list_x_local,
                             pop_seed = pop_seed, save_res = TRUE, 
                             save_path = paste0(wd, "/Code/Nimble/"))

parameter_true <- list(
  c_all = sim_pop$true_Ci,
  g_mat = sim_pop$true_Gij,
  pi_global = as.vector(sim_pop$true_pi_global),
  theta_global = sim_pop$true_global_thetas,
  theta_local = sim_pop$true_local_thetas,
  nu = sim_pop$true_nu
)


#================ (2) Run RPC using NIMBLE =====================================


## Model Code
code <- nimbleCode({
  # Priors
  pi_global[1:K] ~ ddirch(alpha[1:K])
  for (j in 1:J) {
    for (k in 1:K) {
      theta_global[j, k, 1:R] ~ ddirch(eta[1:R])
    }
    for (h in 1:H) {
      theta_local[j, h, 1:R] ~ ddirch(eta[1:R])
      nu[h, j] ~ dbeta(1.0, 1.0)
    }
  }
  
  # MENGBING SAYS: DON'T SPECIFY POSTERIOR!!!!!!!!!!
  
  # Likelihood
  for (i in 1:N) {
    c_all[i] ~ dcat(pi_global[1:K])  # latent class assignment
    for (j in 1:J) {
      g_mat[i, j] ~ dbern(nu[h_all[i], j])  # global-local assignment
      # g_mat[i, j] ~ dbern((nu[h_all[i], j] * theta_global[j, c_all[i], x_mat[i, j]]) / 
      #                       (nu[h_all[i], j] * theta_global[j, c_all[i], x_mat[i, j]]) + 
      #                       ((1 - nu[h_all[i], j]) * theta_local[j, h_all[i], x_mat[i, j]]))
      theta_temp[i, j, 1:R] <- g_mat[i, j] * theta_global[j, c_all[i], 1:R] + 
        (1 - g_mat[i, j]) * theta_local[j, h_all[i], 1:R]
      x_mat[i, j] ~ dcat(theta_temp[i, j, 1:R])
      # if (g_mat[i, j] == 1) {
      #   x_mat[i, j] ~ dcat(theta_global[j, c_all[i], 1:R])
      # } else {
      #   x_mat[i, j] ~ dcat(theta_local[j, h_all[i], 1:R])
      # }
    }
    
    # for (k in 1:K) {
    #   # prod_k <- pi_global[k]
    #   for (j in 1:J) {
    #     theta_temp[i, k, j] <- theta_global[j, k, x_mat[i, j]] * g_mat[i, j]
    #     # if (g_mat[i, j] == 1) {
    #     #   prod_k <- prod_k * theta_global[j, k, x_mat[i, j]] 
    #     # }
    #   }
    #   # probs_c[k] <- prod_k
    #   
    #   probs_c[i, k] <- pi_global[k] * prod(theta_temp[i, k, ])
    # }
    # c_all[i] ~ dcat(probs_c[i, 1:K])
  }
  
  
})


# Model constants
consts <- list(N = N,
               J = J, 
               K = K, 
               R = R, 
               H = H,
               alpha = rep(K, K),
               eta = rep(1.0, R), 
               h_all = sim_pop$true_Hi)

# Data
data_list <- list(x_mat = as.matrix(sim_pop$X_data))

# Initializations
# init_list <- list(
#   c_all = sim_pop$true_Ci,
#   g_mat = sim_pop$true_Gij,
#   pi_global = as.vector(sim_pop$true_pi_global),
#   theta_global = sim_pop$true_global_thetas,
#   theta_local = sim_pop$true_local_thetas,
#   nu = sim_pop$true_nu,
#   theta_temp = array(0, dim = c(N, J, R))
# )
# Try random initializations
init_list <- list(
  c_all = rep(1:K, length.out = N),
  g_mat = matrix(rbinom(n = N, size = 1, prob = 0.5), nrow = N, ncol = J),
  pi_global = rep(1, times = K) / K,
  theta_global = array(1 / R, dim = c(J, K, R)),
  theta_local = array(1 / R, dim = c(J, H, R)),
  nu = matrix(0.5, nrow = H, ncol = J),
  theta_temp = array(0, dim = c(N, J, R))
)


# Model creation and compilation
model <- nimbleModel(code, data = data_list, inits = init_list, constants = consts)

# MCMC configuration, creation, and compilation
conf <- configureMCMC(model,
                      monitors = c("pi_global", "theta_global", "theta_local", 
                                   "nu", "c_all", "g_mat"))
mcmc_model <- buildMCMC(conf)
# Saves mcmc sampler configuration in file 'code_MID_8_nfCode'
compiled_model_initial <- compileNimble(model, showCompilerOutput = TRUE,
                                        dirName = paste0(wd, res_dir))

# save cmcmc for future use
compiled_model <- compileNimble(mcmc_model, project = compiled_model_initial)

end_time_compile <- Sys.time()

# start sampling
start_time_draw <- Sys.time()
samples <- runMCMC(compiled_model, niter = 1000, nburnin = 0, setSeed = 1, nchains = 1)
end_time_draw <- Sys.time()

start_time_draw <- Sys.time() # 6 mins
samples10000 <- runMCMC(compiled_model, niter = 10000, nburnin = 0, setSeed = 1, nchains = 1)
end_time_draw <- Sys.time()

# saveRDS(samples, file = paste0(wd, res_dir, "nimble_toy_randomInit.rds"))
# save(model = model, configured_model = conf, mcmc_model = mcmc_model,
#      compiled_model_initial = compiled_model_initial,
#      compiled_model = compiled_model,
#      samples = samples, file = paste0(wd, res_dir, "nimble_toy.RData"))




#================ (3) Examine NIMBLE output ====================================

## load nimble results
load(paste0(wd, res_dir, "nimble_toy.RData"))
samples <- readRDS(paste0(wd, res_dir, "nimble_toy_randomInit.rds"))

nchains <- length(samples)
# s <- as.mcmc.list(lapply(1:nchains, function(i) coda::mcmc(samples[[i]])))
s <- coda::mcmc(samples)  # columns are the parameters
nsamples <- nrow(s)
colors_plot <- mcmcplotsPalette(nchains)


## plot trace of theta_global
j <- 1
r <- 1
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Item probability: item ", j),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  traplot(s, parms = paste0("theta_global[", j, ", ", k, ", ", r, "]"), 
          main = paste0("theta_global[", j, ", ", k, ", ", r, "]"), 
          greek = FALSE, auto.layout = FALSE, ylim = c(0,1))
  abline(h = parameter_true$theta_global[j, k, r], col = "black")
}

## plot trace of theta_local
j <- 1
r <- 1
mult.fig(mfrow = c(H, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main= paste0("Item probability: item ", j),
         cex = 1, marP = - c(0, 1, 2, 0))
for (h in 1:H) {
  traplot(s, parms = paste0("theta_local[", j, ", ", h, ", ", r, "]"), 
          main = paste0("theta_local[", j, ", ", h, ", ", r, "]"), 
          greek = FALSE, auto.layout = FALSE, ylim = c(0, 1))
  abline(h = parameter_true$theta_global[j, h, r], col = "black")
}


## plot trace of pi
mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main = paste("Trace plots for pi"),
         cex = 1, marP = - c(0, 1, 2, 0))
for (k in 1:K) {
  traplot(s, parms = paste0("pi_global[", k, "]"), 
          main = paste0("pi_global[", k, "]"), 
          greek = FALSE, auto.layout = FALSE, ylim = c(0, 1))
  abline(h = parameter_true$pi_global[k], col = "black")
}


## plot trace of nu
mult.fig(mfrow = c(H, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main = paste("Trace plots for nu"),
         cex = 1, marP = - c(0, 1, 2, 0))
j=1
for (h in 1:H) {
  traplot(s, parms = paste0("nu[", h, ", ", j, "]"), 
          main = paste0("nu[", h, ", ", j, "]"), 
          greek = FALSE, auto.layout = FALSE, ylim = c(0, 1))
  abline(h = parameter_true$nu[h, j], col = "black")
}


## look at c_all --------
index <- c(1, 500, 1000)
mult.fig(mfrow = c(3, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main = paste("c_all assignment"),
         cex = 0.8, marP = - c(0, 1, 2, 0))
for (i in index) {
  traplot(s, parms = paste0("c_all[", i, "]"),
          main = paste0("c_all[", i, "]"),
          greek = FALSE, auto.layout = FALSE)
  abline(h = parameter_true$c_all[i], col = "black")
}
prop.table(table(s[, "c_all[1]"]))
prop.table(table(s[, "c_all[500]"]))
prop.table(table(s[, "c_all[1000]"]))


## look at g_mat --------
i <- 1
var_index <- c(1, 5, 10)
mult.fig(mfrow = c(3, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
         main = paste("g_mat assignment"),
         cex = 0.8, marP = - c(0, 1, 2, 0))
for (j in var_index) {
  traplot(s, parms = paste0("g_mat[", i, ", ", j, "]"),
          main = paste0("g_mat[", i, ", ", j, "]"),
          greek = FALSE, auto.layout = FALSE)
  abline(h = parameter_true$g_mat[i, j], col = "black")
}
prop.table(table(s[, "g_mat[1, 1]"]))
prop.table(table(s[, "g_mat[1, 5]"]))
prop.table(table(s[, "g_mat[1, 10]"]))

# i = 1; j = 1
# index <- c(1, 5, 10)
# mult.fig(mfrow = c(3, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
#          main = paste("c_all assignment"),
#          cex = 0.8, marP = - c(0, 1, 2, 0))
# for (i in index) {
#   traplot(s, parms = paste0("c_all[", i, "]"),
#           main = paste0("c_all[", i, "]"),
#           greek = FALSE, auto.layout = FALSE)
#   abline(h = parameter_true$c_all[i], col = "black")
#   matplot(cumsum(s[, paste0("c_all[", i, "]")]-1) / (1:nsamples), type="n",
#           xaxt="n", yaxt="n", bty="n", ylab = "", ylim = c(0, 1), 
#           main = paste0("Proportion of c_all in class 2, i = ", i))
#   .graypr()
#   matlines(cumsum(s[, paste0("c_all[", i, "]")]-1) / (1:nsamples), 
#            type = "l", col = colors_plot, lty = 1, ylim = c(0, 1))
#   abline(h = sum(parameter_true$c_all[i]==2) / J, col = "black", lty = 2)
# }
# table(out_mcmc[[n]]$Z_samples[i, j, ])
# parameter_true$G[i,]
# table(parameter_true$Z[i,])
# 
# 
# ## look at c_all --------
# i = 1; j = 1
# var_index <- c(1, 15, 30)
# mult.fig(mfrow = c(3, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
#          main = paste("c_all assignment"),
#          cex = 0.8, marP = - c(0, 1, 2, 0))
# for (j in var_index) {
#   matplot(cumsum(s[, paste0("c_all[", i, ", ", j, "]")]-1) / (1:max_iter), type="n",
#           xaxt="n", yaxt="n", bty="n", ylab = "", ylim = c(0, 1), 
#           main = paste0("Proportion of Z_ij in class 2, i = ", i, ", j = ", j))
#   .graypr()
#   matlines(cumsum(s[, paste0("c_all[", i, ", ", j, "]")]-1) / (1:max_iter), type = "l", col = colors_plot, lty = 1, ylim = c(0, 1))
#   abline(h = sum(parameter_true$c_all[i, ]==2) / J, col = "black", lty = 2)
# }
# table(out_mcmc[[n]]$Z_samples[i, j, ])
# parameter_true$G[i,]
# table(parameter_true$Z[i,])



# ## look at alpha_g
# var_index <- 1:K
# mult.fig(mfrow = c(K, 1), mar = c(2, 2, 1.5, 0.25) + 0.1,
#          main = paste("alpha"),
#          cex = 0.8, marP = - c(0, 1, 2, 0))
# for (k in 1:K) {
#   # traplot(s, parms = paste0("gamma[", k, "]"), 
#   #         main = paste0("gamma[", k, "]"), 
#   #         greek = FALSE, auto.layout = FALSE)
#   alpha_k_samples <- s[, paste0("gamma[", k, "]")] * s[,"xi"]
#   matplot(alpha_k_samples, type="n",
#           xaxt="n", yaxt="n", bty="n", ylab = "", ylim = c(0, 1),
#           main = paste0("gamma[", k, "]"))
#   .graypr()
#   matlines(alpha_k_samples, type = "l", col = colors_plot, lty = 1, ylim = c(0, 1))
#   abline(h = parameter_true$alpha_g[k], col = "black", lty = 2)
# }




