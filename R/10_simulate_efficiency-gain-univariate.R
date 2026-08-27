# AR(1)-GARCH(1,1) Efficiency Gain: Grid Simulation
# Produces a heatmap of empirical variance ratio over (alpha, beta) parameter space

library(rugarch)
library(ggplot2)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# ─── Fixed parameters ─────────────────────────────────────────────────────────
set.seed(42)
n     <- 3000   # sample size
nsim  <- 500    # replications per grid point
phi   <- 0.5    # AR(1) coefficient
mu    <- 0      # mean
omega <- 0.5    # GARCH intercept

# ─── Parameter grid ───────────────────────────────────────────────────────────
alpha_grid <- seq(0.025, 0.90, by = 0.025)
beta_grid  <- seq(0.025, 0.90, by = 0.025)

grid <- expand.grid(alpha = alpha_grid, beta = beta_grid) %>%
  filter(alpha + beta < 0.98)

cat(sprintf("Grid has %d valid parameter combinations\n", nrow(grid)))

# ─── Simulation function ───────────────────────────────────────────────────────
simulate_ar1_garch11 <- function(n, phi, omega, alpha, beta, mu = 0) {
  sigma2    <- numeric(n)
  u         <- numeric(n)
  y         <- numeric(n)
  sigma2[1] <- omega / (1 - alpha - beta)
  u[1]      <- rnorm(1, 0, sqrt(sigma2[1]))
  y[1]      <- mu + u[1]
  for (t in 2:n) {
    sigma2[t] <- omega + alpha * u[t-1]^2 + beta * sigma2[t-1]
    u[t]      <- rnorm(1, 0, sqrt(sigma2[t]))
    y[t]      <- mu + phi * y[t-1] + u[t]
  }
  return(y)
}

# ─── GARCH specification ──────────────────────────────────────────────────────
spec_correct <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(1, 0), include.mean = TRUE)
)

# ─── Run one grid point ───────────────────────────────────────────────────────
run_grid_point <- function(alpha, beta, n, phi, omega, nsim, spec_correct) {
  phi_correct <- numeric(nsim)
  phi_arma    <- numeric(nsim)
  
  for (i in 1:nsim) {
    y <- simulate_ar1_garch11(n, phi, omega, alpha, beta)
    
    fit_correct <- tryCatch(
      ugarchfit(spec_correct, data = y, solver = "hybrid"),
      error = function(e) NULL
    )
    fit_arma <- tryCatch(
      arima(y, order = c(1, 0, 0)),
      error = function(e) NULL
    )
    
    if (!is.null(fit_correct) && !is.null(fit_arma)) {
      phi_correct[i] <- coef(fit_correct)["ar1"]
      phi_arma[i]    <- coef(fit_arma)["ar1"]
    } else {
      phi_correct[i] <- NA
      phi_arma[i]    <- NA
    }
  }
  
  phi_correct <- na.omit(phi_correct)
  phi_arma    <- na.omit(phi_arma)
  
  data.frame(
    alpha      = alpha,
    beta       = beta,
    rho        = alpha + beta,
    var_ratio  = var(phi_arma) / var(phi_correct),
    sd_correct = sd(phi_correct),
    sd_arma    = sd(phi_arma),
    n_ok       = length(phi_correct)
  )
}

# ─── Set up parallel backend ──────────────────────────────────────────────────
n_cores <- max(1, detectCores() - 1)
cat(sprintf("Using %d cores\n", n_cores))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export required objects and functions to workers
clusterExport(cl, c("simulate_ar1_garch11", "n", "phi", "omega",
                    "nsim", "spec_correct"))
clusterEvalQ(cl, {
  library(rugarch)
  library(dplyr)
})

# ─── Run over grid in parallel ────────────────────────────────────────────────
cat("Running grid simulation in parallel...\n")
t_start <- Sys.time()

results_list <- foreach(
  i        = seq_len(nrow(grid)),
  .combine = rbind,
  .packages = c("rugarch")
) %dopar% {
  run_grid_point(
    alpha        = grid$alpha[i],
    beta         = grid$beta[i],
    n            = n,
    phi          = phi,
    omega        = omega,
    nsim         = nsim,
    spec_correct = spec_correct
  )
}

stopCluster(cl)

t_end <- Sys.time()
cat(sprintf("Done in %.1f minutes\n", as.numeric(t_end - t_start, units = "mins")))

results_df <- as.data.frame(results_list)

# ─── Figure ───────────────────────────────────────────────────────────────────
p <- ggplot(results_df, aes(x = alpha, y = beta, fill = var_ratio)) +
  geom_tile() +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0.1, 0.9, by = 0.1)) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    name      = "Variance ratio",
    limits    = c(1, NA)
  ) +
  geom_abline(
    slope     = -1,
    intercept = seq(0.3, 0.9, by = 0.2),
    linetype  = "dashed",
    colour    = "grey40",
    linewidth = 0.4
  ) +
  annotate("text", x = 0.42, y = 0.52, label = "rho == 0.9",
           parse = TRUE, colour = "grey40", size = 3) +
  annotate("text", x = 0.42, y = 0.32, label = "rho == 0.7",
           parse = TRUE, colour = "grey40", size = 3) +
  annotate("text", x = 0.42, y = 0.12, label = "rho == 0.5",
           parse = TRUE, colour = "grey40", size = 3) +
  labs(
    x = expression(alpha ~ "(ARCH coefficient)"),
    y = expression(beta ~ "(GARCH coefficient)")
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave("Figures/efficiency_gain_heatmap.pdf", p, width = 7, height = 5.5)
cat("Figure saved.\n")
print(p)
