## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(AlphaSimR)
library(FieldSimR)
library(ggplot2)

## -----------------------------------------------------------------------------
n_envs <- 3 # Number of environments.
n_traits <- 2 # Number of traits.
n_blocks <- c(2, 3, 2) # Number of blocks in the three environments, respectively.
block_dir <- "col" # Layout of blocks (“side-by-side”).
n_cols <- c(10, 15, 10) # Total umber of columns per location.
n_rows <- 20 # Total number of rows per location.
plot_length <- 8 # Plot length; here in meters (column direction).
plot_width <- 2 # Plot width; here in meters (row direction).

## -----------------------------------------------------------------------------
H2 <- c(0.3, 0.3, 0.3, 0.5, 0.5, 0.5) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

## -----------------------------------------------------------------------------
var <- c(0.085, 0.12, 0.06, 15.1, 8.5, 11.7) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

## -----------------------------------------------------------------------------
# Calculation of error variances based on the genetic variance and target heritability vectors.
calc_var_R <- function(var, H2) {
  varR <- (var / H2) - var
  return(varR)
}

var_R <- calc_var_R(var, H2)
var_R # Vector of error variances: c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

## -----------------------------------------------------------------------------
spatial_model <- "Bivariate" # Spatial error model.
prop_spatial <- 0.4 # Proportion of spatial to total error variance.

S_cor_R <- rand_cor_mat(n_traits, min_cor = 0, max_cor = 0.5, pos_def = TRUE)
S_cor_R

## -----------------------------------------------------------------------------
prop_ext <- 0.2
ext_dir <- "row"
ext_ord <- "zig-zag"

E_cor_R <- rand_cor_mat(n_traits, min_cor = 0, max_cor = 0.5, pos_def = TRUE)
E_cor_R

## -----------------------------------------------------------------------------
error_df <- field_trial_error(
  n_envs = n_envs,
  n_traits = n_traits,
  n_blocks = n_blocks,
  block_dir = block_dir,
  n_cols = n_cols,
  n_rows = n_rows,
  plot_length = plot_length,
  plot_width = plot_width,
  var_R = var_R,
  R_cor_R = NULL,
  spatial_model = spatial_model,
  prop_spatial = prop_spatial,
  S_cor_R = S_cor_R,
  prop_ext = prop_ext,
  ext_dir = ext_dir,
  ext_ord = ext_ord,
  E_cor_R = E_cor_R,
  return_effects = TRUE
)

## -----------------------------------------------------------------------------
error_env2 <- error_df$plot_df[error_df$plot_df$env == 2, ]
e_comp_env2 <- error_df$Trait.1[error_df$Trait.1$env == 2, ]

