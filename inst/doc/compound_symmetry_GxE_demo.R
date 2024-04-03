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
ntraits <- 2 # Number of traits.
nenvs <- 3 # Number of environments.
nreps <- c(2, 2, 3) # Number of replicates of each genotype in environments 1, 2, and 3.


nind <- 20 # Number of founder genotypes in the population.
nchr <- 10 # Number of chromosomes.
nseg_sites <- 200 # Number of QTN per chromosome.

## -----------------------------------------------------------------------------
mean <- c(4.9, 5.4, 5.1, 235.2, 228.5, 239.1) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

## -----------------------------------------------------------------------------
var <- c(0.08, 13) # c(grain yield, plant height)

## -----------------------------------------------------------------------------
prop_main <- c(0.4, 0.6) # c(grain yield, plant height)

## -----------------------------------------------------------------------------
corA <- matrix( # Matrix of additive genetic correlations grain yield and plant height.
  c(
    1.0, 0.5,
    0.5, 1.0
  ),
  ncol = 2
)

## ----echo=FALSE---------------------------------------------------------------
corA

## -----------------------------------------------------------------------------
meanDD <- c(0.4, 0.4, 0.4, 0.1, 0.1, 0.1) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

varDD <- c(0.2, 0.2) # c(grain yield, plant height)

prop_mainDD <- 0.4 # Same value set for traits 1 and 2.

corDD <- diag(2)

## ----echo=FALSE---------------------------------------------------------------
corDD

## -----------------------------------------------------------------------------
input_asr <- compsym_asr_input(
  ntraits = ntraits,
  nenvs = nenvs,
  mean = mean,
  var = var,
  prop.main = prop_main,
  corA = corA,
  meanDD = meanDD,
  varDD = varDD,
  prop.mainDD = prop_mainDD,
  corDD = corDD
)

## -----------------------------------------------------------------------------
founders <- runMacs( # Simulation of founder genotypes using AlphaSimR's "MAIZE" presets
  nInd = nind, # to mimic the species' evolutionary history.
  nChr = nchr,
  segSites = nseg_sites,
  species = "MAIZE",
  nThreads = 2
)

SP <- SimParam$new(founders)

## -----------------------------------------------------------------------------
SP$addTraitAD( # Additive + dominance trait simulation.
  nQtlPerChr = nseg_sites,
  mean = input_asr$mean,
  var = input_asr$var,
  corA = input_asr$corA,
  meanDD = input_asr$meanDD,
  varDD = input_asr$varDD,
  corDD = input_asr$corDD,
  useVarA = FALSE
)

founders <- newPop(founders)

## -----------------------------------------------------------------------------
pool_A <- makeDH(founders[1:10], nDH = 1) # Pool A: 1 DH line from founders 1 to 10, respectively.
pool_B <- makeDH(founders[11:20], nDH = 1) # Pool B: 1 DH line from founders 11 to 20, respectively.

dh_lines <- mergePops(list(pool_A, pool_B))

factorial_plan <- as.matrix(expand.grid(A = pool_A@id, B = pool_B@id)) # Factorial crossing plan.

hybrid_pop <- makeCross(pop = dh_lines, crossPlan = factorial_plan, nProgeny = 1) # Hybrid genotypes.

## -----------------------------------------------------------------------------
gv_df <- compsym_asr_output(
  pop = hybrid_pop,
  ntraits = ntraits,
  nenvs = nenvs,
  nreps = nreps
)

## ----echo=FALSE---------------------------------------------------------------
head(gv_df)

## ----echo=TRUE, fig.height = 4, fig.width = 9, fig.align = "center"-----------
ggplot(gv_df, aes(x = gv.Trait1, fill = factor(env))) +
  geom_histogram(color = "#e9ecef", alpha = 0.8, position = "identity", bins = 50) +
  scale_fill_manual(values = c("violetred3", "goldenrod3", "skyblue2")) +
  labs(x = "Genetic values for grain yield (t/ha)", y = "Count", fill = "Environment")

