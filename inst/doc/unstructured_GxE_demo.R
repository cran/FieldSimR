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
nreps <- c(2, 2, 3) # Number of replicates tested within environments 1, 2 and 3.


nind <- 20 # Number of founder genotypes in the population.
nchr <- 10 # Number of chromosomes.
nseg_sites <- 200 # Number of QTN per chromosome.

## -----------------------------------------------------------------------------
mean <- c(4.9, 5.4, 5.1, 235.2, 228.5, 239.1) # c(Yld:E1, Yld:E2, Yld:E3, Prt:E1, Prt:E2, Prt:E3)

## -----------------------------------------------------------------------------
var <- c(0.085, 0.12, 0.06, 15.1, 8.5, 11.7) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)

## -----------------------------------------------------------------------------
meanDD <- c(0.4, 0.4, 0.4, 0.1, 0.1, 0.1) # c(Yld:E1, Yld:E2, Yld:E3, Pht:E1, Pht:E2, Pht:E3)
varDD <- 0.2 # Same value set for all environment-within-trait combinations

## -----------------------------------------------------------------------------
TcorA <- matrix( # Matrix of additive genetic correlations between the two traits.
  c(
    1.0, 0.6,
    0.6, 1.0
  ),
  ncol = 2
)

## ----echo=FALSE---------------------------------------------------------------
TcorA

## -----------------------------------------------------------------------------
EcorA <- matrix(
  c( # Matrix of additive genetic correlations between the three environments.
    1.0, 0.4, 0.6,
    0.4, 1.0, 0.5,
    0.6, 0.5, 1.0
  ),
  ncol = 3
)

## ----echo=FALSE---------------------------------------------------------------
EcorA

## ----message=FALSE------------------------------------------------------------
corA <- rand_cor_mat( # Additive genetic correlation structure.
  (ntraits * nenvs), # Could be used instead of TcorA and EcorA.
  min.cor = 0.1,
  max.cor = 0.9,
  pos.def = TRUE
)

round(corA, 2)

## -----------------------------------------------------------------------------
corDD <- diag(6)

## ----echo=FALSE---------------------------------------------------------------
corDD

## -----------------------------------------------------------------------------
input_asr <- unstr_asr_input(
  ntraits = ntraits,
  nenvs = nenvs,
  mean = mean,
  var = var,
  TcorA = TcorA,
  EcorA = EcorA,
  meanDD = meanDD,
  varDD = varDD,
  corDD = corDD
)

## -----------------------------------------------------------------------------
founders <- runMacs( # Simulation of founder genotypes using AlphaSimR's "MAIZE" presets
  nInd = nind, # to mimic the species' evolutionary history.
  nChr = nchr,
  segSites = nseg_sites,
  inbred = FALSE,
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
gv_df <- unstr_asr_output(
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

