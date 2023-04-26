# PGScpt
A Likelihood-Based Method for Risk Parameter Estimation under Polygenic Models from Case-Parent Trios

## Example Analysis
We provide a simple example for running our proposed method using simulated data. The R function of the proposed method is in [PGScpt](R/PGScpt.R). The R function to simulate data is available here [simulation](R/simulation.R).
```
rm(list = ls())
source("./R/PGScpt.R")
source("./R/simulation.R")
```
First simulate 100000 families based on a population disease model with two independent environmental variables $E_1$ (binary) and $E_2$ (continuous)
```
set.seed(04262023)
dat = sim_prospective_population(
  n_fam=100000, #Number of families in the population
  cor_e_prs=FALSE, #Assuming no correlation between PGS and E
  cor_strat=0.25, 
  rho2=0.25,
  rho_mf=0, #no assortative mating
  alpha_fam=-5.5, #Intercept term
  betaG_normPRS=0.4,
  betaE_bin=0.2, 
  betaE_norm=-0.6,
  betaGE_normPRS_bin=0.2, 
  betaGE_normPRS_norm=-0.4,
  envir=TRUE
)
```
View the simulated data.
```
#Convert the list to a data frame

```
Fit our proposed method.
```
startTime <- Sys.time()


endTime <- Sys.time()
```
Print running time.
```
print(endTime - startTime)

```
Output the final results.


# Questions
Please send your questions/suggestions to zwang389@jhu.edu
