# pgs.cpt
A Likelihood-Based Method for Risk Parameter Estimation under Polygenic Models from Case-Parent Trios

## Example Analysis
We provide a simple example for running our proposed method using simulated data. The R function of the proposed method is in [pgs.cpt](R/PGScpt.R). The R function to simulate data is available here [simulation](R/simulation.R).
```
rm(list = ls())
require(devtools)
source_url("https://github.com/ziqiaow/pgs.cpt/blob/main/R/PGScpt.R?raw=TRUE")
source_url("https://github.com/ziqiaow/pgs.cpt/blob/main/R/simulation.R?raw=TRUE")

#If directly downloaded from Github 
source("./R/PGScpt.R")
source("./R/simulation.R")
```
First simulate 200000 families based on a population disease risk model following a logistic regression with PGS main effect, two independent environmental variables $E_1$ (binary) and $E_2$ (continuous), and their interaction effect with PGS. Assume the marginal disease prevalence is Pr(D=1)=0.01.
```

#example analysis
set.seed(04262023)

#Generate 200000 families with offsprings disease prevalence to be 1%
dat = sim_prospective_population(n_fam=200000, #Number of families in the population
                               cor_e_prs=FALSE, #Assuming no correlation between PGS and E
                               cor_strat=0.25, #Random effect term to create population stratification bias in PGS
                               rho2=0.25, #Random effect term to create population stratification bias in PGSxE
                               rho_mf=0, #No assortative mating
                               alpha_fam=-5.5, #Intercept
                               betaG_normPRS=0.4, #Main effect of PGS
                               betaE_bin=0.2, #Main effect of E1
                               betaE_norm=0.2, #Main effect of E2
                               betaGE_normPRS_bin=0, #Interaction effect of PGSxE1
                               betaGE_normPRS_norm=0, #Interaction effect of PGSxE2
                               envir=TRUE) #Include environmental variables in the disease risk model


```

View the simulated data.
```
head(dat$E_sim)
#     E_sim_bin  E_sim_norm
#[1,]         0 -0.79841223
#[2,]         0 -0.26671287
#[3,]         0 -0.08344659
#[4,]         1  2.05270642
#[5,]         0 -0.51954546
#[6,]         0  0.02629269

head(dat$pgs_fam)
#          pgs_c       pgs_m       pgs_f
#[1,] -0.7506135 -2.04896031 -0.56878873
#[2,] -1.1159313 -0.04741079 -1.00756479
#[3,]  0.4394364 -0.76419757  0.01806875
#[4,]  0.2408536  0.99687805  0.07513997
#[5,]  0.7600675  1.59946273 -0.05412554
#[6,] -0.4180675  0.29314448  0.14294136

table(dat$D_sim)

#     0      1 
#197822   2178 
```

Randomly select 1000 affected probands and their families
```
id = sample(which(dat$D_sim==1),size=1000,replace=F)
PRS_fam = dat$pgs_fam[id,]
envir = dat$E_sim[id,]
```

Fit our proposed method to the randomly selected 1000 trios
```
startTime <- Sys.time()
res_sim = pgs.cpt(pgs_offspring = PRS_fam[,1], 
                 pgs_mother = PRS_fam[,2], 
                 pgs_father = PRS_fam[,3],
                 GxE_int = TRUE, #If GxE_int is FALSE, then fit a model without interaction effect between PGSxE. "formula" and "E" will be ignored in the function.
                 formula = ~ factor(E_sim_bin)+E_sim_norm, #For categorical variables, remember to add factor().
                 E = envir, 
                 side = 2,
                 numDeriv = F)

endTime <- Sys.time()
```

Print the final results.
```
res_sim$res_beta
#                           Estimate   Std.Error   Z.value   Pvalue
# PGS                       0.38144034 0.10852411  3.5147981 0.0004400885
# PGS x factor(E_sim_bin)1 -0.09540741 0.13803817 -0.6911669 0.4894606688
# PGS x E_sim_norm          0.06266047 0.06394153  0.9799652 0.3271033111
```

Print running time.
```
print(endTime - startTime)
#Time difference of 0.005946875 secs
```

# Questions
Please send your questions/suggestions to zwang389@jhu.edu
