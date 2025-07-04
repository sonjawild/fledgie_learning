
# 1) Does parental behaviour predict fledgling strategy -------------------


# 1.1) Load data ----------------------------------------------------------

setwd("Data/")
all.fledgies <- read.csv("fledgling_parents_strategies.csv", row.names = 1)

head(all.fledgies)
# columns:
# - PIT: unique PIT code for each fledgling
# - strategy: whether they learned to 'solve', 'scrounge', 'both', or 'none'
# - strat.solve: whether they learned to solve (yes) or not (no)
# - strat.scrounge: whether they learned to scrounge (yes) or not (no)
# - parents.solve: did fledglings have knowledgeable solving parents (1) or not (0)
# - parents.scorunge: did fledglings have knowledgeable scrounging parents (1) or not (0)
# - parents.num.solves: how many solves have parents produced between fledging and date of acquisition (solving)
# - parents.num.scrounges: how many scrounges have parents produced between fledging and date of acquisition (scrounging)
# - site: either 'Mill' or 'Guett'

# we add a column for log number of solves and scrounges
all.fledgies$parents_solves_log <- log(all.fledgies$parents.num.solves+1)
all.fledgies$parents_scrounges_log <- log(all.fledgies$parents.num.scrounges+1)


# 1.2) Model ------------------------------------------------------------

# load libraries
library(brms)
library(rstan)


set.seed(5)
model.parent.predict <-
  brms::brm(
    strat.solve ~  parents_solves_log + parents_scrounges_log,
    all.fledgies,
    family = bernoulli(),
    chains=4,
    iter=6000,
    cores = 4,
    prior=   c(prior(normal(0, 10), class= Intercept),
               prior(normal(0,  5), class= b)),
     control= list(adapt_delta= .98, max_treedepth= 15)
  )



# save(model.parent.predict, file="Output/model_parent_predict.RData")
# model output can be loaded here directly
load("Output/model_parent_predict.RData")

# confirm stationarity and mixing and fit of model
plot(model.parent.predict)


summary(model.parent.predict)
# Family: bernoulli 
# Links: mu = logit 
# Formula: strat.solve ~ parents_solves_log + parents_scrounges_log 
# Data: all.fledgies (Number of observations: 103) 
# Draws: 4 chains, each with iter = 6000; warmup = 3000; thin = 1;
# total post-warmup draws = 12000
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                -2.11      0.55    -3.30    -1.15 1.00    10036     7497
# parents_solves_log        5.62      1.82     2.54     9.62 1.00     1652     2832
# parents_scrounges_log    -0.71      0.52    -1.99     0.03 1.00     2826     2948
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# convert to odds by exponentiating:
round(exp(fixef(model.parent.predict)),2)
#                          Estimate Est.Error  Q2.5    Q97.5
# Intercept                 0.12      1.73  0.04     0.32
# parents_solves_log      274.84      6.20 12.66 15116.01
# parents_scrounges_log     0.49      1.69  0.14     1.03

# how to interpret? 
# values above 1 increase the probability, while values below 1 decrease the probability
# we also check that the confidence itnervals do not span 1
# e.g. parents_solves_log: the odds for the fledgling to be a solver increase by a factor of 269 [12,14,478] by increase in magnitude of parental solves

# we can also convert to probabilities of an on offspring being a solver (most useful for the intercept):
plogis(fixef(model.parent.predict))
#                        Estimate Est.Error       Q2.5     Q97.5
# Intercept             0.1080875 0.6331581 0.03562384 0.2406519
# parents_solves_log    0.9963747 0.8610887 0.92679934 0.9999338
# parents_scrounges_log 0.3302137 0.6281010 0.12051121 0.5068164

# Intercept: probability of an offspring being a solver when both number of parental scrounges and solves are 0



