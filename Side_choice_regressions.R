# 1) Load data ------------------------------------------------------------

setwd("Data/")
# read side.choice dataframe
load("side.choice.input.RDA")

df.side.choice <- puzzle.data.combined.sub

# subset to fledgling solves only
df.side.choice <- subset(df.side.choice, df.side.choice$age=="Fledgling")

head(df.side.choice)

# subset to those with a minimum of 10 solves
IDs <- c()
for(i in unique(df.side.choice$PIT)){
  if(max(subset(df.side.choice$bout.ind, df.side.choice$PIT==i))>=10){
    IDs <- c(IDs,i)
  }
}

df.side.choice <- subset(df.side.choice, df.side.choice$PIT %in% IDs)

head(df.side.choice)
# column names explained:

# PIT: alphanumeric code of individual
# Event: 'right' for a solve to the right, 'left' for a solve to the left
# Date: yyyy-mm-dd
# Time: HH:MM:SS
# Location: location of puzzle box (Mill_1-Mill_4, Guett_1-Guett_3)
# Date.Time: yymmddHHMMSS
# age: 'Fledgling' for first-year birds
# time.since.fl: number of days since fledging when the solve occurred
# parent.frq.right: number of right solves by parents in the five minutes prior to the respective solve
# parent.frq.left: number of left solves by parents in the five minutes prior to the respective solve
# parents.frq.total.right: number of right solves that parents did between fledging and the fledglings first solves (aka right solves during dependence)
# parents.frq.total.left: number of left solves that parents did between fledging and the fledglings first solves (aka left solves during dependence)
# nonparent.frq.right: number of right solves by non-parent adults in the five minutes prior to the respective solve
# nonparent.frq.left: number of left solves by non-parent adults in the five minutes prior to the respective solve
# siblings.frq.right: number of right solves by siblings in the five minutes prior to the respective solve
# siblings.frq.left: number of left solves by siblings in the five minutes prior to the respective solve
# nonsiblings.frq.right: number of right solves by non-siblings (peers) in the five minutes prior to the respective solve
# nonsiblings.frq.left: number of left solves by non-siblings (peers) in the five minutes prior to the respective solve
# bout.ind: cumulative number of solves by individual
# bout.ind.right: cumulative number of right solves by individual
# bout.ind.left: cumulative number of left solves by individual

# we add a few more columns: 

# for each solve, we extract number of solves that have occurred in the 5 mins prior
# and the solving day for each fledgling (i.e. how many days of solving experience it has acquired)
for(i in 1:nrow(df.side.choice)){
  side <- df.side.choice[i, "Event"]
  num.right <- df.side.choice[i, "frq.right"]
  num.left <- df.side.choice[i, "frq.left"]
  days.solving <- unique(subset(df.side.choice$Date, df.side.choice$PIT==df.side.choice[i, "PIT"]))
  date <- df.side.choice[i, "Date"]
  solving.day <- which(days.solving==date)
  parent.freq.right.dependence <- df.side.choice[i, "parent.frq.total.right"]/ (df.side.choice[i, "parent.frq.total.right"]+ df.side.choice[i, "parent.frq.total.left"])
  
  
  site <- substr(df.side.choice[i, "Location"], 1, 2)
  
  # add whether they had conflicting info from parents
  if(is.na(parent.freq.right.dependence)){
    conflicting <- NA
  } else if (parent.freq.right.dependence>=0.9){
    conflicting <- "same_right"
  } else if (parent.freq.right.dependence<=0.1){
    conflicting <- "same_left"
  } else {
    conflicting <- "opp"
  }
    
  
  if(side=="right"){
    frq <- num.right / (num.right+num.left)
    num <- num.right
  } else { 
    frq <- num.left / (num.right+num.left)
  num <- num.left
  }
  df.side.choice[i, "frequency.side.prior"] <- frq
  df.side.choice[i, "solves.prior"] <- num.right+num.left
  df.side.choice[i, "num.side.prior"] <- num
  df.side.choice[i, "solving.day"] <- solving.day
  df.side.choice[i, "frequency.right.prior"] <- num.right/(num.right+num.left)
  df.side.choice[i, "parent.freq.dependence.right"] <- parent.freq.right.dependence 
  df.side.choice[i, "parent.solves.dependence"] <- df.side.choice[i, "parent.frq.total.right"]+ df.side.choice[i, "parent.frq.total.left"]
  df.side.choice[i, "conflicting"] <- conflicting
  df.side.choice[i, "site"] <- site
  }

head(df.side.choice)
# new columns are:
# frequency.side.prior: the proportion of solves during the last 5 minutes occurring on the same side that the fledgling solved
# solves.prior: number of solves during the prior 5 minutes
# num.side.prio: nuumber of solves during the 5 minutes prior on the same side that the fledgling solved
# solving.day: the number of days of solving experience the fledglings have acquired
# frequency.right.prior: the proportion of solves on the right during the last 5 minutes
# parent.freq.dependence.right: the proportion of solves on the right during dependence (from fledgling till the fledglings first solve) performed by parents
# parent.solves.dependence: the number of solves by parents during dependence
# conflicting: whether the parents showed the same solving behaviour ('no') (i.e. more than 90% of solves occurred on the same side) or conflicting 'yes' if more than 10% of solves occurred on opposite sides
# site: study site - Mi or Gu


# 2) First day - conflicting or non-conflicting ----------------------------

# In a first model, we test whether the fledglings initial solving behaviour was influenced by whether they received conflicting, non-conflicting (RR or LL) or no information from their parents during dependence
# If this is the case, we expect fledglings with parents that solved on opposite sides, or did not solve at all to be more mixed in their side choice initially

# we subset the data to the first solving day
df.side.choice.confl <- subset(df.side.choice, df.side.choice$solving.day==1)

df.side.choice.confl$conflicting
# replace the NAs with "none" (for non-solving parents)

df.side.choice.confl$conflicting[is.na(df.side.choice.confl$conflicting)] <- "none"

# conflicting summary
confl.summary <- table(cbind.data.frame(df.side.choice.confl$PIT, df.side.choice.confl$Event))
confl.summary <- cbind.data.frame(rownames(confl.summary),as.numeric(confl.summary[,1]), as.numeric(confl.summary[,2]))
colnames(confl.summary) <- c("PIT", "left", "right")

confl.summary$conflicting <- NULL

for(i in confl.summary$PIT){
  confl.i <- unique(subset(df.side.choice.confl$conflicting, df.side.choice.confl$PIT==i))
  num.parental.solves.i <-unique(subset(df.side.choice.confl$parent.solves.dependence, df.side.choice.confl$PIT==i))
  confl.summary[which(confl.summary$PIT==i), "conflicting"] <- confl.i
  confl.summary[which(confl.summary$PIT==i), "num.parental.solves"] <- num.parental.solves.i
}

# let's add a total number of solves on the first solving day:
confl.summary$total <- confl.summary$left+confl.summary$right

confl.summary$prop.right <- confl.summary$right/confl.summary$total

# and a column of the log number of parental solves
confl.summary$log_parental_solves <- log(confl.summary$num.parental.solves+1)

length(confl.summary$PIT)
# 53 individuals


# Let's run the model
library(brms)
set.seed(5)
model.conflicting <-
  brms::brm(
    prop.right ~ conflicting,
    confl.summary,
    family=gaussian(),
    chains=4,
    iter=4000,
    cores = 4
  )

#save(model.conflicting, file="Output/model_conflicting.RData")
load("Output/model_conflicting.RData")


# # model performance
plot(model.conflicting)
# we can see that chains have mixed well, and they have reached stationarity

# and posterior predictive checks
pp_check(model.conflicting, ndraws= 1e2)
# the dark blue line shows the distribution of the real data, the light blue line shows 100 draws of the posterior. Our model does not explain the data well - meaning that parental side choices during dependence do not predict side choices on the juveniles' first solving day well. 

# Next we look at the summary
summary(model.conflicting)

# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: prop.right ~ conflicting 
# Data: confl.summary (Number of observations: 53) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                 0.58      0.22     0.15     1.01 1.00     3536     4017
# conflictingopp           -0.17      0.25    -0.67     0.33 1.00     3850     4697
# conflictingsame_left     -0.12      0.24    -0.60     0.36 1.00     3865     4758
# conflictingsame_right     0.09      0.24    -0.39     0.57 1.00     3934     4668
# 
# Further Distributional Parameters:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.44      0.05     0.36     0.54 1.00     4899     4691
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# we transform the log-odds into odds by exponentiating
round(exp(fixef(model.conflicting)), 2)

#                        Estimate Est.Error Q2.5 Q97.5
# Intercept                 1.78      1.25 1.16  2.74
# conflictingopp            0.84      1.29 0.51  1.39
# conflictingsame_left      0.89      1.28 0.55  1.43
# conflictingsame_right     1.10      1.27 0.68  1.76

#finally we transform into probabilities
round(plogis(fixef(model.conflicting)),2)

# Estimate Est.Error Q2.5 Q97.5
# Intercept                 0.64      0.55 0.54  0.73
# conflictingopp            0.46      0.56 0.34  0.58
# conflictingsame_left      0.47      0.56 0.36  0.59
# conflictingsame_right     0.52      0.56 0.40  0.64

# compare means between categories:

library(emmeans)

em <- emmeans(model.conflicting, ~ conflicting)
contrast(em, method = "pairwise")

# contrast               estimate lower.HPD upper.HPD
# none - opp               0.1716    -0.302    0.6876
# none - same_left         0.1129    -0.335    0.6124
# none - same_right       -0.0962    -0.560    0.3950
# opp - same_left         -0.0578    -0.365    0.2624
# opp - same_right        -0.2656    -0.591    0.0545
# same_left - same_right  -0.2092    -0.518    0.0675
# 
# Point estimate displayed: median 
# HPD interval probability: 0.95 

# visualize
plot.conflict <- plot(conditional_effects(model.conflicting), points=TRUE, point_args = list(width = 0.2, height=0.01, alpha=0.2))

library(ggplot2)
plot.confl <- plot.conflict$conflicting+
  theme_bw()+
  ylab("Prop. right solves on first solving day")+
  xlab("Parental solving behavior during dependence")+
  scale_x_discrete(labels=c("none", "conflicting", "both left", "both right"))+
   ggtitle("A)") +
   theme(text = element_text(size=14))+
   theme(plot.title = element_text(hjust = -0.08, size=14))

# what are the numbers?

table(confl.summary$conflicting)

# none        opp  same_left same_right 
# 4         13         18         18 

# 3) First solving day - social influence ------------------------------------------

df.side.choice.first.day <-
  subset(
    df.side.choice,
      df.side.choice$solving.day == 1
  )

dim(df.side.choice.first.day)
max(table(df.side.choice.first.day$PIT))
# 51
mean(table(df.side.choice.first.day$PIT))
# [1] 3.981132

# reduce to events where social information is available (these will be thrown out in the model anyway)
df.side.choice.first.day <- df.side.choice.first.day[(df.side.choice.first.day$solves.prior>0),]

dim(df.side.choice.first.day)
# 196

length(unique(df.side.choice.first.day$PIT))
# 51


set.seed(5)

# model first day
model.first.day <-
  brms::brm(
    Event ~ frequency.right.prior*scale(solves.prior) + (1|PIT),
    df.side.choice.first.day,
    family=bernoulli(link="logit"),
    chains=4,
    iter=4000,
    cores = 4
#    control = list(adapt_delta=0.99)
  )


# save(model.first.day, file="Output/model_first_day.RData")
load("Output/model_first_day.RData")


# model performance
plot(model.first.day)


# and posterior predictive checks
pp_check(model.first.day, ndraws= 1e2)

# Next we look at the summary
summary(model.first.day)
# Family: bernoulli 
# Links: mu = logit 
# Formula: Event ~ frequency.right.prior * scale(solves.prior) + (1 | PIT) 
# Data: df.side.choice.first.day (Number of observations: 196) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Multilevel Hyperparameters:
#   ~PIT (Number of levels: 51) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.15      0.39     0.53     2.04 1.00     3155     4602
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                  -1.37      0.38    -2.14    -0.61 1.00     8984     6298
# frequency.right.prior                       2.98      0.62     1.76     4.20 1.00    11336     6719
# scalesolves.prior                          -0.69      0.36    -1.44    -0.01 1.00     9288     6283
# frequency.right.prior:scalesolves.prior     0.87      0.65    -0.37     2.19 1.00     8510     5376
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# exponentiate to get the odds
round(exp(fixef(model.first.day)), 2)

# Estimate Est.Error Q2.5 Q97.5
# Intercept                                   0.26      1.48 0.12  0.55
# frequency.right.prior                      19.50      1.87 5.84 68.62
# scalesolves.prior                           0.50      1.44 0.24  0.99
# frequency.right.prior:scalesolves.prior     2.39      1.88 0.72  8.62

plogis(fixef(model.first.day))
# Estimate Est.Error      Q2.5     Q97.5
# Intercept                               0.2025221 0.5948698 0.1056177 0.3513971
# frequency.right.prior                   0.9516578 0.6508645 0.8536054 0.9852478
# scalesolves.prior                       0.3349179 0.5901108 0.1917225 0.4968341
# frequency.right.prior:scalesolves.prior 0.7042587 0.6559042 0.4090916 0.8992290


# we then extract the conditional effects for the proportion of right solves of other in the 5 mins prior

first.day.plot <- plot(conditional_effects(model.first.day))

# we translate the 'right' and 'left' in the original data to 0 and 1 for it to be on the same scale as the model output
df.side.choice.first.day$Event.no <- df.side.choice.first.day$Event
df.side.choice.first.day$Event.no[df.side.choice.first.day$Event.no=="right"] <- 0
df.side.choice.first.day$Event.no[df.side.choice.first.day$Event.no=="left"] <- 1
df.side.choice.first.day$Event.no <- as.numeric(df.side.choice.first.day$Event.no)


# we then plot both the model output and raw data
plot.first.day.social <- first.day.plot$frequency.right.prior+
  theme_bw()+
  geom_line(col="black", lwd=1.02)+
  geom_point(df.side.choice.first.day, mapping=aes(x=frequency.right.prior, y=Event.no), inherit.aes = FALSE, position=position_jitter(w = 0.02, h = 0.06), alpha=0.2, size=2, color="grey15") +
  scale_y_continuous(name="Side choice first solving day", limits=c(-0.1, 1.1), breaks=c(0,1), labels =c("left", "right"))+
  scale_x_continuous(name="Prop. of right 5 solves mins prior", limits=c(-0.02, 1.02), breaks = c(0, 0.5, 1), labels=c("0", "0.5", "1"))+
  geom_hline(yintercept=0.5, lty="dashed")+
  ggtitle("B)") +
  theme(text = element_text(size=14))+
  theme(plot.title = element_text(hjust = -0.12, size=14))



# 4) Influence of asocial and social information  --------------------------
# we are now considering all solving days

# we add a column to the data frame adding the number of right and left solves the fledgling produced an hour before solving
df.side.choice$total.self.60 <- df.side.choice$frq.left.self.60+df.side.choice$frq.right.self.60

df.side.choice.soc.vs.asoc <- df.side.choice

# subset to those events where they have asocial and social information available
 df.side.choice.soc.vs.asoc <- subset(df.side.choice, df.side.choice$total.self.60>=1 & df.side.choice$solves.prior>=1)

df.side.choice.soc.vs.asoc$Event[1]
# left is baseline

dim(df.side.choice.soc.vs.asoc)
# 8041    events where both social and asocial information are available

length(unique(df.side.choice.soc.vs.asoc$PIT))
# 53


# make solves a proportion
df.side.choice.soc.vs.asoc$frq.right.self.60 <- df.side.choice.soc.vs.asoc$frq.right.self.60/df.side.choice.soc.vs.asoc$total.self.60

#df.side.choice.soc.vs.asoc$frequency.right.prior <- df.side.choice.soc.vs.asoc$frequency.right.prior

df.side.choice.soc.vs.asoc$bout.ind <- df.side.choice.soc.vs.asoc$bout.ind/100

head(df.side.choice.soc.vs.asoc)


# model soc.vs.asoc
set.seed(5)
model.soc.vs.asoc <-
  brms::brm(
    Event ~ frq.right.self.60*bout.ind +  frequency.right.prior*bout.ind + (1|PIT),
    df.side.choice.soc.vs.asoc,
    family=bernoulli(link="logit"),
    chains=4,
    iter=4000,
    cores = 4
#    control = list(adapt_delta=0.99)
  )

#save(model.soc.vs.asoc, file="Output/model_soc_vs_asoc.RData")
load("Output/model_soc_vs_asoc.RData")


# model performance
plot(model.soc.vs.asoc)

# and posterior predictive checks
pp_check(model.soc.vs.asoc, ndraws= 1e2)

# Next we look at the summary
 summary(model.soc.vs.asoc)
 # Family: bernoulli 
 # Links: mu = logit 
 # Formula: Event ~ frq.right.self.60 * bout.ind + frequency.right.prior * bout.ind + (1 | PIT) 
 # Data: df.side.choice.soc.vs.asoc (Number of observations: 8041) 
 # Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
 # total post-warmup draws = 8000
 # 
 # Multilevel Hyperparameters:
 #   ~PIT (Number of levels: 53) 
 # Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
 # sd(Intercept)     0.93      0.14     0.68     1.23 1.00     1978     3303
 # 
 # Regression Coefficients:
 #   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
 # Intercept                         -2.26      0.17    -2.60    -1.92 1.00     1231     2505
 # frq.right.self.60                  1.99      0.19     1.62     2.35 1.00     4609     5953
 # bout.ind                          -0.04      0.03    -0.11     0.02 1.00     5214     5268
 # frequency.right.prior              1.42      0.19     1.06     1.79 1.00     5219     5888
 # frq.right.self.60:bout.ind         0.20      0.08     0.04     0.36 1.00     5247     5954
 # bout.ind:frequency.right.prior    -0.28      0.06    -0.40    -0.15 1.00     5180     6103
 # 
 # Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
 # and Tail_ESS are effective sample size measures, and Rhat is the potential
 # scale reduction factor on split chains (at convergence, Rhat = 1).
 
 
 
# get odds by exponentiating
round(exp(fixef(model.soc.vs.asoc)),2)
#                                 Estimate Est.Error Q2.5 Q97.5
# Intercept                          0.10      1.19 0.07  0.15
# frq.right.self.60                  7.30      1.20 5.07 10.49
# bout.ind                           0.96      1.03 0.90  1.02
# frequency.right.prior              4.15      1.20 2.89  5.98
# frq.right.self.60:bout.ind         1.22      1.09 1.04  1.44
# bout.ind:frequency.right.prior     0.76      1.07 0.67  0.86

# get probabilities of producing a left solve
plogis(fixef(model.soc.vs.asoc))
# Estimate Est.Error       Q2.5     Q97.5
# Intercept                      0.09413136 0.5435444 0.06914428 0.1283033
# frq.right.self.60              0.87946181 0.5463651 0.83524827 0.9129780
# bout.ind                       0.48935396 0.5079453 0.47327020 0.5046390
# frequency.right.prior          0.80594790 0.5464535 0.74282266 0.8567593
# frq.right.self.60:bout.ind     0.54904130 0.5208303 0.50931125 0.5895810
# bout.ind:frequency.right.prior 0.43101537 0.5161563 0.40120903 0.4621383


# set all conditions to average for plotting
# apart from the variables of itnerest: first, bout.ind. We will set the values corresponding to 0, 300, 600, 900 solves 

bout.ind <- c(0, 3, 6, 9)


int_conditions <- list(
  bout.ind = setNames(c(bout.ind), c(bout.ind)),
  frequency.right.prior = setNames(c(0.5), c(0.5)))


# we then extract the conditional effects for the frequency of right solves (personal info):experience (measured as individual cumulative solving days)

asoc.plot <- plot(conditional_effects(model.soc.vs.asoc, int_conditions = int_conditions))

# we translate the 'right' and 'left' in the original data to 0 and 1 for it to be on the same scale as the model output
df.side.choice.soc.vs.asoc$Event.no <- df.side.choice.soc.vs.asoc$Event
df.side.choice.soc.vs.asoc$Event.no[df.side.choice.soc.vs.asoc$Event.no=="right"] <- 0
df.side.choice.soc.vs.asoc$Event.no[df.side.choice.soc.vs.asoc$Event.no=="left"] <- 1
df.side.choice.soc.vs.asoc$Event.no <- as.numeric(df.side.choice.soc.vs.asoc$Event.no)

my_palette <- rev(RColorBrewer::brewer.pal(name="Blues",n=9)[c(3,5,7,9)])


# we then plot both the model output and raw data

asoc.plot2 <- asoc.plot$`frq.right.self.60:bout.ind`+
  theme_bw()+
  geom_line(lwd=1.02)+
  # geom_point(df.side.choice.soc.vs.asoc, mapping=aes(x=frequency.right.prior, y=Event.no), inherit.aes = FALSE, position=position_jitter(w = 0.01, h = 0.06), alpha=0.2, size=1.5)+
  #  geom_point(df.side.choice.social, mapping=aes(x=frequency.right.prior, y=Event.no), color="grey15", inherit.aes = FALSE, position=position_jitter(w = 0.02, h = 0.08), alpha=0.6, size=2)+
  scale_y_continuous(name="Side choice overall", limits=c(-0.1, 1.1), breaks=c(0,1), labels =c("left", "right"))+
  #    scale_y_continuous(name="", limits=c(0, 1), breaks=c(0,1), labels =c("", ""))+
  #    scale_x_continuous(name=c("Prop. right solves 60 mins prior \n(personal information)"), limits=c(0, 1), breaks = c(0, 0.5, 1), labels=c("0", "0.5", "1"))+
  scale_x_continuous(name=c("Prop. of right solves 60 mins prior \n(personal information)"), limits=c(0, 1), breaks = c(0, 0.5, 1), labels=c("0", "0.5", "1"))+
  theme(legend.position = "none")+
  ggtitle("C)") +
  theme(text = element_text(size=14))+
  geom_hline(yintercept=0.5, lty="dashed")+
  theme(plot.title = element_text(hjust = -0.12, size=14))+
  scale_fill_manual(name = "Personal experience \n[# cumulative solves]", breaks=c(bout.ind), values = rev(my_palette), labels=c("0","300","600", "900"))+
  scale_color_manual(name = "Personal experience \n[# cumulative solves]", breaks=c(bout.ind), values = rev(my_palette), labels=c("0","300","600", "900"))+
  theme(legend.position= c(.3, .75))


# next, we look at the influence of social information with increasing experience
int_conditions <- list(
  bout.ind = setNames(c(bout.ind), c(bout.ind)),
  frq.right.self.60 = setNames(c(0.5), c(0.5)))

soc.plot <- plot(conditional_effects(model.soc.vs.asoc, int_conditions = int_conditions, "frequency.right.prior:bout.ind"))

soc.plot2 <-  soc.plot$`frequency.right.prior:bout.ind`+
  theme_bw()+
  geom_line(lwd=1.02)+
  # geom_point(df.side.choice.soc.vs.asoc, mapping=aes(x=frequency.right.prior, y=Event.no), inherit.aes = FALSE, position=position_jitter(w = 0.01, h = 0.06), alpha=0.2, size=1.5)+
#  scale_fill_manual(breaks=c(bout.ind), values = rev(my_palette))+
#  scale_color_manual(breaks=c(bout.ind), values = rev(my_palette))+
  #  geom_point(df.side.choice.social, mapping=aes(x=frequency.right.prior, y=Event.no), color="grey15", inherit.aes = FALSE, position=position_jitter(w = 0.02, h = 0.08), alpha=0.6, size=2)+
  scale_y_continuous(name="", limits=c(-0.1, 1.1), breaks=c(0,1), labels =c("left", "right"))+
  #scale_y_continuous(name="Probability of producing a right solve", limits=c(0, 1), breaks=c(0,0.5,1), labels =c(0, 0.5, 1))+
  scale_x_continuous(name=c("Prop. of right solves 5 mins prior \n(social information)"), limits=c(0, 1), breaks = c(0, 0.5, 1), labels=c("0", "0.5", "1"))+
  ggtitle("D)") +
  theme(text = element_text(size=14))+
  geom_hline(yintercept=0.5, lty="dashed")+
  theme(plot.title = element_text(hjust = -0.12, size=14))+
  scale_fill_manual(name = "Personal experience \n[# cumulative solves]", breaks=c(bout.ind), values = rev(my_palette), labels=c("0","300","600", "900"))+
  scale_color_manual(name = "Personal experience \n[# cumulative solves]", breaks=c(bout.ind), values = rev(my_palette), labels=c("0","300","600", "900"))+
  theme(legend.position = c(.3, .75))


# 
# # not plot again without the legend
# soc.plot3 <-  soc.plot$`frequency.right.prior:bout.ind`+
#   theme_bw()+
#   geom_line(lwd=1.02)+
#   # geom_point(df.side.choice.soc.vs.asoc, mapping=aes(x=frequency.right.prior, y=Event.no), inherit.aes = FALSE, position=position_jitter(w = 0.01, h = 0.06), alpha=0.2, size=1.5)+
#   #  scale_fill_manual(breaks=c(bout.ind), values = rev(my_palette))+
#   #  scale_color_manual(breaks=c(bout.ind), values = rev(my_palette))+
#   #  geom_point(df.side.choice.social, mapping=aes(x=frequency.right.prior, y=Event.no), color="grey15", inherit.aes = FALSE, position=position_jitter(w = 0.02, h = 0.08), alpha=0.6, size=2)+
#   scale_y_continuous(name="", limits=c(-0.1, 1.1), breaks=c(0,1), labels =c("left", "right"))+
#   #scale_y_continuous(name="Probability of producing a right solve", limits=c(0, 1), breaks=c(0,0.5,1), labels =c(0, 0.5, 1))+
#   scale_x_continuous(name=c("Prop. right solves 5 mins prior \n(social information)"), limits=c(0, 1), breaks = c(0, 0.5, 1), labels=c("0", "0.5", "1"))+
#   ggtitle("D)") +
#    theme(legend.position = "none")+
#   theme(text = element_text(size=14))+
#   geom_hline(yintercept=0.5, lty="dashed")+
#   theme(plot.title = element_text(hjust = -0.12, size=14))+
#   scale_fill_manual(name = "Personal experience \n[# cumulative solves]", breaks=c(bout.ind), values = rev(my_palette), labels=c("0","300","600", "900"))+
#   scale_color_manual(name = "Personal experience \n[# cumulative solves]", breaks=c(bout.ind), values = rev(my_palette), labels=c("0","300","600", "900"))
# 
# library(cowplot)
# 
# # we extract the legend for plot 2 
# legend <- get_legend(soc.plot2)

library(ggpubr)

png(filename="C:/Users/sonja/Desktop/Konstanz/Breeding Season 2021/Fledgie social learning/Figures/side_choices.png",width=1800, height=1500, res=170)

# to left align the legend, we pretend we're plotting a second legend of soc plot 3 (which has legend set to "none")
ggarrange(plot.confl, plot.first.day.social, asoc.plot2, soc.plot2, nrow=2, ncol=2)


dev.off()



# how long from first to third solution
all.d <- NULL

for(i in unique(df.side.choice$PIT)){
  di <- subset(df.side.choice$Date, df.side.choice$PIT==i)[c(1,3)]
  di.data <- c(i, di[1], di[2])
  all.d <- rbind.data.frame(all.d, di.data)
}

colnames(all.d) <- c("PIT", "first", "third")

median(as.POSIXct(all.d$third)- as.POSIXct(all.d$first))/(60*60*24)
