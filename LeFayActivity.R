###### EEB 622 in-class activity ######
###### Solai Le Fay

### load relevant packages ###
library(tidyverse)
library(brms)
library(bayesplot)
library(marginaleffects)
library(performance)
library(bayesplot)


### Load data ###
recdata <- read.csv("recordings.csv")
sensdata <- read.csv("sensorinfo.csv")


### Research Questions ###
# 1. Does sound from boats influence beluga whale song duration and frequency (total calls/hr)? 
# 2. Do boats have other direct impacts on whale calling behavior, outside of their noise impacts?


### Q1 model variables - 2 models ### 
### Main effect of interest is boat NOISE
### totsongs ~ boatnoise + watertemp + boatactivity + (1|sensorid) + (1|dayid)
# Poisson OR negative binomial distribution

### songlength ~ boatnoise + watertemp + boatactivity + (1|sensorid) + (1|dayid)
# Gamma distribution


### Q2 model variables - 2 models ###
### Main effect of interest is boat ACTIVITY
# totsongs ~ boatactivity + waterdepth + distshore + (1|sensorid) + (1|dayid)
# Poisson OR negative binomial distribution

# songlength ~ boatactivity + waterdepth + distshore + (1|sensorid) + (1|dayid)
# Gamma distribution



###### Analyze Data ######
# Merge data frames
# Convert sensorID to categorical


merged_dat <- recdata %>% left_join(sensdata, by = "sensorid")
merged_dat$sensorid <- as.character(merged_dat$sensorid)

merged_dat$watertemp <- as.numeric(merged_dat$watertemp)

#remove one NA row
merged_dat <- merged_dat[-10, ]


#### Run models ####
# Research question 1
# poisson first, check for over dispersion, change to nb if over dispersed
mod1_totsongs <- brm(totsongs ~ scale(boatnoise) + scale(watertemp) + scale(boatactivity) + (1|sensorid) + (1|dayid), data = merged_dat, family = poisson (link = "log"))

#gamma
mod1_songlength <- brm(songlength ~ scale(boatnoise) + scale(watertemp) + scale(boatactivity) + (1|sensorid) + (1|dayid), data = merged_dat, family = Gamma (link = "log"))


# Research question 2
# poisson first, check for over dispersion, change to nb is over dispersed
mod2_totsongs <- brm(totsongs ~ scale(boatactivity) + scale(waterdepth) + scale(distshore) + scale(boatnoise) + (1|sensorid) + (1|dayid), data = merged_dat, family = poisson(link = "log"))

#gamma
mod2_songlength <- brm(songlength ~ scale(boatactivity) + scale(waterdepth) + scale(distshore) + scale(boatnoise) + (1|sensorid) + (1|dayid), data = merged_dat, family = Gamma (link = "log"))


### Check both poisson models for over dispersion using pp checks
#Create a function for dispersion
dispersion <- function(x) {var(x)/mean(x)}

pp_check(mod1_totsongs) + xlim(c(0, 300))
pp_check(mod2_totsongs) + xlim(c(0, 300))

ppc_stat(y = merged_dat$totsongs,
         # Compare the dispersion in the real data (y)
         yrep = posterior_predict(mod1_totsongs,
                                  ndraws = 1000),
         # to dispersion in predictions from our posterior (yrep)
         stat="dispersion")

ppc_stat(y = merged_dat$totsongs,
         # Compare the dispersion in the real data (y)
         yrep = posterior_predict(mod2_totsongs,
                                  ndraws = 1000),
         # to dispersion in predictions from our posterior (yrep)
         stat="dispersion")


### Data over dispersed! Rerun models with negative binomial distribution
nbmod1_totsongs <- brm(totsongs ~ scale(boatnoise) + scale(watertemp) + scale(boatactivity) + (1|sensorid) + (1|dayid), data = merged_dat, family = negbinomial (link = "log"))

nbmod2_totsongs <- brm(totsongs ~ scale(boatactivity) + scale(waterdepth) + scale(distshore) + scale(boatnoise) + (1|sensorid) + (1|dayid), data = merged_dat, family = negbinomial (link = "log"))


# re-check dispersion
pp_check(nbmod1_totsongs) + xlim(c(0, 300))
pp_check(nbmod2_totsongs) + xlim(c(0, 300))

ppc_stat(y = merged_dat$totsongs,
         # Compare the dispersion in the real data (y)
         yrep = posterior_predict(nbmod1_totsongs,
                                  ndraws = 1000),
         # to dispersion in predictions from our posterior (yrep)
         stat="dispersion")

ppc_stat(y = merged_dat$totsongs,
         # Compare the dispersion in the real data (y)
         yrep = posterior_predict(nbmod2_totsongs,
                                  ndraws = 1000),
         # to dispersion in predictions from our posterior (yrep)
         stat="dispersion")


#check convergence for final models
plot(nbmod1_totsongs)
plot(mod1_songlength)
plot(mod2_songlength)
plot(nbmod2_totsongs)
#all converged
#4 chains, 2000 iterations, 1000 are warmup


#### Interpret Results - Research Q 1 ####
#### Main effect of interest is boat NOISE
summary(nbmod1_totsongs)
summary(mod1_songlength)

bayes_R2(nbmod1_totsongs) #0.54
bayes_R2(mod1_songlength) #0.33

performance_mae(nbmod1_totsongs) #21.3
performance_mae(mod1_songlength) #8.63

min(merged_dat$totsongs) #0
max(merged_dat$totsongs) #777
mean(merged_dat$totsongs) #26.99
hist(merged_dat$totsongs) #the max of 777 may be an error that is supposed to say 77, but I do not know whale biology well enough to know if 777 is a red flag or reasonable.

min(merged_dat$songlength) #0.007
max(merged_dat$songlength) #167.5

### plots
# visualize the posteriors
# Q1 totsongs
mcmc_plot(nbmod1_totsongs, type = "areas", 
          prob = .9,
          point_est = "median") + 
  geom_vline(xintercept=0, linetype="dashed", alpha=0.2) +
  theme_bw()

# visualize the posteriors
# Q1 songlength
mcmc_plot(mod1_songlength, type = "areas", 
          prob = .9, 
          point_est = "median") + 
  geom_vline(xintercept=0, linetype="dashed", alpha=0.2) +
  theme_bw()

#plot the posterior - total songs ~ boat noise
posteriorQ1a <- data.frame(nbmod1_totsongs)
ggplot(posteriorQ1a, aes(x = b_scaleboatnoise)) +
  geom_histogram() +
  scale_fill_manual(values=c("#20a198"))+
  theme_bw()

#plot the posterior - song length ~ boat noise
posteriorQ1b <- data.frame(mod1_songlength)
ggplot(posteriorQ1b, aes(x = b_scaleboatnoise)) +
  geom_histogram() +
  scale_fill_manual(values=c("#20a198"))+
  theme_bw()


#plot predictions - by default all other vairables are held at their means
plot_predictions(nbmod1_totsongs, condition="boatnoise") + theme_bw()
plot_predictions(mod1_songlength, condition="boatnoise") + theme_bw()

#random effects mcmc_plot
mcmc_plot(nbmod1_totsongs, pars = "^r")


#### Interpret Results - Research Q 2 ####
#### Main effect of interest is boat ACTIVITY
summary(nbmod2_totsongs)
summary(mod2_songlength)

bayes_R2(nbmod2_totsongs) #0.58
bayes_R2(mod2_songlength) #0.34

performance_mae(nbmod2_totsongs) #23.9
performance_mae(mod2_songlength) #10.0

min(merged_dat$totsongs) #0
max(merged_dat$totsongs) #777

### plots
# visualize the posteriors
# Q1 totsongs
mcmc_plot(nbmod2_totsongs, type = "areas", 
          prob = .9,
          point_est = "median") + 
  geom_vline(xintercept=0, linetype="dashed", alpha=0.2) +
  theme_bw()

# visualize the posteriors
# Q1 songlength
mcmc_plot(mod2_songlength, type = "areas", 
          prob = .9, 
          point_est = "median") + 
  geom_vline(xintercept=0, linetype="dashed", alpha=0.2) +
  theme_bw()

#plot the posterior - total songs ~ boat noise
posteriorQ1a <- data.frame(nbmod2_totsongs)
ggplot(posteriorQ1a, aes(x = b_scaleboatactivity)) +
  geom_histogram() +
  scale_fill_manual(values=c("#20a198"))+
  theme_bw()

#plot the posterior - song length ~ boat noise
posteriorQ1b <- data.frame(mod2_songlength)
ggplot(posteriorQ1b, aes(x = b_scaleboatnoise)) +
  geom_histogram() +
  scale_fill_manual(values=c("#20a198"))+
  theme_bw()


#plot predictions - by default all other vairables are held at their means
plot_predictions(nbmod2_totsongs, condition="boatactivity") + theme_bw()
plot_predictions(mod2_songlength, condition="boatactivity") + theme_bw()

