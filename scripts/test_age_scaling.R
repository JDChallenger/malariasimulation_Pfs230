#suppressPackageStartupMessages(library(ggplot2))
library(malariasimulation)
library(malariaEquilibrium)
library(cowplot)
library(reshape2)
library(ggplot2)
library(gridExtra)

year <- 365
month <- 30 # ish
sim_length <- 10 * year
human_population <- 19000
starting_EIR <- 70

#One way to link seasonality to admin units- using input from deterministic package
#seas <- read.csv('C:\\Users\\jchallen\\Dropbox\\NC_data_code\\deterministic-malaria-model-master-NC\\data-raw\\admin_units_seasonal.csv')
#lst <- as.numeric(seas[seas$NAME_0=='Burkina Faso' & seas$NAME_1=='Sahel',5:11])

seas <- read.csv('/Users/jchallen/Documents/Git/deterministic-malaria-model/data-raw/admin_units_seasonal.csv')
lst <- as.numeric(seas[seas$NAME_0=='Burkina Faso' & seas$NAME_1=='Cascades',5:11])

simparams <- get_parameters(list(
  human_population = human_population,
  model_seasonality = T, 
  individual_mosquitoes = F,
  phi_bednets = 0.79, #proportion of bites taken in bed. Check with Ellie's review!
  species = 1,
  g0 = lst[1],
  g = c(lst[2], lst[4], lst[6]), 
  h = c(lst[3], lst[5], lst[7]),
  #tbv_md = 0.00028045,
  #tbv_mt = 0.00028045,
  #tbv_ma = 0.00028045,
  #tbv_mu = 0.00028045,
  prevalence_rendering_min_ages = 0 * year,
  prevalence_rendering_max_ages = 100 * year, 
  clinical_incidence_rendering_min_ages = 0,
  clinical_incidence_rendering_max_ages = 100 * year,
  incidence_rendering_min_ages = 0,
  incidence_rendering_max_ages = 100 * year
)
)

simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params))
simparams <- set_clinical_treatment(simparams, 1, 1, 0.5 )

simparams <- set_equilibrium(simparams, starting_EIR)

peak <- peak_season_offset(simparams)
peak

out1 <- run_simulation(sim_length, simparams)

ggplot(out1) + geom_line(aes(x = timestep, y = n_detect_lm_0_36500/ n_age_0_36500)) + 
  theme_classic()

######################### 2. PEV only (R21) ##########################

# Mass PEV
R21params <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, #
  max_ages = 63 * year, #
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 1,
  booster_spacing = c(1*year), 
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),
  booster_profile = list(r21_booster_profile)
)

out2 <- run_simulation(sim_length, R21params)

R21params_A <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 0.6,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2A <- run_simulation(sim_length, R21params_A)

R21params_B <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 0.4,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2B <- run_simulation(sim_length, R21params_B)

R21params_C <- set_mass_pev(
  simparams,
  profile = r21_profile,
  coverages = c(0.9),
  min_wait = 0,
  min_ages = 5*year, # The minimum age for the target population to be vaccinated. 
  max_ages = 63 * year, # The maximum age for the target population to be vaccinated.
  timesteps = c(7)*year + peak - 3.5*month,
  adult_scaling = 0.2,
  booster_spacing = c(1*year), # The booster is given 10 years after?
  booster_coverage = matrix(rep(0.95,1), nrow = 1, ncol= 1),#, # Coverage of the booster dose is 5%.
  booster_profile = list(r21_booster_profile)
)

out2C <- run_simulation(sim_length, R21params_C)

df2 <- data.frame('t' = out1$timestep/365, 
                  'pr1'=out1$n_detect_lm_0_36500/out1$n_age_0_36500,
                  'pr2'=out2$n_detect_lm_0_36500/out2$n_age_0_36500,
                  'pr2A'=out2A$n_detect_lm_0_36500/out2A$n_age_0_36500,
                  'pr2B'=out2B$n_detect_lm_0_36500/out2B$n_age_0_36500,
                  'pr2C'=out2C$n_detect_lm_0_36500/out2C$n_age_0_36500)
ggplot(df2) + geom_line(aes(x = t, y = pr1)) + theme_classic() + 
  geom_line(aes(x = t, y = pr2), color = 'purple') + xlim(c(7,10)) + 
  geom_line(aes(x = t, y = pr2A), color = 'magenta') +
  geom_line(aes(x = t, y = pr2B), color = 'orange') +
  geom_line(aes(x = t, y = pr2C), color = 'darkgreen') 

df2i <- data.frame('t' = out1$timestep/365, 
                  'inc1' = out1$n_inc_clinical_0_36500/out1$n_age_0_36500,
                  'inc2' = out2$n_inc_clinical_0_36500/out2$n_age_0_36500,
                  'inc2A' = out2A$n_inc_clinical_0_36500/out2A$n_age_0_36500,
                  'inc2B' = out2B$n_inc_clinical_0_36500/out2B$n_age_0_36500,
                  'inc2C' = out2C$n_inc_clinical_0_36500/out2C$n_age_0_36500)

ggplot(df2i) + geom_line(aes(x = t, y = inc1)) + theme_classic() + 
  geom_line(aes(x = t, y = inc2), color = 'purple') + xlim(c(7,10)) + 
  geom_line(aes(x = t, y = inc2A), color = 'magenta') +
  geom_line(aes(x = t, y = inc2B), color = 'orange') +
  geom_line(aes(x = t, y = inc2C), color = 'darkgreen')

#cumulative inc
tt <- R21params$mass_pev_timesteps + 90

ci <- out1$n_inc_clinical_0_36500[(tt):sim_length]
ci2 <- out2$n_inc_clinical_0_36500[(tt):sim_length]
ci2A <- out2A$n_inc_clinical_0_36500[(tt):sim_length]
ci2B <- out2B$n_inc_clinical_0_36500[(tt):sim_length]
ci2C <- out2C$n_inc_clinical_0_36500[(tt):sim_length]

for(i in 2:(length(ci))){
  ci[i] <- ci[i] + ci[i-1]
  ci2[i] <- ci2[i] + ci2[i-1]
  ci2A[i] <- ci2A[i] + ci2A[i-1]
  ci2B[i] <- ci2B[i] + ci2B[i-1]
  ci2C[i] <- ci2C[i] + ci2C[i-1]
}
gh <- data.frame('t' = seq((tt):sim_length)/365, 
                 'ci' = ci, 'ci2' = ci2, 'ci2A' = ci2A, 'ci2B' = ci2B, 'ci2C' = ci2C)#, 'prop' = (1-ci3/ci2))
ggplot(gh) + #geom_line(aes(x=t, y=ci)) + 
  geom_line(aes(x=t, y=ci2), color = 'purple') + theme_classic() + 
  geom_line(aes(x=t, y=ci2A),color = 'magenta') +  
  geom_line(aes(x=t, y=ci2B),color = 'orange') +
  geom_line(aes(x=t, y=ci2C),color = 'darkgreen') 


######################### 3. TBV only ##########################

tbv_params <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  ages = seq(1,99), #not c(1,99)
  adult_scaling = 1
)

out3 <- run_simulation(sim_length, tbv_params)
table(out3$n_vaccinated_tbv) #???

tbv_paramsA <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.7,
  ages = seq(1,99)
)

out3A <- run_simulation(sim_length, tbv_paramsA)

tbv_paramsB <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.4,
  ages = seq(1,99)
)

out3B <- run_simulation(sim_length, tbv_paramsB)

tbv_paramsC <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.2,
  ages = seq(1,99)
)

out3C <- run_simulation(sim_length, tbv_paramsC)

df3 <- data.frame('t' = out1$timestep/365, 
                  'pr1'=out1$n_detect_lm_0_36500/out1$n_age_0_36500,
                  'pr3'=out3$n_detect_lm_0_36500/out3$n_age_0_36500,
                  'pr3A'=out3A$n_detect_lm_0_36500/out3A$n_age_0_36500,
                  'pr3B'=out3B$n_detect_lm_0_36500/out3B$n_age_0_36500,
                  'pr3C'=out3C$n_detect_lm_0_36500/out3C$n_age_0_36500)
ggplot(df3) + geom_line(aes(x = t, y = pr1)) + theme_classic() + 
  geom_line(aes(x = t, y = pr3), color = 'purple') + xlim(c(7,10)) + 
  geom_line(aes(x = t, y = pr3A), color = 'magenta') +
  geom_line(aes(x = t, y = pr3B), color = 'orange') +
  geom_line(aes(x = t, y = pr3C), color = 'darkgreen') 

df3i <- data.frame('t' = out1$timestep/365, 
                   'inc1' = out1$n_inc_clinical_0_36500/out1$n_age_0_36500,
                   'inc3' = out3$n_inc_clinical_0_36500/out3$n_age_0_36500,
                   'inc3A' = out3A$n_inc_clinical_0_36500/out3A$n_age_0_36500,
                   'inc3B' = out3B$n_inc_clinical_0_36500/out3B$n_age_0_36500,
                   'inc3C' = out3C$n_inc_clinical_0_36500/out3C$n_age_0_36500)

ggplot(df3i) + geom_line(aes(x = t, y = inc1)) + theme_classic() + 
  geom_line(aes(x = t, y = inc3), color = 'purple') + xlim(c(7,10)) + 
  geom_line(aes(x = t, y = inc3A), color = 'magenta') +
  geom_line(aes(x = t, y = inc3B), color = 'orange') +
  geom_line(aes(x = t, y = inc3C), color = 'darkgreen')

#cumulative inc
tt <- tbv_params$tbv_timesteps[1]

ci <- out1$n_inc_clinical_0_36500[(tt):sim_length]
ci3 <- out3$n_inc_clinical_0_36500[(tt):sim_length]
ci3A <- out3A$n_inc_clinical_0_36500[(tt):sim_length]
ci3B <- out3B$n_inc_clinical_0_36500[(tt):sim_length]
ci3C <- out3C$n_inc_clinical_0_36500[(tt):sim_length]

for(i in 2:(length(ci))){
  ci[i] <- ci[i] + ci[i-1]
  ci3[i] <- ci3[i] + ci3[i-1]
  ci3A[i] <- ci3A[i] + ci3A[i-1]
  ci3B[i] <- ci3B[i] + ci3B[i-1]
  ci3C[i] <- ci3C[i] + ci3C[i-1]
}
gh3 <- data.frame('t' = seq((tt):sim_length)/365, 
                 'ci' = ci, 'ci3' = ci3, 'ci3A' = ci3A, 'ci3B' = ci3B, 'ci3C' = ci3C)#, 'prop' = (1-ci3/ci2))
ggplot(gh3) + geom_line(aes(x=t, y=ci)) + 
  geom_line(aes(x=t, y=ci3), color = 'purple') + theme_classic() + 
  geom_line(aes(x=t, y=ci3A),color = 'magenta') +  ylim(c(0,33000)) +
  geom_line(aes(x=t, y=ci3B),color = 'orange') +
  geom_line(aes(x=t, y=ci3C),color = 'darkgreen') + xlim(c(0,2))

#################### TBV, but on non-scaled age groups


tbv_params <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 1,
  ages = seq(0,15)
)

out4 <- run_simulation(sim_length, tbv_params)

tbv_paramsA <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.7,
  ages = seq(0,15)
)

out4A <- run_simulation(sim_length, tbv_paramsA)

tbv_paramsB <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.4,
  ages = c(0,15)
)

out4B <- run_simulation(sim_length, tbv_paramsB)

tbv_paramsC <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.2,
  ages = seq(0,15)
)

out4C <- run_simulation(sim_length, tbv_paramsC)


#cumulative inc (test)
tt <- tbv_params$tbv_timesteps[1]

ci <- out1$n_inc_clinical_0_36500[(tt):sim_length]
ci4 <- out4$n_inc_clinical_0_36500[(tt):sim_length]
ci4A <- out4A$n_inc_clinical_0_36500[(tt):sim_length]
ci4B <- out4B$n_inc_clinical_0_36500[(tt):sim_length]
ci4C <- out4C$n_inc_clinical_0_36500[(tt):sim_length]

for(i in 2:(length(ci))){
  ci[i] <- ci[i] + ci[i-1]
  ci4[i] <- ci4[i] + ci4[i-1]
  ci4A[i] <- ci4A[i] + ci4A[i-1]
  ci4B[i] <- ci4B[i] + ci4B[i-1]
  ci4C[i] <- ci4C[i] + ci4C[i-1]
}
gh4 <- data.frame('t' = seq((tt):sim_length)/365, 
                  'ci' = ci, 'ci4' = ci4, 'ci4A' = ci4A, 'ci4B' = ci4B, 'ci4C' = ci4C)#, 'prop' = (1-ci3/ci2))
ggplot(gh4) + geom_line(aes(x=t, y=ci)) + 
  geom_line(aes(x=t, y=ci4), color = 'purple') + theme_classic() + 
  geom_line(aes(x=t, y=ci4A),color = 'magenta') +  ylim(c(0,33000)) +
  geom_line(aes(x=t, y=ci4B),color = 'orange') +
  geom_line(aes(x=t, y=ci4C),color = 'darkgreen') + xlim(c(0,2))

#################### TBV, but ONLY on scaled age groups

tbv_params <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 1,
  ages = seq(16,99)
)

out5 <- run_simulation(sim_length, tbv_params)

tbv_paramsA <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.7,
  ages = seq(16,99)
)

out5A <- run_simulation(sim_length, tbv_paramsA)

tbv_paramsB <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.4,
  ages = seq(16,99)
)

out5B <- run_simulation(sim_length, tbv_paramsB)

tbv_paramsC <- set_tbv(
  simparams,
  timesteps = c(7,8)*year + peak - 2*month,
  coverages = rep(0.9, 2),
  adult_scaling = 0.2,
  ages = c(16,99)
)

out5C <- run_simulation(sim_length, tbv_paramsC)


#cumulative inc
tt <- tbv_params$tbv_timesteps[1]

ci <- out1$n_inc_clinical_0_36500[(tt):sim_length]
ci5 <- out5$n_inc_clinical_0_36500[(tt):sim_length]
ci5A <- out5A$n_inc_clinical_0_36500[(tt):sim_length]
ci5B <- out5B$n_inc_clinical_0_36500[(tt):sim_length]
ci5C <- out5C$n_inc_clinical_0_36500[(tt):sim_length]

for(i in 2:(length(ci))){
  ci[i] <- ci[i] + ci[i-1]
  ci5[i] <- ci5[i] + ci5[i-1]
  ci5A[i] <- ci5A[i] + ci5A[i-1]
  ci5B[i] <- ci5B[i] + ci5B[i-1]
  ci5C[i] <- ci5C[i] + ci5C[i-1]
}
gh5 <- data.frame('t' = seq((tt):sim_length)/365, 
                  'ci' = ci, 'ci5' = ci5, 'ci5A' = ci5A, 'ci5B' = ci5B, 'ci5C' = ci5C)#, 'prop' = (1-ci3/ci2))
ggplot(gh5) + geom_line(aes(x=t, y=ci)) + 
  geom_line(aes(x=t, y=ci5), color = 'purple') + theme_classic() + 
  geom_line(aes(x=t, y=ci5A),color = 'magenta') +  ylim(c(0,35000)) +
  geom_line(aes(x=t, y=ci5B),color = 'orange') +
  geom_line(aes(x=t, y=ci5C),color = 'darkgreen') + xlim(c(0,2))

