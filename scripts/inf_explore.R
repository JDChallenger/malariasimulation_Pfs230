dinf1 <- data.frame('t' = (out1$timestep - tbv_events$timestep[1])/365, #'ibm4'=out4$infectivity,
                     'ibm1_u5'=out1$infectivity_under5, 
                    'ibm1_SAC'=out1$infectivity_SAC, 
                     'ibm1_16p'=out1$infectivity_16plus)
dinf1m <- melt(dinf1, id.vars = 't')
dinf1m$value2 <- dinf1m$value / human_population
head(dinf1m)
pl3 <- ggplot() + geom_area(data = dinf1m, aes(x=t, y=value2, fill=variable),
                            position = position_stack(reverse = T)) + theme_bw() +
  ylab('Infectivity of the human population') + labs(fill = 'Age group') + 
  xlab('Time since TB31F introduced (years)') + 
  scale_fill_manual(values = c('orange','purple','blue'), 
                    labels = c('Under 5s','S.A.C.','Adults')) #+
pl3


dinf3 <- data.frame('t' = (out3$timestep)/365, #'ibm4'=out4$infectivity,
                    'ibm_u5'=out3$infectivity_under5, 
                    'ibm_SAC'=out3$infectivity_SAC, 
                    'ibm_16p'=out3$infectivity_16plus)
dinf3m <- melt(dinf3, id.vars = 't')
dinf3m$value2 <- dinf3m$value / human_population
head(dinf3m)
ggplot() + geom_area(data = dinf3m, aes(x=t, y=value2, fill=variable),
                     position = position_stack(reverse = T)) + theme_bw() +
  ylab('Infectivity of the human population') + labs(fill = 'Age group') + 
  xlab('Time since TB31F introduced (years)') + xlim(6,9.2) + 
  scale_fill_manual(values = c('orange','purple','blue'), 
                    labels = c('Under 5s','S.A.C.','Adults')) +
  geom_vline(xintercept = tbv_params$tbv_timesteps[1]/365)
