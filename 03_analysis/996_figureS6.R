# fire characteristics in Siberia
library(tidyverse)
library(sf)


wdir = "D:/waves/"
setwd(wdir)

### specs & load data ---------------------

# load fire/weather/wave data
load(paste0('2_pipeline/99_siberia/data_larch_bor_JJA_0.25_56to74_fire.Rdata'))
load(paste0('2_pipeline/99_siberia/data_larch_bor_JJA_0.25_56to74_lightning.Rdata'))

# set plotting paths
plotpath = paste0('3_visuals/0.25deg/fig4/') # general
options(scipen = 999)
theme.size = 7


## preprocessing -------------------------

### lightning
light_reg_tot = lightdat_reg %>% group_by(week, year, reg) %>%
  summarise(lightn = sum(lightn)) %>% filter(week > 13, week < 27)
tundra_light = filter(light_reg_tot, reg == 'tundra') %>% full_join(select(fire_3cl, fire, year, week)) %>% mutate(reg = 'tundra') %>%
  rbind(filter(light_reg_tot, reg == 'boreal') %>% full_join(select(fire_3cl, fire, year, week)) %>% mutate(reg = 'boreal')) %>%
  mutate(lightn = ifelse(is.na(lightn), 0, lightn))
tundra_light_anom = group_by(tundra_light, week, reg) %>% summarise(meanlight = mean(lightn), sdlight = sd(lightn)) %>%
  full_join(tundra_light) %>% mutate(lightanom = (lightn - meanlight)/sdlight)

### fire
burn_tot = firedat %>% group_by(week, year, reg) %>% 
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to turn area from % in km2
  summarise(burned = sum(perc_burned))
tundra_burn = filter(burn_tot, reg == 2) %>% full_join(select(fire_3cl, fire, year, week)) %>% mutate(reg = 2) %>%
  rbind(filter(burn_tot, reg == 255) %>% full_join(select(fire_3cl, fire, year, week))) %>%
  mutate(burned = ifelse(is.na(burned), 0, burned), reg = ifelse(reg == 255, 'boreal', 'tundra'))
tundra_burn_anom = group_by(tundra_burn, reg, week) %>% summarise(meanburn = mean(burned), sdburn = sd(burned)) %>%
  full_join(tundra_burn) %>% mutate(burnanom = (burned - meanburn)/sdburn)

# both with only 2 classes: extreme weeks/all other weeks
tundra_light_2cl = tundra_light %>% 
  mutate(fire2 = ifelse(fire == 'more', 'more', 'other')) %>%
  group_by(fire2, reg, week, year) %>% summarise(lightn = sum(lightn))
tundra_burn_2cl = tundra_burn %>%
  mutate(fire2 = ifelse(fire == 'more', 'more', 'other')) %>%
  group_by(fire2, reg, week, year) %>% summarise(burned = sum(burned))



## Fig. S6A&B ----------------------------------

ggplot() + theme_bw() + 
  labs(y = 'Fraction of data', x = 'Weekly lightning strikes') +
  geom_density(data = filter(tundra_light_2cl, reg == 'boreal', fire2 == 'other'), 
               aes(x = lightn, y = ..count../sum(..count..)), 
               alpha = 0.5, fill = 'lightgray', col = 'black', lty = 2) +
  geom_density(data = filter(tundra_light_2cl, reg == 'boreal', fire2 == 'more'), 
               aes(x = lightn, y = ..count../sum(..count..)), 
               alpha = 0.5, fill = 'firebrick', col = NA) +
  scale_x_continuous(trans = 'pseudo_log', breaks = 10^(0:5)) + #,
  # labels = trans_format("log10", math_format(10^.x))) + 
  geom_text(aes(x = 0, y = 0.0065), label = 'A', size = 9/.pt, fontface = "bold") +
  theme(text = element_text(size = 7))
ggsave(paste0(plotpath, 'fig4_a.svg'), width = 6, height = 5.8, dpi = 300, units = 'cm')

ggplot() + theme_bw() + 
  labs(y = '', x = 'Weekly lightning strikes') +
  geom_density(data = filter(tundra_light_2cl, reg == 'tundra', fire2 == 'other'), 
               aes(x = lightn, y = ..count../sum(..count..)), #), 
               alpha = 0.5, fill = 'lightgray', col = 'black', lty = 2) +
  geom_density(data = filter(tundra_light_2cl, reg == 'tundra', fire2 == 'more'), 
               aes(x = lightn, y = ..count../sum(..count..)), #), 
               alpha = 0.5, fill = 'firebrick', col = NA) +
  scale_x_continuous(trans = 'pseudo_log', breaks = 10^(0:6)) + #coord_flip() +
  geom_text(aes(x = 0, y = 0.012), label = 'B', size = 9/.pt, fontface = "bold") +
  theme(text = element_text(size = 7))
ggsave(paste0(plotpath, 'fig4_b.svg'), width = 6, height = 5.8, dpi = 300, units = 'cm')


### Fig. S6C: Table comparing directly to Arctic front jet and snowmelt --------------
light_reg_tot = light_reg_tot %>% 
  full_join(data.frame(year = rep(2012:2020, each = 26),
                       week = rep(rep(14:26, 9), each = 2),
                       reg = rep(c('tundra', 'boreal'), 13*9))) %>%
  mutate(lightn = ifelse(is.na(lightn), 0, lightn)) 
test = light_reg_tot %>% group_by(week, reg) %>% 
  summarise(light_mean = mean(lightn), light_sd = sd(lightn)) %>%
  full_join(light_reg_tot) %>% mutate(light_anom = (lightn-light_mean)/light_sd) %>%
  left_join(mean_sm) %>% 
  left_join(V250_highfire_cor)

# test for normality: non-normal!
ks.test(test$light_anom, 'pnorm')
ks.test(test$snow, 'pnorm')
ks.test(test$corV, 'pnorm')
shapiro.test(test$light_anom)
shapiro.test(test$snow)
shapiro.test(test$corV)

cormethod = 'spearman'

# overall correlations
test %>% group_by(reg) %>%
  summarise(snow_cor = cor.test(light_anom, snow, method = cormethod)$estimate,
            snow_p = cor.test(light_anom, snow, method = cormethod)$p.value,
            V_cor = cor.test(light_anom, corV, method = cormethod)$estimate,
            V_p = cor.test(light_anom, corV, method = cormethod)$p.value)

# monthly correlations
test %>% mutate(mon = ifelse(week %in% 14:17, 'jun', 
                             ifelse(week %in% 18:22, 'jul', 'aug'))) %>%
  group_by(reg, mon) %>%
  summarise(snow_cor = cor.test(light_anom, snow, method = cormethod)$estimate,
            snow_p = cor.test(light_anom, snow, method = cormethod)$p.value,
            V_cor = cor.test(light_anom, corV, method = cormethod)$estimate,
            V_p = cor.test(light_anom, corV, method = cormethod)$p.value)

# end-june
test %>% filter(week %in% 16:17) %>% group_by(reg) %>%
  summarise(snow_cor = cor.test(light_anom, snow, method = cormethod)$estimate,
            snow_p = cor.test(light_anom, snow, method = cormethod)$p.value)

# weekly timeseries of correlation AFJ - lightning
when = test %>% group_by(week, reg) %>% 
  summarise(cor = cor.test(lightn, corV)$estimate,
            p = cor.test(lightn, corV)$p.value)

ggplot(when, aes(week, cor)) + theme_bw() + facet_wrap(~reg) +
  geom_hline(yintercept = 0, lty = 2, col = 'gray') +
  geom_point() + geom_line() + geom_point(data = filter(when, p < 0.1), aes(week, cor), col = 'red') +
  labs(x = 'Week since March 1', y = 'Cor (No lightning strikes ~ AFJ pattern cor)') +
  theme(text = element_text(size = theme.size))
ggsave(paste0(plotpath, 'ts_lightn_V250_cor.png'), 
       width = 9, height = 5.5, dpi = 300, units = 'cm')

# weekly timeseries of correlation SNOW - lightning
when = test %>% group_by(week, reg) %>% 
  summarise(cor = cor.test(lightn, snow)$estimate,
            p = cor.test(lightn, snow)$p.value)

ggplot(when, aes(week, cor)) + theme_bw() + facet_wrap(~reg) +
  geom_hline(yintercept = 0, lty = 2, col = 'gray') +
  geom_point() + geom_line() + geom_point(data = filter(when, p < 0.1), aes(week, cor), col = 'red') +
  labs(x = 'Week since March 1', y = 'Cor (No lightning strikes ~ snow anomaly)') +
  theme(text = element_text(size = theme.size))
ggsave(paste0(plotpath, 'ts_lightn_snow_cor.png'), 
       width = 9, height = 5.5, dpi = 300, units = 'cm')


# snow lightning cor June-July only
ungroup(test) %>% filter(week < 23) %>% summarise(cor = cor.test(light_anom, snow)$estimate,
                                         p = cor.test(light_anom, snow)$p.value)
ungroup(test) %>% mutate(mon = ifelse(week %in% 14:17, 'jun',
                               ifelse(week %in% 18:22, 'jul', 'aug'))) %>%
  ggplot() + theme_bw() + facet_wrap(reg~mon, scales = 'free_y') +
  geom_point(aes(snowanom, lightn, col = mon))


