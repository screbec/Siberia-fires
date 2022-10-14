library(tidyverse)
library(sf)
library(trend) # for sens slope
library(viridis)
library(ggpubr)
library(lme4)
library(lubridate)
library(caret)

wdir = "D:/waves/"
setwd(wdir)

theme.size = 7 # text size for plotting
geom.text.size = theme.size / (14/5)
plotpath = paste0('3_visuals/0.25deg/')

### load data and process ------------------------------------

# load the wave and fire data
load(paste0('2_pipeline/99_siberia/data_larch_bor_JJA_0.25_56to74_fire.Rdata'))

### load modis snow melt files
snowmelt = readRDS(file = '2_pipeline/05_snowmelt/sm_gee.RDS') %>% rename(lon = x, lat = y)
snowmelt_025deg = readRDS(file = '2_pipeline/05_snowmelt/sm_gee_025deg.RDS') %>% rename(lon = x, lat = y)
mean_sm = readRDS(file = '2_pipeline/05_snowmelt/sm_gee_mean.RDS') %>% add_column(type = 'MODIS')
sm_full = readRDS(file = '2_pipeline/05_snowmelt/sm_nsidc_larch_full.RDS') %>%
  dplyr::select(c(lat, lon, year, snowmelt)) %>%
  mutate(lat = floor(lat / 0.25) * 0.25 + (0.25/2), lon = floor(lon / 0.25) * 0.25 + (0.25/2))

# load nsidc sm 
mean_sm_full = readRDS(file = '2_pipeline/05_snowmelt/sm_nsidc_mean.RDS') %>%
  add_column(type = 'NSIDC')# yearly, full timeseries (always NSIDC)

# load shapefile of study region
boreal = st_read('2_pipeline/00_ecoregions/eastern_siberia.shp')

# load pattern correlation (full time series, 1979-2021)
V250_highfire_cor = readRDS('2_pipeline/99_siberia/0.25deg/V250_highfire_cor.rds')

# load tree line data
tl = read_csv('2_pipeline/01_modis/distance_treeline_centroids.csv') %>% select(-X1) # treeline data without tundra differentiation

## prep data for compound plot ---------------------

# merge data with treeline distance
firedat = firedat %>% left_join(tl) %>% mutate(dist_tl = ifelse(is.na(dist_tl), 0, dist_tl),
                                               dist_tl = ifelse(dist_tl == 0 & lat < 60, NA, dist_tl))
snowmelt_025deg = snowmelt_025deg %>% left_join(tl) %>%
  mutate(dist_tl = ifelse(is.na(dist_tl), 0, dist_tl),dist_tl = ifelse(dist_tl == 0 & lat < 60, NA, dist_tl))
sm_full = sm_full %>% left_join(tl) %>% mutate(dist_tl = ifelse(is.na(dist_tl), 0, dist_tl),
                                               dist_tl = ifelse(dist_tl == 0 & lat < 60, NA, dist_tl))

# merge snowmelt, V and fire
plotdata = V250_highfire_cor %>%
  # join with burned area anomaly
  right_join(fire_anom) %>%
  # join with snowmelt day anomaly
  left_join(mean_sm)


### Fig 3A: Weekly compound plot -----------------------

sm_19_20_21 = filter(mean_sm, year %in% 2019:2021)$snowanom_days
max_V_19_20_21 = (filter(plotdata, year %in% 2019:2021) %>% group_by(year) %>% summarise(corV = max(corV)))$corV
min_V_19_20_21 = (filter(plotdata, year %in% 2019:2021) %>% group_by(year) %>% summarise(corV = min(corV)))$corV

## compound plot
ggplot(plotdata) + theme_bw() + 
  labs(x = 'Snowmelt anomaly (days)', y = 'Correlation with Arctic front jet', col = 'Burn anomaly') +
  geom_vline(xintercept = 0, col = 'gray', lty = 2, size = 0.3) + 
  geom_hline(yintercept = 0, col = 'gray', lty = 2, size = 0.3) +
  geom_point(aes(x = snowanom_days, y = corV, col = burnanom), size = 0.3) +
  scale_color_steps2(low = 'dodgerblue', high = 'darkred', mid = 'gray90', 
                     n.breaks = 8, midpoint = 0, limits = c(-1.18, 4.37),
                     guide=guide_colourbar(title.position="top", 
                                           title.theme = element_text( size = theme.size))) +
  geom_text(aes(x = -6, y = 0.52), label = 'A', size = 10/.pt, fontface = "bold") + 
  theme(text = element_text(size = theme.size), plot.margin = unit(c(0,1,1,1), "mm"), aspect.ratio = 9/10,
        legend.position="top", legend.title.align=0.5, legend.margin=margin(0,0,0,0),
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.05, "cm"),
        legend.box.margin=margin(0,-5,-5,-5)) + 
  geom_text(data = filter(plotdata, year %in% 2019:2021, week == 14), 
            aes(x = snowanom_days, y = max_V_19_20_21+0.06, label = year), size = 6/.pt, col = 'gray22') +
  annotate("rect", xmin = sm_19_20_21-0.14, xmax = sm_19_20_21+0.14, ymin = min_V_19_20_21-0.02, ymax = max_V_19_20_21+0.02, 
           fill = NA, col = 'gray22', size = 0.3/.pt)

ggsave(paste0(plotpath, 'compound/snow_fire_V250_all_modis_', firelatlim, '-74.svg'), 
       width = 5.5, height = 6, dpi = 300, units = 'cm')


## quadrants probabilities
burnthresh = 1
dim(filter(plotdata, burnanom > burnthresh))[1]/dim(plotdata)[1]
# upper left (compound effect)
dim(filter(plotdata, burnanom > burnthresh, snowanom < 0, corV > 0))[1]/
  dim(filter(plotdata, snowanom < 0, corV > mean(plotdata$corV)))[1]
# lower left (snow only)
dim(filter(plotdata, burnanom > burnthresh, snowanom < 0, corV < 0))[1]/
  dim(filter(plotdata, snowanom < 0, corV < mean(plotdata$corV)))[1]
# upper right (corV only)
dim(filter(plotdata, burnanom > burnthresh, snowanom > 0, corV > 0))[1]/
  dim(filter(plotdata, snowanom > 0, corV > mean(plotdata$corV)))[1]
# lwoer right (none)
dim(filter(plotdata, burnanom > burnthresh, snowanom > 0, corV < 0))[1]/
  dim(filter(plotdata, snowanom > 0, corV < mean(plotdata$corV)))[1]
# late snowmelt
dim(filter(plotdata, burnanom > burnthresh, snowanom > 0))[1]/
  dim(filter(plotdata, snowanom > 0))[1]
dim(filter(plotdata, burnanom > burnthresh, snowanom < 0))[1]/
  dim(filter(plotdata, snowanom < 0))[1]
# no arctic front jet
dim(filter(plotdata, burnanom > burnthresh, corV < 0))[1]/
  dim(filter(plotdata, corV < 0))[1]

## linear model with crossvalidation
library(boot)
# plotdata_scale = plotdata %>% mutate(corV = scale(corV))
data_ctrl <- trainControl(method = "boot", number = 100)
model_caret <- train(burnanom ~ snowanom + corV, data = plotdata,                        
                     trControl = data_ctrl, method = "lm", na.action = na.pass) 
model_caret
model_caret$finalModel
sd(model_caret$resample$Rsquared)

model_coef <- function(data, index){
  coef(lm(burnanom ~ snowanom + corV, data = data, subset = index))
}
boot(plotdata, model_coef, 100)


### Fig 3B: Weekly compound plot tundra ---------------------------------------------------

burn_tot = firedat %>% filter(reg == 2) %>% group_by(week, year) %>% 
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to turn area from % in km2
  summarise(burned = sum(perc_burned)) %>% filter(week > 13, week < 27)
tundra_burn = burn_tot %>% full_join(select(fire_3cl, fire, year, week)) %>%
  mutate(burned = ifelse(is.na(burned), 0, burned))
tundra_burn_anom = group_by(tundra_burn, week) %>% summarise(meanburn = mean(burned), sdburn = sd(burned)) %>%
  full_join(tundra_burn) %>% mutate(burnanom = (burned - meanburn)/sdburn)

plotdata_tundra = V250_highfire_cor %>%
  # join with burned area anomaly
  inner_join(tundra_burn) %>%
  # join with snowmelt day anomaly
  left_join(mean_sm)
plotdata_tundra = plotdata_tundra %>% group_by(week) %>%
  summarise(meanburn =mean(burned), sdburn = sd(burned)) %>%
  full_join(plotdata_tundra) %>% mutate(burnanom = (burned-meanburn)/sdburn)

# compound plot with tundra burning
ggplot(plotdata_tundra) + theme_bw() +
  labs(x = 'Snowmelt anomaly (days)', y = 'Correlation with Arctic front jet', col = 'Burn anomaly') +
  geom_vline(xintercept = 0, col = 'gray', lty = 2) + 
  geom_hline(yintercept = 0, col = 'gray', lty = 2) +
  geom_point(aes(x = snowanom_days, y = corV, col = burnanom), size = 0.3) +
  scale_color_steps2(low = 'dodgerblue', high = 'darkred', mid = 'gray90', 
                     n.breaks = 8, midpoint = 0, limits = c(-1.18, 4.37),
                     guide=guide_colourbar(title.position="top", 
                                           title.theme = element_text( size = theme.size))) +
  geom_text(aes(x = -6, y = 0.52), label = 'B', size = 10/.pt, fontface = "bold") + 
  theme(text = element_text(size = theme.size), plot.margin = unit(c(0,1,1,1), "mm"), aspect.ratio = 9/10,
        legend.position="top", legend.title.align=0.5, legend.margin=margin(0,0,0,0),
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.05, "cm"),
        legend.box.margin=margin(0,-5,-5,-5)) + 
  geom_text(data = filter(plotdata_tundra, year %in% 2019:2021, week == 14), 
            aes(x = snowanom_days, y = max_V_19_20_21+0.06, label = year), size = 6/.pt, col = 'gray22') +
  annotate("rect", xmin = sm_19_20_21-0.14, xmax = sm_19_20_21+0.14, ymin = min_V_19_20_21-0.02, ymax = max_V_19_20_21+0.02, 
           fill = NA, col = 'gray22', size = 0.3/.pt)

ggsave(paste0(plotpath, 'compound/snow_V250_tundra_all_MODIS_rel.svg'), 
       width = 5.5, height = 6, dpi = 300, units = 'cm')


## Fig 3C: Latitudinal plot of snowmelt vs AFJ influence on fire activity -------------------------

# summarise fire data per latitude
latsum_fire = mutate(firedat, group = ceiling(dist_tl/100000)*100) %>% 
  group_by(group, week, year) %>% 
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to turn area from % in km2
  summarise(burned = sum(perc_burned), n = n())

# same for snowmelt
latsum_snowmelt = mutate(snowmelt_025deg, group = ceiling(dist_tl/100000)*100) %>% 
  group_by(group, year) %>% summarise(snowanom_px = sum(snowanom),
                                      snowsd = sd(snow), snow = mean(snow))

# regional snow climatology
latsum_snowclim = group_by(latsum_snowmelt, group) %>% summarise(meanclim = mean(snow), sdclim = sd(snow)) %>%
  right_join(latsum_snowmelt) %>%
  mutate(snowanom = (snow-meanclim)/sdclim)

# join with snowmelt and corV
latsum_join = V250_highfire_cor %>%
  # join with snowmelt day anomaly
  right_join(select(latsum_snowclim, group, year, snowanom)) %>%
  # join with burned area anomaly
  left_join(select(latsum_fire, group, week, year, burned)) %>%
  # where burnanom is NA it is really 0
  mutate(burned = ifelse(is.na(burned), 0, burned)) %>%
  # join with regional snow anomaly
  left_join(rename(mean_sm, snowanom_reg = snowanom) %>% select(year, snowanom_reg), by = 'year')

# compute burn anomaly
latsum_fireclim = group_by(latsum_join, group, week) %>% summarise(meanclim = mean(burned), sdclim = sd(burned)) %>%
  right_join(latsum_join) %>% mutate(burnanom = (burned-meanclim)/sdclim) %>%
  filter(!is.na(burnanom)) # filter out latitudes with no burning in any week

# correlations
latsum_cor = latsum_fireclim %>% group_by(group) %>%
  summarise(sm_burn = cor(snowanom, burnanom), sm_burn_p = cor.test(snowanom, burnanom)$p.value,
            V_burn = cor(corV, burnanom), V_burn_p = cor.test(corV, burnanom)$p.value, n = n()) %>% 
  mutate(V_burn_sig = ifelse(V_burn_p < 0.05, V_burn_p, NA), sm_burn_sig = ifelse(sm_burn_p < 0.05, sm_burn_p, NA))
latsum_cor$V_burn_adj = p.adjust(latsum_cor$V_burn_sig, "BH")
latsum_cor$sm_burn_sig = p.adjust(latsum_cor$sm_burn_sig, "BH")

# plot
p1 = ggplot(latsum_cor) + theme_bw() + 
  theme(text = element_text(size = 7), plot.margin = unit(c(0,1,1,1), "mm"), aspect.ratio = 9/10) +
  geom_hline(yintercept = 0, col = 'darkgray', lty = 2) +
  # plot trends/lines
  geom_smooth(aes(group, sm_burn, weight = n), method = 'lm', se = F, col = '#2171b5', size = 0.3) +
  geom_smooth(aes(group, V_burn, weight = n), method = 'loess', span = 5, se = F, col = 'darkred', size = 0.3) +
  # plot the data points
  geom_point(aes(group, sm_burn), col = '#9ecae1', size = 0.3) +
  geom_point(aes(group, V_burn), col = 'mistyrose2', size = 0.3) +
  # plot dark for significant datapoints
  geom_point(data = filter(latsum_cor, sm_burn_p < 0.05), aes(group, sm_burn), col = '#2171b5', size = 0.5) +
  geom_point(data = filter(latsum_cor, V_burn_p < 0.05), aes(group, V_burn), col = 'darkred', size = 0.5) +
  geom_text(aes(x = -1, y = -0.45), label = 'C', size = 10/.pt, fontface = "bold") +
    scale_x_reverse(breaks = rev(seq(0, 2000, by = 500)), expand = c(0.05,0)) +
    labs(x = 'Distance from tree line (km)', y = 'Correlation') + coord_flip() 
ggsave(paste0(plotpath, 'compound/lat_snow_vs_V_treeline.svg'), width = 5.5, height = 5, dpi = 300, units = 'cm')


### Fig 4A Snowmelt trend (NSIDC and MODIS) ----------------
ggplot(mean_sm_full, aes(year, snowmelt)) + theme_bw() + 
  geom_point(aes(shape = type), col = '#2171b5', size = 0.5) + 
  geom_smooth(method = 'lm', col = '#2171b5', fill = '#2171b5', size = 0.5) +
  geom_point(data = mean_sm, aes(year, snow, shape = type), col = '#2171b5', size = 0.5) + 
  scale_shape_manual(values = c(2,16)) +
  labs(x = '', y = 'Regional snowmelt day') +
  theme(text = element_text(size = 7), plot.margin = unit(c(0,1,0,0), "mm"), aspect.ratio = 1,
        legend.margin=margin(0,0,0,0), legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.05, "cm"),
        legend.position = c(0.15, 0.1), legend.title = element_blank())
ggsave(paste0(plotpath, 'snowmelt/snow_trend_NSIDC.svg'), width = 5.5, height = 5.5, dpi = 300, units = 'cm')
summary(lm(snowmelt ~ year, mean_sm_full))
summary(lm(snow ~ year, mean_sm))

# agreement between both datasets
inner_join(mean_sm, mean_sm_full, by = 'year') %>% ungroup() %>%
  summarise(cor.test(snow, snowmelt)$estimate, cor.test(snow, snowmelt)$p.value)

### Fig. 4B: Snowmelt trends per treeline band ----------------------------------------
sm_trend = mutate(sm_full, group = ceiling(dist_tl/100000)*100) %>% 
  group_by(group) %>% summarise(slope = lm(snowmelt ~ year)$coefficients[2], 
                                p = summary(lm(snowmelt ~ year))$coefficients[8],
                                n = n())

sm_trend$p_adj = p.adjust(sm_trend$p, method = 'BH')
p1 = ggplot(sm_trend) + theme_bw() + 
  theme(text = element_text(size = 7), plot.margin = unit(c(0,1,0,0), "mm"), aspect.ratio = 1) +
  geom_smooth(aes(group, slope*10, weight = n), method = 'lm', se = F, col = '#2171b5', fill = '#2171b5', size = 0.3) +  # plot trends/lines
  geom_point(aes(group, slope*10), col = '#9ecae1', size = 0.5) +  # plot the datapoints
  geom_point(data = filter(sm_trend, p < 0.05), aes(group, slope*10), col = '#2171b5', size = 0.5) +  # plot dark for significant datapoints
  labs(x = '', y = 'Trend (days/decade)') + coord_flip() +
  scale_x_reverse(breaks = rev(seq(0, 2000, by = 500))) +
  labs(x = 'Distance from tree line (km)')

ggsave(paste0(plotpath, 'snowmelt/snow_trend_NSIDC_lat_labB_', option, '.svg'), width = 5.5, height = 5, dpi = 300, units = 'cm')


## Fig S7 persistence/frequency trend -------------------------------------------------------
# as a threshold, we use the correlation between the pattern and extreme fire weeks
thresh_groups = inner_join(fire_3cl, V250_highfire_cor)
median(filter(thresh_groups, fire == 'more')$corV)
mean(filter(thresh_groups, fire == 'more')$corV)
quantile(filter(thresh_groups, fire == 'more')$corV)
quantile(filter(thresh_groups, fire == 'more')$corV, probs = c(0.25, 0.5, 0.6, 0.66, 0.75))
ggplot(thresh_groups) + theme_bw() +
  geom_rect(aes(xmin = 0.097, xmax = 0.25, ymin = -Inf, ymax = Inf), fill = 'gray88') +
  geom_rect(aes(xmin = 0.19, xmax = Inf, ymin = -Inf, ymax = Inf), fill = 'gray70') +
  stat_ecdf(data = filter(thresh_groups, fire == 'more'), aes(corV, y = 1 - ..y..), 
            geom = 'step', col = 'firebrick') +#, fill = 'firebrick', alpha = 0.5) + 
  stat_ecdf(data = filter(thresh_groups, !fire == 'more'), aes(x = corV, y = 1 - ..y..), 
            geom = 'step', col = 'black', fill = NA) +
  geom_text(aes(x = -0.3, y = 1.05), label = 'A', size = 10/.pt, fontface = "bold") +
  geom_vline(xintercept = 0.18, lty = 2, size = 0.3) +
  scale_y_continuous(limits = c(0,1.1), expand = c(0,0)) +
  labs(x = 'Correlation with Arctic front jet pattern', y = 'Cumulative fraction of data') +
  theme(text = element_text(size = theme.size))
ggsave(paste0(plotpath, 'v_trend/v250_thresh_groups_hist.svg'), width = 5.5, height = 5, dpi = 300, units = 'cm')

# compare trends for different estimates
quantile(filter(thresh_groups, fire == 'more')$corV, probs = c(0.25, 0.6))
res_pers = list()
for (thresh in seq(0.10, 0.19 , by = 0.01)){
  pers_weeks = filter(V250_highfire_cor, corV > thresh, week < 27) %>% group_by(year) %>% summarise(pers = n()) %>%
    full_join(data.frame(year = c(1979:2021))) %>% mutate(pers = ifelse(is.na(pers), 0, pers))
  perstest = summary(lm(pers ~ year, pers_weeks))
  res_pers = append(res_pers, list(c(thresh, perstest$coefficients[1], perstest$coefficients[2], perstest$coefficients[8])))
}
res_pers = as.data.frame(res_pers, row.names = c('thresh', 'intcpt', 'slope', 'p')) %>% t() %>% as.data.frame()
rownames(res_pers) <- c()

# compute frequency for 1980 and 2020
res_pers = res_pers %>% mutate(val_1980 = round((1980*slope)+intcpt, 2),
                               val_2020 = round((2020*slope)+intcpt, 2),
                               ratio = round(val_2020/val_1980, 2),
                               slope = round(slope, 3),
                               intcpt = round(intcpt, 2),
                               p = round(p, 3))

mean(res_pers$slope)
sd(res_pers$slope)
mean(res_pers$val_1980)
sd(res_pers$val_1980)
mean(res_pers$val_2020)
sd(res_pers$val_2020)
mean(res_pers$ratio)
sd(res_pers$ratio)

# Fig. 4C Trend in AFJ persistence ---------------------------
thresh = 0.18
paste0('prob for more class: ', 1-ecdf(filter(thresh_groups, fire == 'more')$corV)(thresh))
paste0('prob for avg class: ', 1-ecdf(filter(thresh_groups, fire == 'avg')$corV)(thresh))
pers_weeks = filter(V250_highfire_cor, corV > thresh, week < 27) %>% group_by(year) %>% summarise(pers = n()) %>%
  full_join(data.frame(year = c(1979:2021))) %>% mutate(pers = ifelse(is.na(pers), 0, pers))
summary(lm(pers ~ year, pers_weeks))

ggplot(pers_weeks, aes(year, pers)) + theme_bw() + 
  geom_point(col = '#ef3b2c', size = 0.3) +
  geom_smooth(method = 'lm', col = '#ef3b2c', fill = '#ef3b2c', size = 0.3) +
  xlim(c(1967, 2021)) +
  labs(x = '', y = 'Frequency of Arctic front jet pattern') +
  theme(text = element_text(size = theme.size), 
        plot.margin = unit(c(0,1,0,0), "mm"), aspect.ratio = 1)
ggsave(paste0(plotpath, 'v_trend/V250_trend_pers0', thresh*100, '_labC.svg'), width = 5.5, height = 5, dpi = 300, units = 'cm')



## Fig. S4: Influence of terrain type on fire - climate - snowmelt interactions -------------------------------------------

terrain = read_csv('2_pipeline/00_ecoregions/terrain_types_25.csv') %>% select(-X1) %>%
  filter(terrain > 0)

terr_orig = data.frame(terr_name = 
                         c('Eluvial', 'Colluvial', 'Diluvial colluvial', 'Boggy', 'Glacial Valley',
                           'Diluvial solufluction', 'Mid terrace', 'Inter-ridge lowland', 'Old Terrace',
                           'Marshes',  'Outwash', 'Low terrace','Inter-alas', 'Alas', 'Moraine', 
                           'High terrace', 'Valley-sea terrace'))
terr_trans = read_csv('2_pipeline/00_ecoregions/terrain_types_lookup.csv') %>% select(-X1) %>%
  cbind(terr_orig)

terrain = terr_trans %>% right_join(terrain, by = c('TM_ID_num'='terrain'))

fire_terrain_yak = inner_join(firedat, terrain) %>% filter(!is.na(terr_name))

# percentage of burning in study area included
burn_tot = firedat %>% mutate(burned = perc_burned*n_pix*0.25) %>% summarise(sum(burned))
burn_yak = fire_terrain_yak %>% mutate(burned = perc_burned*n_pix*0.25) %>% summarise(sum(burned))
100/burn_tot*burn_yak

# total burned area per terrain type
options(scipen=999)
terr_sum = fire_terrain_yak %>% mutate(burned = perc_burned*n_pix*0.25) %>%
  group_by(terr_name) %>% summarise(burned = sum(burned))
ggplot(terr_sum) + theme_bw() + labs(x = '', y = bquote('Burned area'~(km^2))) + 
  geom_col(aes(terr_name, burned)) + 
  geom_col(data = filter(terr_sum, terr_name %in% c('Alas','Boggy','Eluvial')),
           aes(terr_name, burned), fill = 'Firebrick') +
  geom_text(data = data.frame(x = 15.3, y = 140000, label = "A"), 
            aes(x = x, y = y, label = label), size = 8/.pt, fontface = 'bold') +
  coord_flip() +
  theme(text = element_text(size = theme.size))
ggsave(paste0(plotpath, 'supplement/terrain_types.svg'),  
       width = 7, height = 12, units = 'cm')

# percentage burned in three most common terrain types
100/burn_tot*sum(filter(terr_sum, terr_name %in% c('Eluvial', 'Boggy','Alas'))$burned)
100/burn_yak*sum(filter(terr_sum, terr_name %in% c('Eluvial', 'Boggy','Alas'))$burned)

# weekly compound plots per terrain type
firesum_terr = fire_terrain_yak %>% filter(terr_name %in% c('Alas','Boggy','Eluvial')) %>%
  group_by(terr_name, week, year) %>%
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to compute area (km2) from %
  summarise(burned = sum(perc_burned)) %>%
  # add missing weeks (when no burning happened)
  right_join(data.frame(week = rep(rep(14:26, 21),3),
                        year = rep(rep(2001:2021, each = 13),3),
                        terr_name = rep(c('Alas','Boggy','Eluvial'), each = 273)))%>%
  mutate(burned = ifelse(is.na(burned), 0, burned))
fire_anom_terr = group_by(firesum_terr, week, terr_name) %>% summarise(meanburn = mean(burned), sdburn = sd(burned)) %>%
  right_join(firesum_terr) %>% mutate(burnanom = (burned-meanburn)/sdburn,
                                 burnanom = ifelse(is.na(burnanom), 0, burnanom))
plotdata_terr = V250_highfire_cor %>%
  # join with burned area anomaly
  right_join(fire_anom_terr) %>%
  # join with snowmelt day anomaly
  left_join(mean_sm)

# plot
ggplot(plotdata_terr) + theme_bw() + 
  facet_wrap(~terr_name, ncol = 1,labeller = as_labeller(c('Alas'="", 'Boggy'="", 'Eluvial'=""))) +
  labs(x = 'Snowmelt anomaly (days)', y = 'Correlation with Arctic front jet', col = 'Burn \nanomaly') +
  geom_vline(xintercept = 0, col = 'gray', lty = 2) + 
  geom_hline(yintercept = 0, col = 'gray', lty = 2) +
  geom_point(aes(x = snowanom_days, y = corV, col = burnanom), size = 0.7) +
  scale_color_steps2(low = 'dodgerblue', high = 'darkred', mid = 'gray90', 
                     n.breaks = 7, midpoint = 0, 
                     guide=guide_colourbar(title.position="top", 
                                           title.theme = element_text( size = theme.size))) +
  geom_text(data = data.frame(x = 7, y = 0.4, terr_name = 'Alas', label = "B"), 
            aes(x = x, y = y, label = label), size = 8/.pt, fontface = 'bold') +
  geom_text(data = data.frame(x = 7, y = 0.4, terr_name = 'Boggy', label = "C"), 
            aes(x = x, y = y, label = label), size = 8/.pt, fontface = 'bold') +
  geom_text(data = data.frame(x = 7, y = 0.4, terr_name = 'Eluvial', label = "D"), 
            aes(x = x, y = y, label = label), size = 8/.pt, fontface = 'bold') +
  theme(text = element_text(size = theme.size), 
        plot.margin = unit(c(0,0,0.1,0.1), "mm"), aspect.ratio = 6/7,
        strip.background = element_blank(), strip.text.y = element_blank(),
        legend.position='right', legend.title.align=0.5, legend.margin=margin(0,0,0,0),
        legend.key.width = unit(0.05, "cm"), legend.key.height = unit(0.5, "cm"),
        legend.box.margin=margin(0,-5,-5,-5))
ggsave(paste0(plotpath, 'supplement/compound_terrain_types.svg'),  
       width = 7, height = 15, units = 'cm')

# annual models
fire_anom_year_terr = fire_anom_terr %>% group_by(year, terr_name) %>% summarise(burned = sum(burned))
fire_anom_year_terr = fire_anom_year_terr %>% 
  mutate(burnanom = (burned - mean(fire_anom_year_terr$burned))/sd(fire_anom_year_terr$burned))

# join plotdata yearly
plotdata_year_terr = V250_highfire_year %>%
  # join with burned area anomaly
  right_join(fire_anom_year_terr) %>%
  # join with snowmelt day anomaly
  left_join(mean_sm) 

# linear model
for (terr in unique(plotdata_year_terr$terr_name)){
  print(terr)
  dat = filter(plotdata_year_terr, terr_name == terr)
  fit = lm(burnanom ~ snowanom + meancorV, data = dat)
  print(summary(fit))
  print(cor.test(dat$burnanom, dat$snowanom, use = 'complete.obs', method = 'spearman'))
  print(cor.test(dat$burnanom, dat$meancorV, method = 'spearman'))
}



### when is there a relationship between snowmelt and jet? -----------------
snow65_mean = snow65 %>% summarise(mean = mean(week_med), first = min(week_med), last = max(week_med))
when = plotdata %>% group_by(week) %>% summarise(cor = cor.test(snow, corV, method = 'spearman')$estimate,
                                                 p = cor.test(snow, corV, method = 'spearman')$p.value)

ggplot(when, aes(week, cor)) + theme_bw() + 
  geom_hline(yintercept = 0, lty = 2, col = 'gray') +
  geom_vline(data = snow65_mean, aes(xintercept = mean), col = 'deepskyblue') +
  geom_vline(data = snow65_mean, aes(xintercept = first), col = 'deepskyblue', lty = 2) +
  geom_vline(data = snow65_mean, aes(xintercept = last), col = 'deepskyblue', lty = 2) +
  geom_point() + geom_line() + geom_point(data = filter(when, p < 0.1), aes(week, cor), col = 'red') +
  labs(x = 'Week since March 1', y = 'Cor (snow anomaly ~ AFJ pattern cor)') +
  theme(text = element_text(size = theme.size))
ggsave(paste0(plotpath, 'compound/ts_snow_V250_cor_all_MODIS.png'), 
       width = 5.5, height = 5.5, dpi = 300, units = 'cm')

# per month
plotdata %>% mutate(mon = ifelse(week %in% 14:17, 6, ifelse(week %in% 18:22, 7, 8))) %>%
  group_by(mon) %>% summarise(cor = cor.test(snowanom, corV, method = 'spearman')$estimate,
                              p = cor.test(snowanom, corV, method = 'spearman')$p.value)

# only beginning of June
plotdata %>% filter(week %in% 14:16) %>% ungroup() %>% 
  summarise(cor = cor.test(snowanom, corV, method = 'spearman')$estimate, 
            p = cor.test(snowanom, corV)$p.value, method = 'spearman')

