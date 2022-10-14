# local analysis for NE Siberia (larch ecotone) based on fire climatology
## === script for loading and computing climatologies === ##

# load packages
library(terra)
library(tidyverse)
library(sf)
library(lubridate)

# workdir
wdir = "D:/waves/"
setwd(wdir)

### specs
firelatlim = 56   ## firelatlim is the lower boundary for the fire group calculation (55 or 59)
firelatlim2 = 74  ## firelatlim2 is the upper boundary (fixed to 74 for inclusion of active fires)
endyear = 2021
# set start and end week to JJA
startweek = 13
endweek = 27


## Process LIGHTNING data -----------------------------------
lightdat = read_csv('2_pipeline/04_lightning/lightning_eastern_siberia.csv')

# filter lightning and fire data by minimum latitude and time frame
lightdat = dplyr::filter(lightdat, lat > firelatlim)

# compute weekly sums of lightning
week_df = data.frame(year = rep(2012:2021, each = 30*7),
                     day = rep(1:(30*7), 10)) %>%
  mutate(startday = yday(ymd(paste(year, 3, 2, sep = '-'))),
         doy = day + startday-1) %>% group_by(year) %>% arrange(year, doy) %>%
  bind_cols(week = rep(rep(1:30, each = 7), 10))
lightdat_week_latlon = lightdat %>% mutate(year = year(date), doy = yday(date),
                                    startday = yday(ymd(paste(year, 3, 2, sep = '-'))),
                                    endday = startday + 30*7) %>%
  dplyr::filter(doy > startday, doy < endday) %>%
  left_join(select(week_df, year, doy, week)) %>%
  filter(week > startweek, week < endweek)
lightdat_week = lightdat_week_latlon %>% group_by(year, week) %>% summarise(lightn = n()) %>%
  # add weeks without lightning strikes for climatology
  right_join(data.frame(week = rep((startweek+1):(endweek-1), 21),
                        year = rep(2001:2021, each = endweek-startweek-1))) %>%
  mutate(lightn = ifelse(is.na(lightn), 0, lightn))

# check which strikes are in tundra region
polar = st_read('0_data/ecoregions/cav/test_wgs84.shp') %>% add_column(reg = 'tundra') %>%
  select(reg)
lightdat_reg = lightdat_week_latlon %>% st_as_sf(coords = c('lon','lat'), crs = 4326) %>%
  st_join(polar) %>% mutate(reg = ifelse(is.na(reg), 'boreal', reg)) %>% st_drop_geometry() %>%
  group_by(week, year, reg) %>% summarise(lightn = n())

# lightning climatology
light_anom = lightdat_week %>% group_by(week) %>% summarise(meanlight = mean(lightn), sdlight = sd(lightn)) %>%
  right_join(lightdat_week) %>% mutate(lightanom = (lightn-meanlight)/sdlight)

save(list = c('light_anom', 'lightdat_reg'),
     file = paste0('2_pipeline/99_siberia/data_larch_bor_JJA_0.25_', firelatlim, 'to', firelatlim2, '_lightning.Rdata'))




## fire data and extreme weeks --------------------------------------------------
# load FIRE data
fire_week_1deg = read_csv(paste0('2_pipeline/01_modis/perc_ba_0.25deg_week_70N_larch.csv')) %>% select(-1)

# filter lightning and fire data by minimum latitude and time frame
firedat = fire_week_1deg %>% filter(lat > firelatlim, week > startweek, week < endweek)

# summarise fire data per week and year
firesum = group_by(firedat, week, year) %>%
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to compute area (km2) from %
  summarise(burned = sum(perc_burned)) %>%
  # add missing weeks (when no burning happened)
  right_join(data.frame(week = rep((startweek+1):(endweek-1), endyear-2000),
                        year = rep(2001:endyear, each = endweek-startweek-1))) %>%
  mutate(burned = ifelse(is.na(burned), 0, burned))

# compute fire climatology and anomalies
fire_anom = group_by(firesum, week) %>% summarise(meanburn = mean(burned), sdburn = sd(burned)) %>%
 right_join(firesum) %>% mutate(burnanom = (burned-meanburn)/sdburn,
                                burnanom = ifelse(is.na(burnanom), 0, burnanom))

# create groups of fire activity based on sd
fire_3cl = fire_anom %>%
  mutate(fire = ifelse(burnanom >1, 'more', ifelse(burnanom >= (-1), 'avg', 'less')),
         fire = fct_relevel(factor(fire), 'less'))
saveRDS(fire_3cl, paste(wdir, '2_pipeline/99_siberia/0.25deg/fire_3cl.rds', sep = ''))
save(list = c('fire_anom', 'fire_3cl', 'firedat', 'firesum', 'fire_week_1deg', 'firelatlim'),
     file = paste0('2_pipeline/99_siberia/data_larch_bor_JJA_0.25_', firelatlim, 'to', firelatlim2, '_fire.Rdata'))




### Compute AFJ pattern correlation ------------------------------------------------
V250_anom = readRDS(paste(wdir, '2_pipeline/02_era5_processing/data_totalwind.rds', sep = ''))
meanV = readRDS('2_pipeline/99_siberia/0.25deg/mean_V_df.rds') %>% rename(Vhighfire = anom)

# correlation with meanV
V250_highfire_cor = meanV %>%
  # join with all weeks
  left_join(V250_anom) %>%
  # compute correlation of each week with the mean high fire V
  group_by(year, week) %>% summarise(corV = cor(anom, Vhighfire, method = 'spearman'))

saveRDS(V250_highfire_cor, paste(wdir, '2_pipeline/99_siberia/0.25deg/V250_highfire_cor.rds', sep = ''))
