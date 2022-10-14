# fire characteristics in Siberia
library(tidyverse)
library(sf)
library(scales) # for log trans axes

wdir = "D:/waves/"
setwd(wdir)

### specs & load data ---------------------

# load fire/weather/wave data
load(paste0('2_pipeline/99_siberia/data_larch_bor_JJA_', res, '_56to74_fire.Rdata'))

# load frp data
frp = readRDS(paste('2_pipeline/01_modis/frp_larch_JJA_56_006.rds', sep = ''))
st_geometry(frp) = NULL

# aggregate 0.25 degree fire data to lon/lat over full time
firetot_025deg_3cl = firedat %>% mutate(burned = perc_burned*n_pix*0.25) %>%
  right_join(select(fire_3cl, week, year, fire)) %>%
  group_by(lat, lon, fire) %>% summarise(burned = sum(perc_burned))

# load 0.25 degree Arctic Circle fire data
fire_AC = read_csv('2_pipeline/01_modis/perc_ba_0.25deg_week_70N_AC.csv') %>%
  select(-1) %>% mutate(perc_burned = perc_burned*n_pix*0.25) 

# load larch shapefile for overlay
boreal = st_read('2_pipeline/00_ecoregions/eastern_siberia.shp')

# set plotting paths
plotpath = paste0('3_visuals/', res, 'deg/') # general
options(scipen = 999)
theme.size = 7

## Table S1 ----------------------------
test = fire_3cl %>% mutate(burn_frac = round(burned/sum(fire_3cl$burned)*100,2)) %>% 
  filter(fire == 'more')
View(test)
sum(test$burn_frac)

## stats for Arctic circle and larch forests ------------------------

# summary stats for AC, yearly
firetot_AC = fire_AC %>% group_by(year) %>% summarise(burned = sum(perc_burned))
firetot_AC = firetot_AC %>% mutate(burned_perc = burned/sum(firetot_AC$burned))

# only siberia
firetot_AC_larch = fire_AC %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  st_intersection(boreal) %>% as.data.frame() %>%
  group_by(year) %>% summarise(burned_larch = sum(perc_burned)) %>%
  left_join(firetot_AC) %>% mutate(perc_larch = burned_larch/burned)
((firetot_AC_larch[19,4]*firetot_AC_larch[19,5])+
  (firetot_AC_larch[20,4]*firetot_AC_larch[20,5])+
  (firetot_AC_larch[21,4]*firetot_AC_larch[21,5]))/
  (firetot_AC_larch[19,4]+firetot_AC_larch[20,4]+firetot_AC_larch[21,4])

# what percentage is burned in each year
group_by(fire_3cl, year) %>% summarise(sum(burned), sum(burned)/sum(fire_3cl$burned), n()) %>% View()
# filter(fire_3cl, burned > 2000) %>% summarise(sum(burned)/sum(fire_3cl$burned), n())
# arrange(fire_3cl, desc(burned)) %>% mutate(csum = cumsum(burned)/sum(fire_3cl$burned)) %>% View()

# how much is burned in fire class 'more'
group_by(fire_3cl, fire) %>% summarise(sum(burned), sum(burned)/sum(fire_3cl$burned), n())

## burning above treeline per year (%) -------------------
tl_ba_tot = filter(firedat, reg == 2) %>% mutate(perc_burned = perc_burned*n_pix*0.25) %>% summarise(burned_tot = sum(perc_burned))
firedat %>% filter(reg == 2) %>% group_by(year) %>% 
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% summarise(burned = sum(perc_burned), perc = burned/tl_ba_tot)

## Fig. S1: FRP stats ----------------------------------------------------------------------

# frp per fire class
frp_2cl = left_join(frp, fire_3cl) %>% mutate(fire = ifelse(fire == 'more', 'more', 'other'))

# differences in mean between high fire and all others
ttest = t.test(frp_norm ~ fire, data = frp_2cl)
med = group_by(frp_2cl, fire) %>% summarise(median = median(frp_norm))

# plots
ggplot()+ theme_bw() + 
  stat_bin(data = filter(frp_2cl, fire == 'other'),
           aes(x = frp_norm, y = ..count../sum(..count..)),
           binwidth = 0.1, geom = 'line', col = 'darkgray') +
  stat_bin(data = filter(frp_2cl, fire == 'more'),
           aes(x = frp_norm, y = ..count../sum(..count..)),
           binwidth = 0.1, geom = 'line', col = 'Firebrick') +
  geom_vline(xintercept = ttest$estimate[1], col = 'Firebrick', lty = 2) +
  geom_vline(xintercept = ttest$estimate[2], col = 'darkgray', lty = 2) +
  # geom_vline(xintercept = med$median[1], col = 'Firebrick', lty = 3) +
  # geom_vline(xintercept = med$median[2], col = 'darkgray', lty = 3) +
  scale_x_continuous(trans = 'log10') +
  labs(x = bquote('Fire radiative power (MW/'*~ km^2*')'), y = 'Fraction of data') +
  theme(text = element_text(size = theme.size))
ggsave(paste0(plotpath, 'supplement/frp_3cl_power.png'), width = 5.5, height = 5.5, units = 'cm', dpi = 300)

## Fig.S9: compare our data to Xu et al 2022 -------------------------------------------

# split up wenxuans data between tundra and taiga
tundra = st_read('0_data/ecoregions/cav/test_wgs84.shp') %>% add_column(reg = 'tundra') %>% select(reg)
wenxuan_fire = read_csv('D:/waves/2_pipeline/01_modis/wenxuan_erl_ba.csv') %>% select(-geometry)
wenxuan_fire = wenxuan_fire %>% st_as_sf(coords = c('lon','lat'),crs=4326) %>% 
  st_join(tundra) %>% mutate(reg = ifelse(is.na(reg),'boreal',reg))
wenxuan_sum = wenxuan_fire %>% st_drop_geometry() %>% group_by(reg,year) %>% summarise(burned = n()*0.25)

# clip our data to Yakutia
yak = st_read('0_data/ecoregions/yakutia.shp') %>% select(NAME_1)
fire_yak = fire_week_1deg %>% st_as_sf(coords=c('lon','lat'),crs=4326) %>% # here we use data from March to end August!
  st_join(yak) %>% filter(NAME_1 == 'Sakha')
yak_sum = fire_yak %>% st_drop_geometry() %>% mutate(burned = perc_burned*n_pix*0.25,
                                                     reg = ifelse(reg == 255, 'boreal','tundra')) %>%
  group_by(reg,year) %>% summarise(burned = sum(burned))

# combine both datasets
combine = yak_sum %>% add_column(origin = 'MCD64') %>% filter(year > 2011, year < 2021) %>%
  bind_rows(wenxuan_sum %>% add_column(origin='MCD64+VNP14'))

ggplot(combine) + theme_bw() + 
  facet_wrap(~reg, scales = 'free_y', labeller = as_labeller(c('boreal'="", 'tundra'=""))) + 
  labs(x = '', y = bquote('Burned area'~(km^2)), fill = '') +
  geom_col(aes(year,burned, fill = origin), position=position_dodge()) +
  scale_x_continuous(expand=c(0.02,0.02)) + scale_y_continuous(expand=c(0.02,0.05)) +
  geom_text(data = data.frame(x = 2011.8, y = 70000, reg = 'boreal', label = "A"), 
            aes(x = x, y = y, label = label), size = 8/.pt, fontface = 'bold') +
  geom_text(data = data.frame(x = 2011.8, y = 920, reg = 'tundra', label = "B"), 
            aes(x = x, y = y, label = label), size = 8/.pt, fontface = 'bold') +
  theme(text = element_text(size = theme.size), 
        plot.margin = unit(c(0,0,0.1,0.1), "mm"), aspect.ratio = 9/10,
        strip.background = element_blank(), strip.text.y = element_blank(),
        legend.position = c(0.1, 0.9),legend.margin=margin(0,0,0,0), 
        legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.1, "cm"),
        legend.box.margin=margin(0,-5,-10,-5), legend.background = element_blank())
ggsave(paste0(plotpath, 'fire_climatology/comparison_erl_yak.svg'),  
       width = 12, height = 5.5, dpi = 300, units = 'cm')

# compute overall increase when adding VIIRS
combine_long = pivot_wider(combine, names_from = 'origin', values_from = burned) %>%
  mutate(perc_increase = 100/MCD64*`MCD64+VNP14`-100)
combine_long %>% group_by(reg) %>% mutate(MCD64 = ifelse(is.na(MCD64), 0, MCD64),
                                          `MCD64+VNP14` = ifelse(is.na(`MCD64+VNP14`), 0, `MCD64+VNP14`)) %>%
  summarise(MCD64 = sum(MCD64), `MCD64+VNP14` = sum(`MCD64+VNP14`)) %>%
  mutate(perc_increase = 100/MCD64*`MCD64+VNP14`-100)

# interannual correlation
combine_long %>% group_by(reg) %>% summarise(cor = cor.test(MCD64, `MCD64+VNP14`)$estimate,
                                             p = cor.test(MCD64, `MCD64+VNP14`)$p.value)


