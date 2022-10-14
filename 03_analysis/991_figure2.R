# local analysis for NE Siberia based on fire climatology
## === script for composite plots in Fig. 2 === ##
library(tidyverse)
library(viridis)
library(sf)
library(stars)
library(scales)

# workdir 
wdir = "D:/waves/"
setwd(wdir)

## load data ------------------

# load(paste0('2_pipeline/99_siberia/data_larch_bor_JJA_', res, '_56to74_newvars.Rdata'))
states = st_read('0_data/naturalearth/states/ne_110m_admin_0_map_units.shp') %>%
  st_crop(xmin=-180, xmax=180, ymin=30, ymax=90)
states2 = filter(states, ADM0_A3 %in% c('CAN', 'USA', 'RUS', 'NOR', 'GRL', 'SWE', 'FIN', 'DNK', 'ISL', 
                                        'BLR', 'POL','LTU', 'LVA', 'EST', 'DEU', 'BEL', 'NLD', 'GBR'))

# set plotting paths
plotpathhf = paste0('3_visuals/0.25deg/high_fire/') # only high fire group
theme.size = 7
comp_new = F # compute significance anew?


### compute significance -----------------------------------------------
comp_sig_1var = function(variable){
  ### variables:
  # variable: variable to process
  
  ## 1) read in the datasets
  wdir = "D:/waves/"
  
  if (variable == 'V'){
    data = readRDS(paste(wdir, '2_pipeline/02_era5_processing/data_totalwind.rds', sep = ''))
  }else if (variable == 'fwi'){
    data = readRDS(paste(wdir, '2_pipeline/02_era5_processing/fwi_data_raw.rds', sep = ''))
  }else{
    data = readRDS(paste(wdir, '2_pipeline/02_era5_processing/data', variable, '.rds', sep = ''))
  }
  
  fire_3cl = readRDS(paste(wdir, '2_pipeline/99_siberia/0.25deg/fire_3cl.rds', sep = '')) %>%
    dplyr::select(week, year, fire) %>% filter(fire == 'more')
  
  # 2) create groups: high-fire weeks and climatology
  # select only high-fire weeks
  data_hf = data %>% right_join(fire_3cl)
  # combine with all data (climatology)
  data_group = data_hf %>% bind_rows(data %>% add_column(fire = 'all'))
  # compute composite plot (mean of high fire weeks)
  hf_mean = data_hf %>% group_by(lon, lat) %>% summarise(anom = mean(anom))
  
  ## ---
  ## 3) compute t-test and field significance
  data_sig = data_group %>% unite('loc', lat, lon) %>% group_by(fire, loc) %>%
    summarise(anom = list(anom)) %>% spread(fire, anom) %>% group_by(loc)
  data_sig = data_sig %>% mutate(p_value = t.test(unlist(more), unlist(all))$p.value,
                                 p_val_sig = ifelse(p_value < 0.05, p_value, NA)) %>% ungroup()
  data_sig$p_adj = p.adjust(data_sig$p_val_sig, "BH")
  data_sig = dplyr::select(data_sig, loc, p_value, p_adj)
  
  
  ## 4) save significance result
  data_sig = data_sig %>% separate(loc, into = c('lat', 'lon'), sep = '_') %>%
    mutate_at(1:2, as.numeric)
  saveRDS(data_sig, paste0(wdir, '2_pipeline/99_siberia/0.25deg/sig_', variable, '.rds'))
  
  ## save composite plot
  saveRDS(hf_mean, paste0(wdir, '2_pipeline/99_siberia/0.25deg/mean_', variable, '_df.rds'))
  
  # convert composite plot to stars raster and save
  mean_out_rast = st_as_stars(hf_mean, coords = c('lon','lat'), crs = 4326)
  st_crs(mean_out_rast) = st_crs(4326)
  mean_out_rast = mean_out_rast %>% st_transform_proj(3995)
  saveRDS(mean_out_rast, paste0(wdir, '2_pipeline/99_siberia/0.25deg/mean_', variable, '.rds'))
}

if (comp_new){
  for (variable in c('t2m','tp','cape','z500','V','vpd')){
    comp_sig_1var(variable)
    gc()
  }
}

## 2A plot Meteorology only high fire with stippling -----------------------

# background map of states
ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm")) +
  geom_sf(data = states2, fill = NA, col = 'gray20', size = 0.5/.pt) +
  coord_sf(crs = 3995, xlim = c(-4200000, 4200000), ylim = c(-4200000, 4200000))
ggsave(paste0(plotpathhf, 'countries.svg'),
       width = 5.8, height = 5.8, units = 'cm')
  
for (var0 in c('t2m','tp','fwi','z500','V','vpd')){ 
  met_V_anom_rast = readRDS(paste0('2_pipeline/99_siberia/0.25deg/mean_', var0, '.rds'))
  sig_out = readRDS(paste0('2_pipeline/99_siberia/0.25deg/sig_', var0, '.rds')) # for field significance
  
  print(var0)
  # extract maximum absolute value
  vmax = ifelse(abs(min(met_V_anom_rast$anom, na.rm = TRUE)) > max(met_V_anom_rast$anom, na.rm = TRUE),
                abs(min(met_V_anom_rast$anom, na.rm = TRUE)), max(met_V_anom_rast$anom, na.rm = TRUE))
  pal = ifelse(var0 == 't2m', 'RdBu', ifelse(var0 == 'vpd', 'BrBG', 
        ifelse(var0 == 'z500', 'PiYG', 'PuOr')))
  lab = ifelse(var0 == 't2m', 'A', ifelse(var0 == 'tp', 'B', ifelse(var0 == 'vpd', 'C', 
        ifelse(var0 == 'fwi', 'D', ifelse(var0 == 'z500', 'E', 'F')))))
  test = met_V_anom_rast %>% dplyr::select(anom) %>% adrop()
  
  # create base plot & legend (separately)
  prast = ggplot() + theme_void() +
    geom_stars(data = test[1]) +
    coord_sf(crs = 3995, xlim = c(-4200000, 4200000), ylim = c(-4200000, 4200000)) +
    labs(col = '', fill = '') +
    theme(text = element_text(size = theme.size), plot.margin = unit(c(0,0,0,0), "mm"))
  plegend = ggplot() + theme_void() +
    geom_stars(data = test[1]) + labs(col = '', fill = '') +
    theme(text = element_text(size = theme.size), plot.margin = unit(c(0,0,0,0), "mm"),
          legend.position="bottom", legend.title.align=0.5,
          legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.1, "cm"))
  # add color palette
  if (var0 %in% c('t2m','vpd','V','z500')){
    prast = prast + scale_fill_distiller(palette = pal, na.value="white", 
                                   limits = c(vmax*(-1), vmax), guide = F)
    plegend = plegend + scale_fill_distiller(palette = pal, na.value="white", 
                                           limits = c(vmax*(-1), vmax))
  }else if (var0 == 'fwi'){
    prast = prast + scale_fill_gradient2(high = '#ae017e', low = '#41b6c4', mid = 'white', 
                                   na.value="white", midpoint = 0, 
                                   limits = c(vmax*(-1), vmax), guide = F)
    plegend = plegend + scale_fill_gradient2(high = '#ae017e', low = '#41b6c4', mid = 'white', 
                                         na.value="white", midpoint = 0, 
                                         limits = c(vmax*(-1), vmax))
  }else{
    prast = prast + scale_fill_gradient2(high = '#081d58', low = '#662506', mid = 'white', 
                                   na.value="white", midpoint = 0, guide = F)
    plegend = plegend + scale_fill_gradient2(high = '#081d58', low = '#662506', mid = 'white', 
                                         na.value="white", midpoint = 0)
  }
  
  # stippling vector img
  sig_contour = sig_out %>% filter(!is.na(p_adj)) %>% 
    st_as_sf(coords = c('lon','lat'), crs = 4326) %>% st_transform(3995) 
  sig_downscale = cbind(as.data.frame(sig_contour), as.data.frame(sf::st_coordinates(sig_contour))) %>%
    mutate(X = round(X/100000, 0)*100000, Y = round(Y/100000, 0)*100000) %>% dplyr::select(X, Y) %>% distinct() %>% 
    st_as_sf(coords = c('X','Y'), crs = 3995)
  p1 = ggplot() + theme_void() + theme(plot.margin = unit(c(0,0,0,0), "mm")) +
    geom_sf(data = sig_downscale, col = 'gray40', size = 0.3/.pt) +
    coord_sf(crs = 3995, xlim = c(-4200000, 4200000), ylim = c(-4200000, 4200000)) +
    labs(col = '', fill = '')
    # geom_text(aes(x = -4200000, y = 4200000), label = lab, size = 10/.pt, fontface = "bold") + 
  
  # colorbar only 
  leg <- get_legend(plegend)
  as_ggplot(leg)
  ggsave(paste0(plotpathhf, var0, '_hf_longclim_legend.svg'),
         width = 5.5, height = 2, units = 'cm')
  
  ggsave(paste0(plotpathhf, var0, '_hf_longclim.png'),
         prast, width = 5.8, height = 5.8, dpi = 300, units = 'cm')
  ggsave(paste0(plotpathhf, var0, '_hf_longclim_sig2.svg'),
         p1, width = 5.8, height = 5.8, units = 'cm')
  
  
}

### Fig. S5A: v250 pattern in 2021/2020/2019 ----------------------------------
fire_3cl = readRDS(paste(wdir, '2_pipeline/99_siberia/0.25deg/fire_3cl.rds', sep = ''))
V250_anom = readRDS(paste(wdir, '2_pipeline/02_era5_processing/data_totalwind.rds', sep = ''))

yr = 2021

hf_weeks_yr = filter(fire_3cl, year == yr, fire == 'more')$week
V250_yr = filter(V250_anom, year == yr, week %in% hf_weeks_yr) %>%
  group_by(lon, lat) %>% summarise(V250anom = mean(anom)) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)

vmax = ifelse(abs(min(V250_yr$V250anom, na.rm = TRUE)) > max(V250_yr$V250anom, na.rm = TRUE),
              abs(min(V250_yr$V250anom, na.rm = TRUE)), max(V250_yr$V250anom, na.rm = TRUE))
ggplot() + theme_void() +
  geom_sf(data = V250_yr, aes(col = V250anom)) + 
  geom_sf(data = states2, fill = NA, col = 'gray40', size = 0.3) +
  coord_sf(crs = 3995, xlim = c(-4200000, 4200000), ylim = c(-4200000, 4200000)) +
  labs(col = '', fill = '')+ 
  scale_color_distiller(palette = 'PuOr', na.value="white", limits = c(vmax*(-1), vmax)) +
  theme(text = element_text(size = theme.size), plot.margin = unit(c(0,1,1,1), "mm"),
        legend.position="bottom", legend.title.align=0.5,
        legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.1, "cm"))
ggsave(paste0(plotpathhf, 'V250anom_hf_', yr, '.png'),
       width = 5.5, height = 5.5, dpi = 300, units = 'cm')

### Fig. S5B: geopotential height anomalies in 2021/2020/2019 ----------------
z500 = readRDS(paste(wdir, '2_pipeline/02_era5_processing/dataz500.rds', sep = ''))

yr = 2021
hf_weeks_yr = filter(fire_3cl, year == yr, fire == 'more')$week
z500_yr = filter(z500, year == yr, week %in% hf_weeks_yr) %>%
  group_by(lon, lat) %>% summarise(z500 = mean(anom)) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326)
vmax = ifelse(abs(min(z500_yr$z500, na.rm = TRUE)) > max(z500_yr$z500, na.rm = TRUE),
              abs(min(z500_yr$z500, na.rm = TRUE)), max(z500_yr$z500, na.rm = TRUE))
ggplot() + theme_void() +
  geom_sf(data = z500_yr, aes(col = z500)) + 
  geom_sf(data = states2, fill = NA, col = 'gray40', size = 0.3) +
  coord_sf(crs = 3995, xlim = c(-4200000, 4200000), ylim = c(-4200000, 4200000)) +
  labs(col = '', fill = '')+ 
  scale_color_distiller(palette = 'PiYG', na.value="white", limits = c(vmax*(-1), vmax)) +
  theme(text = element_text(size = theme.size), plot.margin = unit(c(0,1,1,1), "mm"),
        legend.position="bottom", legend.title.align=0.5,
        legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.1, "cm"))
ggsave(paste0(plotpathhf, 'z500anom_hf_', yr, '.png'),
       width = 5.5, height = 5.5, dpi = 300, units = 'cm')

