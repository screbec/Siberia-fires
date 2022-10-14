library(tidyverse)
library(sf)
library(ggpubr)
wdir = "D:/waves/"
setwd(wdir)
theme.size = 7 # text size for plottingF
geom.text.size = theme.size / (14/5)

# set plotting path
plotpath = paste0('3_visuals/0.25deg/')
process = TRUE

## load data and compute weekly means over study region
boreal = st_read('2_pipeline/00_ecoregions/eastern_siberia.shp')
fire_3cl = readRDS(paste(wdir, '2_pipeline/99_siberia/0.25deg/fire_3cl.rds', sep = ''))
V250_highfire_cor = readRDS('2_pipeline/99_siberia/0.25deg/V250_highfire_cor.rds')

if (process == TRUE){
  z500 = readRDS('2_pipeline/02_era5_processing/dataz500.rds') %>% rename(z500 = anom) %>% select(-value) %>%
    filter(lon > 79, lon < 170, lat < 78, lat > 56) # rough prefilter
  t2m = readRDS('2_pipeline/02_era5_processing/datat2m.rds') %>% rename(t2m = anom) %>% select(-value) %>%
    filter(lon > 79, lon < 170, lat < 78, lat > 56) # rough prefilter
  vpd = readRDS('2_pipeline/02_era5_processing/datavpd.rds') %>% rename(vpd = anom) %>% select(-value) %>%
    filter(lon > 79, lon < 170, lat < 78, lat > 56) # rough prefilter
  tp = readRDS('2_pipeline/02_era5_processing/datatp.rds') %>% rename(tp = anom) %>% select(-value) %>%
    filter(lon > 79, lon < 170, lat < 78, lat > 56) # rough prefilter
  
  # weekly mean meteorology in study region
  all_met_week = full_join(z500, t2m) %>% full_join(vpd) %>% full_join(tp) %>% 
    st_as_sf(coords = c('lon','lat'), crs = 4326) %>% st_intersection(boreal) %>%
    as.data.frame() %>% group_by(week, year) %>% 
    summarise(temp = mean(t2m), vpd = mean(vpd), z500 = mean(z500), tp = mean(tp)) 
  
  # add FWI
  fwi = readRDS('2_pipeline/02_era5_processing/datafwi.rds') %>%
    group_by(week, year) %>% summarise(fwi = mean(anom))
  all_met_week = all_met_week %>% full_join(fwi)
  
  # save as rds
  saveRDS(all_met_week, '2_pipeline/99_siberia/0.25deg/all_met_week.rds')
  
}else{
  all_met_week = readRDS('2_pipeline/99_siberia/0.25deg/all_met_week.rds')
}



### Fig S3 is arctic front jet pattern correlated to temperature and vpd
afj_meteo = full_join(fire_3cl, V250_highfire_cor) %>% 
  full_join(pivot_longer(all_met_week, cols = 3:7, names_to = 'var', values_to = 'anom')) 
options(scipen = 0)
cnt = 0
for (var0 in unique(afj_meteo$var)){
  cnt = cnt + 1
  p1 = ggplot(filter(afj_meteo, var == var0), aes(corV, anom)) + theme_bw() +
    geom_smooth(method = 'lm', col = 'gray60', size = 0.6) +
    geom_point(aes(col = burnanom), size = 0.6) + 
    scale_color_steps(breaks = c(0, 1, 2, 3), low = '#fee5d9', high = '#a50f15', guide = F) +
    labs(x = 'Correlation with Arctic front jet', y = 'Anomaly') +
    # geom_text(aes(x = -0.42, y = 1.4), label = toupper(letters[cnt]), size = 9/.pt, fontface = "bold", hjust = 0.2) + 
    theme(text = element_text(size = theme.size), 
          plot.margin = unit(c(0,1,1,1), "mm"), aspect.ratio = 6/7)
  ggsave(paste0(plotpath, 'supplement/afj_meteo_scatter_', toupper(letters[cnt]), '.svg'), p1,
         width = 5.8, height = 5.8, dpi = 300, units = 'cm')
}
for (var0 in unique(afj_meteo$var)){
  cnt = cnt + 1
  p1 = ggplot(filter(afj_meteo, var == var0), aes(anom, burnanom)) + theme_bw() + 
    geom_point(col = 'gray', size = 0.6) + 
    geom_smooth(method = 'lm', col = 'gray60', size = 0.6) +
    stat_cor(label.y = 2.9, size = geom.text.size, p.accuracy = 0.001, method = 'spearman') + 
    labs(x = 'Anomaly', y = 'Burned area anomaly') +
    geom_text(aes(x = -Inf, y = 3.5), label = toupper(letters[cnt]), 
              size = 9/.pt, fontface = "bold", hjust = -1) + 
    theme(text = element_text(size = theme.size), 
          plot.margin = unit(c(0,1,1,1), "mm"), aspect.ratio = 6/7)
  # ggsave(paste0(plotpath, 'supplement/afj_meteo_scatter_', toupper(letters[cnt]), '.png'), p1,
  #        width = 5.5, height = 5.5, dpi = 300, units = 'cm')
}
# only colorbar
p = ggplot(filter(afj_meteo, var == var0), aes(corV, anom, col = burnanom)) + 
  geom_point() + theme_minimal() + labs(color = 'Burned area\nanomaly') +
  scale_color_steps(breaks = c(0, 1, 2, 3), low = '#fee5d9', high = '#a50f15')
leg <- get_legend(p)
as_ggplot(leg)
ggsave(paste0(plotpath, 'supplement/afj_meteo_scatter_cb.svg'),
       width = 5.8, height = 5.8, dpi = 300, units = 'cm')

options(scipen=999)
group_by(afj_meteo, var) %>% 
  summarise(p = cor.test(anom, burnanom,method='spearman')$p.value,
            rho = cor.test(anom, burnanom,method='spearman')$estimate)
