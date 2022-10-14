# fire characteristics in Siberia
library(tidyverse)
library(viridis)
library(sf)
library(stars)

# workdir on laptop
wdir = "D:/waves/"
setwd(wdir)

### specs & load data ---------------------

# load fire/weather/wave data
load(paste0('2_pipeline/99_siberia/data_larch_bor_JJA_', res, '_56to74_fire.Rdata'))
states = st_read('0_data/naturalearth/states/ne_110m_admin_0_map_units.shp') %>%
  st_crop(xmin=-180, xmax=180, ymin=30, ymax=90)

# set plotting paths
plotpath = paste0('3_visuals/', res, 'deg/fire_climatology/')
options(scipen = 999)
theme.size = 7

## compute climatology ----------------------------------

# summarise fire data per latitude
latsum_fire = group_by(fire_week_1deg, lat, week, year) %>% 
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to turn area from % in km2
  summarise(burned = sum(perc_burned))

# compute climatology
latfire_clim = group_by(latsum_fire, lat, week) %>% summarise(meanclim = mean(burned), sdclim = sd(burned))

# BINNING to distance from treeline
tl = read_csv('2_pipeline/01_modis/distance_treeline_0.25deg.csv') %>% select(-X1)
latsum_fire_1deg = fire_week_1deg %>% left_join(tl) %>% mutate(group = ceiling(dist_tl/4)) %>% 
  group_by(group, week, year) %>% 
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to turn area from % in km2
  summarise(burned = sum(perc_burned))

# climatology aggregated
latfire_clim_1deg = group_by(latsum_fire_1deg, group, week) %>% summarise(meanclim = mean(burned), sdclim = sd(burned))

# compute individual years and weeks in comparison to climatology
latsum_fireclim = latsum_fire_1deg %>% left_join(latfire_clim_1deg) %>%
  mutate(burnclim = (burned-meanclim)/sdclim)

### Summing by tundra/taiga ------------------------
week_df = data.frame(year = rep(rep(2001:2021, 13),each = 2),
                     week = rep(rep(14:26, each = 7), each = 2),
                     reg = rep(c(255,2), 273))

fire_week_reg = firedat %>% group_by(week, year, reg) %>% 
  mutate(perc_burned = perc_burned*n_pix*0.25) %>% # this is to compute area (km2) from %
  summarise(burned = sum(perc_burned)) %>%
  # fill up missing weeks
  full_join(week_df, by = c('week','year','reg')) %>% mutate(burned = ifelse(is.na(burned), 0, burned)) %>%
  # compute cumulative sums over the year
  group_by(year, reg) %>% arrange(week) %>% mutate(cum_burn = cumsum(burned))

### Fig.1: Regional climatology taiga/tundra -----
mon_breaks = c(14, 18+3/7, 22+6/7)
mon_labs = c( 'Jun', 'Jul', 'Aug')
ggplot(fire_week_reg) + theme_bw() + 
  labs(x = '', y = bquote('Cumulative burned area'~(km^2)), fill = '') + 
  facet_wrap(~factor(reg, levels = c('255','2')), scales = 'free', 
             labeller = as_labeller(c(`255`="", `2`=""))) +
  theme(text = element_text(size = 7), 
        plot.margin = unit(c(0,0,0,0), "mm"), aspect.ratio = 9/10,
        strip.background = element_blank(), strip.text.y = element_blank(),
        legend.position = c(0.05, 0.9),legend.margin=margin(0,0,0,0), 
        legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.1, "cm"),
        legend.box.margin=margin(0,-5,-10,-5), legend.background = element_blank()) +
  scale_x_continuous(limits = c(13.5,26.5), breaks = mon_breaks-0.5,
                     labels = mon_labs, expand = c(0.01,0)) +
  geom_col(aes(x = week, y = cum_burn), 
           fill = 'lightgray', col = 'lightgray', position = 'stack', width = 0.5) +
  geom_col(data = filter(fire_week_reg, year %in% 2019:2021), 
           aes(week, cum_burn, fill = as.factor(year)), width = 0.5) +
  scale_fill_manual(values = c('#FEA873FF', '#E95562FF', '#67000d')) +
  geom_text(data = data.frame(x = 14, y = 1450, reg = 2, label = "B"), 
            aes(x = x, y = y, label = label), size = 10/.pt, fontface = 'bold') +
  geom_text(data = data.frame(x = 14, y = 620000, reg = 255, label = "A"), 
            aes(x = x, y = y, label = label), size = 10/.pt, fontface = 'bold')
ggsave(paste0(plotpath, 'fireclim_reg_cum.svg'),  
       width = 12, height = 5.5, dpi = 300, units = 'cm')

### Fig. S8: Fire climatology and line plots --------------------------
mon_breaks = c(1, 5+3/7, 9+5/7, 14, 18+3/7, 22+6/7, 27+2/7)
mon_labs = c('1 Mar', '1 Apr', '1 May', '1 Jun', '1 Jul', '1Aug', '1 Sep')

p1 = ggplot(filter(latfire_clim_1deg, meanclim > 0)) + theme_bw() + 
  labs(x = '', y = '', fill = bquote('Burned area'~(km^2))) +
  geom_raster(aes(week, group, fill = meanclim)) +
  scale_fill_viridis(option = 'A', direction = -1,
                     guide=guide_colourbar(title.position="top", 
                                           title.theme = element_text( size = 7))) + 
  coord_fixed() + 
  theme(text = element_text(size = 7), plot.margin = unit(c(0,0,0,0), "mm"), aspect.ratio = 9/10,
        legend.position="top", legend.title.align=0.5, legend.margin=margin(0,0,0,0), 
        legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.05, "cm"),
        legend.box.margin=margin(0,-5,-10,-5))
if (option == 'treeline'){
  p1 + geom_text(aes(x = 14, y = 1), label = 'A', size = 9/.pt, fontface = "bold") + 
    scale_x_continuous(breaks = mon_breaks[4:6]-0.5, limits = c(13.5, 26.5),
                       labels = mon_labs[4:6], expand = c(0.01,0)) +
    scale_y_reverse(breaks = rev(seq(0, 20, by = 5)), expand = c(0,0), limits = c(20, 0),
                    labels = c('20\u00B0', '15\u00B0', '10\u00B0', '5\u00B0', '0\u00B0'))
}else{
  p1 + geom_rect(xmin = 13.5, xmax = 26.5, ymin = 56, ymax = 75, col = 'gray40', fill = NA, size = 0.2) +
    scale_x_continuous(breaks = mon_breaks-0.5, #c(0.5, 4.5, 8.5, 13.5, 17.5, 21.5, 26.5), 
                       labels = mon_labs, expand = c(0.01,0)) +
    scale_y_continuous(limits = c(52, 75), breaks = seq(55, 75, by = 5), expand = c(0.01,0), 
                       labels = c('55\u00B0N', '60\u00B0N', '65\u00B0N', '70\u00B0N', '75\u00B0N')) +
    geom_text(aes(x = 2, y = 74), label = 'A', size = 9/.pt, fontface = "bold")
}
ggsave(paste0(plotpath, 'fireclimatology_70_A_', option, '.svg'),  
       width = 5.5, height = 5.5, dpi = 300, units = 'cm')

# local minimum in latitude burning over summer
fire_line_allyears = fire_week_1deg %>% filter(week > 13, week < 27) %>%
  mutate(perc_burned = perc_burned*n_pix*0.25,
         lat = floor(lat) + 0.5) %>%
  group_by(lat) %>% summarise(firesummer_tot = sum(perc_burned))
fire_line_1920 = fire_week_1deg %>% filter(year  > 2018, week > 13, week < 27) %>%
  mutate(perc_burned = perc_burned*n_pix*0.25,
         lat = floor(lat) + 0.5) %>%
  group_by(lat, year) %>% summarise(perc_burned = sum(perc_burned)) %>%
  left_join(fire_line_allyears) %>% mutate(firesummer = perc_burned/firesummer_tot*100)
fire_week_1deg %>% filter(week < 27) %>%
  mutate(seas = ifelse(week < 14, 'Spring', 'Summer'),
         lat = floor(lat) + 0.5) %>%
  group_by(lat, seas) %>% mutate(perc_burned = perc_burned*n_pix*0.25) %>%
  summarise(burned = sum(perc_burned)) %>% 
  ggplot() + theme_bw() + 
  labs(x = '', y = bquote('Seasonal burned area ('*~ km^2*')'), fill = '', lty = '') +
  geom_vline(xintercept = firelatlim + 0.5, col = 'gray', lty = 2) +
  geom_area(data = fire_line_1920, aes(lat, firesummer*1000, fill = as.factor(year)),
            col = NA, alpha = 0.5) + #, position = 'identity') +
  scale_fill_manual(values = c('#FEA873FF', '#E95562FF', '#67000d'), guide = F) +
  # geom_line(data = fire_line_1920, aes(lat, firesummer*600, lty = as.factor(year)), col = 'firebrick') + 
  scale_y_continuous(limits = c(0,100000), expand = c(0.01,0), 
                     sec.axis = sec_axis(~ ./1000, name = 'Cumulative burned area 2019-2021 (%)')) +
  scale_x_continuous(limits = c(52, 75), breaks = seq(55, 75, by = 5), expand = c(0.01,0), 
                     labels = c('55\u00B0N', '60\u00B0N', '65\u00B0N', '70\u00B0N', '75\u00B0N')) +
  geom_line(aes(lat, burned, lty = seas)) + 
  scale_linetype_manual(values = c(2,1), guide = F) + 
  geom_text(aes(y = 5000, x = 74), label = 'B', size = 9/.pt, fontface = "bold") +
  geom_text(aes(y = 15000, x = 71), label = '2020', size = 5/.pt, col = '#E95562FF') +
  geom_text(aes(y = 45000, x = 69), label = '2019', size = 5/.pt, col = '#FEA873FF') +
  geom_text(aes(y = 10000, x = 62.5), label = '2021', size = 5/.pt, col = '#67000d') +
  theme(text = element_text(size = theme.size), aspect.ratio = 9/10, plot.margin = unit(c(0,0,0,0), "mm"),
        legend.position = c(0.9, 0.9), legend.background = element_rect(fill = NA, colour = NA), legend.key.size = unit(4, 'point'),
        axis.line.x.top = element_line(color = "#E95562FF"), axis.ticks.x.top = element_line(color = "#E95562FF"),
        axis.text.x.top = element_text(color = "#E95562FF"), axis.title.x.top = element_text(color = "#E95562FF")) +
   coord_flip()

ggsave(paste0(plotpath, 'fireclim_summer_labels.svg'),  
       width = 5.8, height = 5.5, dpi = 300, units = 'cm')

### Fig. S8C: Study area circumpolar ------------------------------------------
studyarea = st_read('2_pipeline/00_ecoregions/eastern_siberia.shp')
# states2 = filter(states, ADM0_A3 %in% c('CAN', 'USA', 'RUS', 'NOR', 'GRL', 'SWE', 'FIN', 'DNK', 'ISL', 
#                                         'BLR', 'POL','LTU', 'LVA', 'EST', 'DEU', 'BEL', 'NLD', 'GBR'))
states2 = st_transform(states, 3995) %>%
  st_crop(xmin=-4700000, xmax=4700000, ymin=-4200000, ymax=4500000)
taiga_circum = st_read('D:/waves/0_data/ecoregions/borealbiome_wwf_te/boreal.shp') %>% st_union()
taiga_simple = st_simplify(taiga_circum, dTolerance = 0.25)
tundra_circum = st_read('D:/waves/0_data/ecoregions/cav/cp_veg_la.shp') %>% st_union()
tundra_simple = st_simplify(tundra_circum, dTolerance = 0.25)
permafrost = st_read('D:/waves/0_data/geographic/UiO_PEX_PERZONES_5.0_20181128_2000_2016_NH/UiO_PEX_PERZONES_5.0_20181128_2000_2016_NH.shp')
perma_simple = st_simplify(permafrost, dTolerance = 0.01) %>%
  mutate(name = ifelse(GRIDCODE == 1, 'Isolated', ifelse(GRIDCODE == 2, 'Sporadic', 
                ifelse(GRIDCODE == 3, 'Discontinuous', 'Continuous'))),
         name = factor(name, levels = c('Isolated', 'Sporadic', 'Discontinuous', 'Continuous')))
fire1920 = filter(fire_week_1deg, year %in% 2019:2020) %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326)

p1 = ggplot() + theme_bw() + labs(y = '', x = '') +
  geom_sf(data = perma_simple, aes(fill = name), col = NA) +
  scale_fill_brewer(palette = 'Blues', direction = 1) +
  geom_sf(data = states2, fill = NA, col = 'darkgray', size = 0.1) +
  geom_sf(data = taiga_simple, col = 'chartreuse4', fill = NA, size = 0.2) +
  geom_sf(data = tundra_simple, col = 'aquamarine3', fill = NA, alpha = 0.5, size = 0.3) +
  geom_sf(data = studyarea, col = 'firebrick', fill = NA, size = 0.3) +
  coord_sf(crs = 3995, xlim = c(-4200000, 4200000), ylim = c(-3800000, 4000000)) +
  scale_x_continuous(labels = function(x) paste0(x, '\u00B0', "W")) +
  scale_y_continuous(labels = function(x) paste0(x, '\u00B0', "N")) +
  geom_text(aes(y = Inf, x = -Inf), label = 'C', size = 10/.pt, fontface = "bold", hjust = -1, vjust = 1.5) +
  theme(text = element_text(size = 7), plot.margin = unit(c(0,0,0,0), "mm"), aspect.ratio = 9/10,
        panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.1),
        legend.margin=margin(0,0,0,0), legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.05, "cm"),
        legend.position = c(0.23, 0.13), legend.title = element_blank()) 
# ggsave(paste0(plotpath, 'studymap_AGU.png'), p1,  
#        width = 7, height = 7 , dpi = 300, units = 'cm')
ggsave(paste0(plotpath, 'studymap_S5.svg'), p1,  
       width = 6.1, height = 5.8, units = 'cm')

### Fig.1C: Zoomed study area ------------------------------
studyarea = st_read('2_pipeline/00_ecoregions/eastern_siberia.shp')

mat <- list(matrix(c(60, 56, 65, 56, 70, 56, 75, 56, 80, 56, 85, 56, 
                     90, 56, 95, 56, 100, 56, 105, 56, 110, 56, 115, 56, 
                     120, 56, 125, 56, 130, 56, 135, 56, 140, 56, 145, 56, 
                     150, 56, 155, 56, 160, 56, 165, 56, 170, 56, 175, 56, 
                     180, 56, 180, 90, 0, 90,0, 56, 60, 56), 
                   ncol = 2, byrow = TRUE))
box <- st_polygon(mat)
box <- st_geometry(box)
box <- st_set_crs(box, "+proj=longlat +datum=WGS84 +no_defs")
study_true = st_intersection(studyarea, box)


states2 = st_read('0_data/naturalearth/states/ne_110m_admin_0_map_units.shp') %>%
  st_crop(xmin=60, xmax=180, ymin=30, ymax=90) %>% st_transform(3576) %>%
  st_crop(xmin=-1000000, xmax=3800000, ymin=-4500000, ymax=100000)
ru_regions = st_read('D:/shared_data/basemaps/russia_admin_princeton/RUS_adm1_fix.shp') %>%
  st_crop(xmin=60, xmax=180, ymin=30, ymax=77.8) %>% st_transform(3576) %>%
  st_crop(xmin=-1000000, xmax=3800000, ymin=-4500000, ymax=100000)
  # filter(!NAME_1 %in% c('Krasnoyarsk','Chukot',"Arkhangel'sk", 'Yamal-Nenets'))
ru_clean = ru_regions
for (i in 1:dim(ru_regions)[1]){ # cleaned version without islands
  all_polys = st_cast(ru_regions$geometry[i], "POLYGON")
  largest = all_polys[st_area(all_polys) == max(st_area(all_polys))]
  ru_clean$geometry[i] = largest
}
tundra_circum = st_read('D:/waves/0_data/ecoregions/cav/cav_russia_line_diss.shp') %>% st_union()
tundra_simple = st_simplify(tundra_circum, dTolerance = 5000)

library(caTools)
library(reshape)

bartalev_path =  'D:/shared_data/bartalev/Binary/'
bartalev = read.ENVI(paste0(bartalev_path, 'ne_v4_2000.bil'), 
                     headerfile=paste0(bartalev_path, "ne_v4_2000_envi.hdr"))
lons = seq(from = -180.00446428, by = 0.00892857, length.out = dim(bartalev)[2])
lats = seq(from = 75.00446428, by = -0.0089285, length.out = dim(bartalev)[1])
dimnames(bartalev) <- list(lat = lats, lon = lons)
barta_tidy = melt(bartalev) %>% filter(lon > 80, lon < 180, lat < 72.2) %>%
  mutate(lat = round(lat*10)/10, lon = round(lon*10)/10) %>%
  group_by(lat,lon) %>% count(value) %>% slice(which.max(n))
barta_tidy = barta_tidy %>% mutate(value = ifelse(value == 7, 1, 0)) %>% 
  group_by(lat, lon) %>% summarise(value = max(value))
barta_sf = st_as_stars(barta_tidy, coords = c('lon','lat'), crs = 4326)
st_crs(barta_sf) = st_crs(4326)
barta_sf = barta_sf %>% st_transform(3576)


prast = ggplot() + theme_minimal() + labs(y = '', x = '', fill = '') +
  geom_stars(data = barta_sf, aes(fill = as.factor(value)), col=NA, alpha = 0.5) +
  scale_fill_manual(values=c('white','mediumseagreen'), guide = F) +
  coord_sf(crs = 3576, xlim = c(-750000,3500000), ylim = c(-4200000, -200000)) +
  theme(plot.margin = unit(c(0,0,0,0), "mm"), aspect.ratio = 9/10,
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) 
p1 = ggplot() + theme_bw() + labs(y = '', x = '', fill = '') +
  geom_sf(data = states2, fill = 'white', col = NA) +
  geom_sf(data = tundra_simple, col = 'gray20 ', lty = 1, size = 0.3) +
  geom_sf(data = ru_clean, fill = NA, col = 'gray80', size = 0.05) +
  geom_sf(data = states2, fill = NA, col = 'gray40', size = 0.1) +
  geom_sf(data = study_true, col = 'firebrick', fill = NA, size = 0.3) +
  coord_sf(crs = 3576, xlim = c(-750000,3500000), ylim = c(-4200000, -200000)) +
  scale_x_continuous(breaks = seq(50, 150, by = 10), labels = function(x) paste0(x, '\u00B0', "E")) +
  scale_y_continuous(breaks = seq(50, 90, by = 10), labels = function(x) paste0(x, '\u00B0', "N")) +
  geom_text(aes(y = Inf, x = -Inf), label = 'C', size = 10/.pt, fontface = "bold", hjust = -1, vjust = 1.5) +
  theme(text = element_text(size = theme.size), 
        plot.margin = unit(c(0,0,0,0), "mm"), aspect.ratio = 9/10,
        panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.1),
        legend.margin=margin(0,0,0,0), legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.05, "cm"),
        legend.position = c(0.15, 0.85), legend.title = element_blank()) 
ggsave(paste0(plotpath, 'studymap_detail_rast.png'), prast,  
       width = 5.5, height = 5.5 , dpi = 300, units = 'cm')
ggsave(paste0(plotpath, 'studymap_detail.svg'), p1,
       width = 5.5, height = 5.5 , units = 'cm')
