## extracting average FRP from MCD14

library(tidyverse)
library(sf)
library(lubridate)

setwd('D:/waves')

MODIS_scan_angle = function(samplei){
  # Retrieve MODIS 1 km active fire scan pixel sizes.
  # param samplei: scan sample number
  # return: dS (along-scan pixel size), dT (along-track pixel size)
  # pixdim_MODIS = np.array(MODIS_scan_angle(np.arange(1354))) * 463.312716528 * 2  # dS, dT in meters.
  if ((0 > samplei) | (samplei > 1353)){
    print('One or multiple sample number(s) not in range 0-1353')
  }
  Re = 6378.137  # earth radius in km.
  h = 705  # altitude of satellite in km.
  r = Re + h  # satellite track radius in km.
  s = 1.0 / h  # = p/h. p = pixel nadir resolution, h = altitude of satellite.
  N = 1354  # number of samples
  phi = -0.5 * N * s + 0.5 * s + (samplei - 1) * s  # zenith angle in radians. Ichoku & Kaufman (2005)
  # phi = abs(s * (i - 676.5))  # zenith angle in radians. MODIS fire documentation.
  dS = Re * s * ((cos(phi) / (sqrt((Re / r) ** 2 - (sin(phi)) ** 2))) - 1)  # dS, along-scan ground pixel size.
  dT = r * s * (cos(phi) - sqrt((Re / r) ** 2 - (sin(phi)) ** 2))  # dT, along-track ground pixel size.
  
  return(dS * dT)
}


larch = st_read('2_pipeline/00_ecoregions/eastern_siberia.shp')
dat_sub = data.frame(year = numeric(),
                     doy=numeric(),
                     frp_norm=numeric())

# first read ML active fires text files
filenames = list.files('D:/shared_data/MCD14MLv006_txt', '*txt', full.names = TRUE)
coltypes = rep('numeric', 12)
coltypes[12] = 'character'
for (file in filenames){
  data = read_delim(file, delim = ' ', col_names = TRUE,
                    col_types = cols_only(YYYYMMDD = 'i', lat = 'c', lon = 'c',
                              sample = 'c', FRP = 'c', conf = 'c', type = 'c')) %>% 
    mutate_at(c(2:7), as.numeric) %>%
    filter(type == 0, conf >= 30) %>%
    st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
    st_intersection(larch) %>% select(c(1:3))
  if (dim(data)[1]  == 0) next
  data =  mutate(data, area = MODIS_scan_angle(sample),
                 date = ymd(YYYYMMDD), year = year(date), doy = yday(date),
                 frp_norm = FRP/area) %>%
    select(c(7:9))
  dat_sub = rbind(dat_sub, data) # probably would be better to either process this in python or replace this loop with list.append

# add week number and filter for summer months
nyears = max(dat_sub$year)-min(dat_sub$year)+1
week_df = data.frame(year = rep(min(dat_sub$year):max(dat_sub$year), each = 30*7),
                     day = rep(1:(30*7), nyears)) %>%
  mutate(startday = yday(ymd(paste(year, 3, 2, sep = '-'))),
         doy = day + startday-1) %>% group_by(year) %>% arrange(year, doy) %>%
  bind_cols(week = rep(rep(1:30, each = 7), nyears))
datsub_week = dat_sub %>% mutate(startday = yday(ymd(paste(year, 3, 2, sep = '-'))),
                                 endday = startday + 30*7) %>%
  dplyr::filter(doy > startday, doy < endday) %>%
  left_join(select(week_df, year, doy, week)) %>% 
  filter(week > 13, week < 27)

# save to rds
saveRDS(datsub_week, paste('2_pipeline/01_modis/frp_larch_JJA_56_', collection, '.rds', sep = ''))
