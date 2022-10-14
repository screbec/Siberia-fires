 # snow composite plots
library(tidyverse)
library(raster)
library(lubridate)
library(sf)

setwd('D:/waves/')

### Process nsidc snow data ----------------------
path_to_data = '0_data/snow/nsidc/data/data_bsq/'
files = list.files(path = path_to_data, pattern = '*.bsq')
for (filename in files[-1]){
  date = ymd(substr(filename, 21, 28))
  cond1 = TRUE
  cond2 = month(date) > 1
  cond3 = month(date) < 9
  if (cond1 & cond2 & cond3){
    # read file
    data = raster(paste0(path_to_data, filename))
    
    # read to tidy df and throw out nodata values
    data_df = as.data.frame(data, xy=TRUE) %>% 
      rename(snow = 3) %>% mutate(date = date) %>%
      filter(snow < 254)
    if( !exists('out') ){
      out = data_df
    }else{
      out = bind_rows(out, data_df)
    }
  }
}

# 0 means no snow on land, 1 means snow on land, all else means other
out = out %>% mutate(year = year(date), mon = month(date),
                     snow_on_off = ifelse(snow == 0, 0, ifelse(snow == 1, 1, NA)),
                     doy_off = ifelse(snow_on_off == 0, yday(date), NA)) 

snowmelt = out %>% 
  # filter all days with no snow on land
  filter(snow_on_off == 0) %>%
  # take the first yday for each pixel and year
  group_by(year, x, y) %>% slice(which.min(doy_off)) %>%
  # all areas without winter snow (snowmelt day before March 1) and having still snow on July 2nd (glaciers) are set to NA
  mutate(snowmelt = ifelse((mon > 2) & (mon < 8), doy_off, NA))

# compute anomalies
snowmelt_anom = group_by(snowmelt, x, y) %>% filter(!is.na(snowmelt)) %>% summarise(snow_clim = mean(snowmelt)) %>%
  right_join(snowmelt) %>% mutate(snow_anom = snowmelt-snow_clim) %>% filter(!is.na(snowmelt))

### clipping ------
# larch ecoregion shapefile for clipping
larch = st_read('2_pipeline/00_ecoregions/eastern_siberia.shp')

# transform points to geographic for filtering
data_sf_larch = st_as_sf(snowmelt_anom, coords = c('x', 'y')) %>% st_set_crs(st_crs(6931)) %>%
  # intersect with larch ecoregion
  st_intersection(st_transform(larch, 6931)) %>% st_transform(4326)

# convert back to simple dataframe with lat lon
data_larch = data_sf_larch %>% mutate(lon = sf::st_coordinates(.)[,1],
                                      lat = sf::st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL)

# yearly average in same region as fire (larch, >56 N)
sm_yak = data_larch %>% filter(lat > 56) %>%
  group_by(year) %>% summarise(snowmelt = mean(snowmelt))
sm_yak_anom = mutate(sm_yak, snow_anom = snowmelt-mean(sm_yak$snowmelt))

# 2001 - 2019 dataset only
sm_yak_2001 = filter(sm_yak, year > 2000)
sm_yak_anom_2001 = sm_yak_2001 %>% # this uses only 2001-2019 as reference for anomaly calculation
  mutate(sm_yak_2001, snow_anom = snowmelt-mean(sm_yak_2001$snowmelt))
data_larch_2001 = filter(data_larch, year > 2000)
  


# save to file
saveRDS(snowmelt_anom, file = '2_pipeline/05_snowmelt/sm_nsidc_full.RDS')
saveRDS(data_larch, file = '2_pipeline/05_snowmelt/sm_nsidc_larch_full.RDS')
saveRDS(sm_yak_anom, file = '2_pipeline/05_snowmelt/sm_nsidc_mean.RDS')

saveRDS(sm_yak_anom_2001, file = '2_pipeline/05_snowmelt/sm_nsidc_mean_2001.RDS')
saveRDS(data_larch_2001, file = '2_pipeline/05_snowmelt/sm_nsidc_larch_2001.RDS')
