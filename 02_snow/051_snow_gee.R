library(raster)
library(tidyverse)

setwd('D:/waves')


### Create sm file ------------------------
path_to_file = '0_data/snow/gee_ldos/'
# open your raster with n bands
r_deg <- stack(paste0(path_to_file, 'ldos_larch_15_2021_mosaic_clip_025deg.tif'))
r <- stack(paste0(path_to_file, 'ldos_larch_15_2021_mosaic_clip.tif'))

for(i in 1:nlayers(r)){
  print(i+2000)
  # access band
  band <- r[[i]]
  band_agg <- r_deg[[i]]

  # compute regional average in >56N
  data_avg = as.data.frame(band, xy=TRUE) %>% rename(snow = 3) %>%
    # only include snowmelt days since beginning fo march and end of july
    filter(snow > 59, snow < 213, x > 56) %>% 
    summarise(snow = mean(snow)) %>% mutate(year = 2000 + i)

  # aggregate for map
  # band_agg2 = aggregate(band, 2, fun=mean)
  # band_agg10 = aggregate(band, 10, fun=mean)

  # read aggregated data to tidy df
  data_map = as.data.frame(band_agg, xy=TRUE) %>% rename(snow = 3) %>%
    # only include snowmelt days between March 1 and July 1
    filter(snow > 59, snow < 183) %>%
    mutate(year = 2000 + i)

  # add to outfiles
  if(i == 1){
    mean_sm = data_avg
    sm_map = data_map
  }else{
    mean_sm = bind_rows(mean_sm, data_avg)
    sm_map = bind_rows(sm_map, data_map)
  }

  # save raster in a separate file
  # writeRaster(band, paste0(path_to_file, 'ldos_sib_mosaic_larch_band',i,'.tif', sep=''))
}
rm(band_agg, band, data_avg, data_map, r, i)

# compute anomalies
sm_map = sm_map %>% group_by(x, y) %>% summarise(snowmean = mean(snow)) %>%
  full_join(sm_map) %>% mutate(snowanom = snow-snowmean)
mean_sm = mean_sm %>% mutate(snowanom = (snow - mean(mean_sm$snow))/sd(mean_sm$snow),
                             snowanom_days = snow - mean(mean_sm$snow))

# save to file
saveRDS(sm_map, file = '2_pipeline/05_snowmelt/sm_gee_025deg.RDS')
saveRDS(mean_sm, file = '2_pipeline/05_snowmelt/sm_gee_mean.RDS')

