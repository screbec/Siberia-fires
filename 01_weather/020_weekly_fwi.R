library(ncdf4)
library(reshape2)
library(tidyverse)
library(lubridate)
library(tictoc)
library(abind)
library(sf)

### read and process netCDFs into tidy dataframe

wdir = "D:/waves/" # working directory
var = 'fwi' # v250 u250 t2m tp cape rh1000


#### netcdf operations --------------------------------------------------------------------

readNCs <- function(files) {
  "loops through list of nc files, extracts date from filename and 
  returns a tidy dataframe "
  # iterate through the ncs
  for (i in 1:length(files)) {
    
    # open a connection to the ith nc file
    nc <- nc_open(files[i])
    
    # store to variables
    lon <- ncvar_get(nc, attributes(nc$dim)$names[2])
    lat <- ncvar_get(nc, attributes(nc$dim)$names[3])
    date = substr(files[i],15,22)
    
    # read variable and close connections
    temp <- ncvar_get(nc)
    
    # write lat lon and time into array dimnames
    dimnames(temp) <- list(lon = lon, lat = lat)
    # melt to tidy dataframe
    dat <- melt(temp, value.name = 'value') %>% na.omit() %>% 
      # filter(lon > 79, lon < 170, lat < 78) %>% 
      filter(lat > 50) %>% add_column(date = date)
    # close netcdf 
    nc_close(nc)
    
    # set the name of output variable and bind the new data to it
    if (exists("result")){
      result = bind_rows(result, dat)
    }else{
      result = dat
    }
  }
  return(result)
}

### create list of netcdf files
setwd(paste(wdir, '0_data/era5/', var, '/daily', sep = ''))
flist <- list.files(pattern = "*.nc$")


### read ncs into tidy dataframe of weekly averages over latitude band
tic()
data <- readNCs(flist)
toc()


# summarise to weekly
week_df = data.frame(mon = c(rep(6, 30), rep(7,31), rep(8,30)),
                     day = c(1:30,1:31,1:30),
                     week = rep(14:26, each = 7))
data = data %>% mutate(date = ymd(date), year = year(date), mon = month(date), day = day(date))%>% 
  left_join(week_df)%>% group_by(lon, lat, year, week) %>% summarise(value = mean(value))
data = data %>% group_by(lat, lon, week) %>% summarise(mean = mean(value), sd = sd(value)) %>%
  right_join(data, by = c('lat', 'lon', 'week')) %>% mutate(anom = (value-mean)/sd) %>%
  select(lat, lon, week, year, value, anom)
saveRDS(data, paste(wdir, '2_pipeline/02_era5_processing/fwi_data_raw.rds', sep = ''))

# clip by boreal area and summarise by area
boreal = st_read(paste0(wdir, '2_pipeline/00_ecoregions/eastern_siberia.shp'))
data_study = filter(data, lat > 56) %>%
  st_as_sf(coords = c('lon','lat'), crs = 4326) %>% 
  st_intersection(boreal) %>% select(year,week,value)
data_study = as.data.frame(st_coordinates(data_study)) %>% 
  bind_cols(as.data.frame(data_study)) %>% rename(lat = Y, lon = X)
  
# simplify dataframe and compute anomalies
data_study = data_study %>% group_by(lat, lon, week) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>%
  right_join(data_study, by = c('lat', 'lon', 'week')) %>% 
  mutate(anom = (value-mean)/sd) %>%
  select(lat, lon, week, year, value, anom) %>% na.omit()


# write to rds
saveRDS(data_study, paste(wdir, '2_pipeline/02_era5_processing/data', var, '.rds', sep = ''))

