library(ncdf4)
library(reshape2)
library(tidyverse)
library(lubridate)
library(tictoc)
library(abind)


### read and process netCDFs into tidy dataframe

wdir = "D:/waves/" # working directory
var = 'tp' # v250 u250 t2m tp cape rh1000


#### netcdf operations --------------------------------------------------------------------

readNCs <- function(files, var) {
  "loops through list of nc files, calculates the weekly mean of each variable 
  for each lon/lat gridcell and returns a tidy dataframe with all variables as columns"
  # iterate through the ncs
  for (i in 1:length(files)) {
    print(files[i])
    
    # open a connection to the ith nc file
    nc <- nc_open(files[i])
    
    # store to variables
    lon <- ncvar_get(nc, attributes(nc$dim)$names[3])
    lat <- ncvar_get(nc, attributes(nc$dim)$names[4])
    time <- ncvar_get(nc, attributes(nc$dim)$names[1])
    yr0 = year(hours(time[1]) + ymd_hms("1900-01-01 00:00:0.0"))
    yr = year(hours(tail(time, n = 1)) + ymd_hms("1900-01-01 00:00:0.0"))
    varname <- attributes(nc$var)$names[2]
    
    # read variable and close connections
    temp <- ncvar_get(nc, varname)
    
    ## --- for single level variables (v250/u250)
    if (yr0 == yr){
      ## --- only when one nc file per year
        # write lat lon and time into array dimnames
        dimnames(temp) <- list(lon = lon, lat = lat, time = time)
        # clip first day of March and last three days of September to have 30 weeks
        temp = temp[ , ,1:91]
        # melt to tidy dataframe
        dat <- melt(temp, value.name = varname) %>%
          # add a column for week
          arrange(time) %>% add_column(week = rep(14:26, each = 7*length(lon)*length(lat)*length(levels))) %>%
          # take mean per week
          group_by(lon, lat, week) %>% summarise(time = first(time), value = mean(!!sym(varname))) %>%
          mutate(year = yr) %>% select(-time)
      
      ## multiple years in one file
      }else{
        # write lat lon and time into array dimnames
        dimnames(temp) <- list(lon = lon, lat = lat, time = time)
        year_real = year(hours(time) + ymd_hms("1900-01-01 00:00:0.0"))
        for (nyear in unique(year_real)){
          # split by year
          tempyr = temp[ , ,which(year_real == nyear)]
          # clip last day of august to have 13 weeks
          if (dim(tempyr)[3] == 92){
            tempyr = tempyr[ , ,1:91]
          }else{
            print('wrooooong!!')
          }
          
          # reshape to long format
          tempyr = melt(tempyr, value.name = varname) %>%
            # add a column for week
            arrange(time) %>% add_column(week = rep(14:26, each = 7*length(lon)*length(lat)*length(levels))) %>%
            # take mean per week
            group_by(lon, lat, week) %>% summarise(time = first(time), value = mean(!!sym(varname))) %>%
            mutate(year = nyear) %>% select(-time)
            
          if (exists('datyrl')) {datyrl = bind_rows(datyrl, tempyr)} else {datyrl = tempyr}
          }
        dat = datyrl
        rm(datyrl)
      }
    
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
data <- readNCs(flist, var = var) #c(1:3,6:17)
toc()

# data = readRDS(paste(wdir, '2_pipeline/02_era5_processing/data', var, '.rds', sep = ''))
# simplify dataframe and compute anomalies
if (var %in% c('tp','cape','z500','rh1000','t2m')){
  data = data %>% group_by(lat, lon, week) %>% summarise(mean = mean(value), sd = sd(value)) %>%
    right_join(data, by = c('lat', 'lon', 'week')) %>% mutate(anom = (value-mean)/sd) %>%
    select(lat, lon, week, year, value, anom)
}

# write to rds
saveRDS(data, paste(wdir, '2_pipeline/02_era5_processing/data', var, '.rds', sep = ''))

# optional: clip to MODIS time frame (for faster reading if only parts of data are needed)
datMODIS = filter(data, year > 2000)
saveRDS(datMODIS, paste(wdir, '2_pipeline/02_era5_processing/data', var, '_MODIS.rds', sep = ''))


### combine v250 and u250 to total wind and write out -------------------------------------------------
data = readRDS(paste(wdir, '2_pipeline/02_era5_processing/datav250.rds', sep = '')) %>% 
  rename(v250 = varname) %>% ungroup() %>% filter(lat > 50) %>%
  full_join(readRDS(paste(wdir, '2_pipeline/02_era5_processing/datau250.rds', sep = '')) %>% 
              rename(u250 = varname) %>% ungroup() %>% filter(lat > 50)) %>%
  mutate(value = sqrt(u250^2 + v250^2)) %>% select(lon, lat, week, year, value)

data = data %>% group_by(lat, lon, week) %>% summarise(mean = mean(value), sd = sd(value)) %>%
  right_join(data, by = c('lat', 'lon', 'week')) %>% mutate(anom = (value-mean)/sd) %>%
  select(lat, lon, week, year, value, anom)

saveRDS(data, paste(wdir, '2_pipeline/02_era5_processing/data_totalwind.rds', sep = ''))
datMODIS = filter(data, year > 2000)
saveRDS(datMODIS, paste(wdir, '2_pipeline/02_era5_processing/data_totalwinds_MODIS.rds', sep = ''))

### combine t and rh to vpd and write out ---------------------------------------
data = readRDS(paste(wdir, '2_pipeline/02_era5_processing/datat2m.rds', sep = '')) %>% 
  rename(t2m = value) %>% select(-anom) %>%
  full_join(readRDS(paste(wdir, '2_pipeline/02_era5_processing/datarh1000.rds', sep = '')) %>% 
              rename(r = value)) %>% select(-anom) %>%
  mutate(value = (100 - r) * (610.7 * 10^(7.5 * t2m / (237.3 + t2m)))) 

data = data %>% group_by(lat, lon, week) %>% summarise(mean = mean(value), sd = sd(value)) %>%
  right_join(data, by = c('lat', 'lon', 'week')) %>% mutate(anom = (value-mean)/sd) %>%
  select(lat, lon, week, year, value, anom)
  
saveRDS(data, paste(wdir, '2_pipeline/02_era5_processing/datavpd.rds', sep = ''))
datMODIS = filter(data, year > 2000)
saveRDS(datMODIS, paste(wdir, '2_pipeline/02_era5_processing/datavpd_MODIS.rds', sep = ''))

