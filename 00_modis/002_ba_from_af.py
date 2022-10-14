# -*- coding: utf-8 -*-
"""
Uses linear regression to approximate MODIS burned area from 
MODIS active fire counts within 0.25 degree grid cells
for areas North of 70N
@author: Rebecca Scholten
"""

wd = 'D:/waves/' # work directory of current project
sdir = 'D:/shared_data/' # folder with data shared between projects

import numpy as np
import datetime
import gdal, ogr, osr
from glob import glob
import time
from scipy import stats

'''
prerequisites:
    active fires 2002 - 2021
    burned area 2002 - 2021
'''

#%% function for reading mcd14 txt files
def world2Pixel(gt, Xgeo, Ygeo):
    ''' Uses a geomatrix (gdal.GetGeoTransform()) to calculate the pixel location of a geospatial coordinate'''
    gt = list(gt)
    Xpx = int((Xgeo - gt[0]) / gt[1])
    Ypx = int((Ygeo - gt[3]) / gt[5])
    return (Xpx, Ypx)

def readTif(filename, band, lims = None, return_ds = True):
    '''reads a tif raster band into numpy array using gdal,
    optionally the bounds can be specified in this order:
    (ulX, ulY, lrX, lrY)
    and the dataset or only the array can be returned'''
    ds0 = gdal.Open(filename)
    bd0 = ds0.GetRasterBand(band)
    if lims:
        minX, maxY, maxX, minY = lims
        geoTransImg = ds0.GetGeoTransform()
        ulX, ulY = world2Pixel(geoTransImg, minX, maxY)
        lrX, lrY = world2Pixel(geoTransImg, maxX, minY)
        arr0 = bd0.ReadAsArray(ulX, ulY, lrX - ulX, lrY - ulY)
    else:
        arr0 = bd0.ReadAsArray()
    if return_ds:
        return ds0, bd0, arr0
    else:
        ds0 = bd0 = None
        return arr0

def writeTif(filename, arr, epsg, geotrans, dtype, metadata = None, scale_factor = None):
    '''writes a 3D array to geotiff, first dimension should be time
    automativally creates a nodata value based on data type
    compulsory inputs: filename, array, projection in epsg code, geotransform, data type
    optional inputs: metadata, scale-factor'''
    cols = arr.shape[2]
    rows = arr.shape[1]
    bands = len(arr)
        
    # set projection and geotransform
    dsSRS = osr.SpatialReference()
    dsSRS.ImportFromEPSG(epsg)
        
    # set data type and nodatavalue
    if dtype in ['Byte', 'Int16']:
        arr = np.round(arr)
    
    if dtype == 'Byte':         NoDataValue = np.iinfo(np.uint8).max
    elif dtype == 'Int16':      NoDataValue = np.iinfo(np.int16).max
    elif dtype == 'Float32':    NoDataValue = -1   # Max float32 doesnt work, and danger of rounding errors: np.finfo(np.float32).max
    arr[np.isnan(arr)] = NoDataValue
    
    if dtype == 'Byte':         gdal_dtype = gdal.GDT_Byte  # Alternatively use: gdal_array.NumericTypeCodeToGDALTypeCode()
    elif dtype == 'Int16':      gdal_dtype = gdal.GDT_Int16
    elif dtype == 'Float32':    gdal_dtype = gdal.GDT_Float32
    
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(filename, cols, rows, bands, gdal_dtype, options=['COMPRESS=LZW', 'INTERLEAVE=BAND', 'TILED=YES'])
    ds.SetGeoTransform(geotrans)
    ds.SetProjection(dsSRS.ExportToWkt())
    
    if metadata is not None:
        ds.SetMetadata(metadata)
    
    for i in range(bands):
        band = ds.GetRasterBand(i + 1)
        if NoDataValue:
            band.SetNoDataValue(NoDataValue)
        if scale_factor is not None:
            band.SetScale(scale_factor)
        band.WriteArray(arr[i, :, :])
    
    ds = None   # close file.
    return None

def read_MCD14ML_clip(fnmFC, ext = False):
    ''' Read monthly MCD14 fire locations (txt files)
    '''
    import pandas as pd
    # read and extract
    usecols = ['YYYYMMDD','HHMM','sat','lat','lon','T21', 'T31','sample','FRP','conf','type','dn']
    # read dataframe
    df = pd.read_csv(fnmFC,delim_whitespace=True,skiprows=1,header=None,
                     names = usecols,parse_dates=['YYYYMMDD'])

    # spatial filter and quality filter
    if ext:
        out = df.loc[(df['lat'] > ext[1]) & (df['lat'] < ext[3]) &
                     (df['lon'] > ext[0]) & (df['lon'] < ext[2]) &
                     (df['type'] == 0) & (df['conf'] > 29)]  
    else:
        out = df.loc[(df['type'] == 0) & (df['conf'] > 29)]
    
    # convert date in doy and filter on day
    year = int(str.split(file,'.')[-4][:4])
    out['doy'] = out['YYYYMMDD'].dt.dayofyear
    startday = int(datetime.datetime(year, 3, 2, 12).strftime('%j')) # start with March 1
    endday = startday + 7*30 # end day is 30 weeks later
    out = out.loc[(out['doy'] >= startday) & (out['doy'] < endday)] # endday includes cutoff day!
    
    # drop unused columns
    out = out.drop(out.columns[list(range(3))+list(range(5,12))], axis=1)
    
    return out


def addAF2grid(arr, lon_grid, lat_grid, date_index):
    '''takes lists of lat, lon and doy indices and adds them to a
    array containing active fire counts'''
    for i in range(len(lon_grid)): # loop through all fire pixels
        if len(lon_grid[i]) > 1 or len(lat_grid[i]) > 1:
            # if a fire pixel is at the border we assign a fraction of the value to all fields evenly
            # i.e. if it is at the border between 2 pixels, both get the value 0.5
            val = 1 / (len(lon_grid[i]) * len(lat_grid[i]))
            if val == 0.5:
                if len(lon_grid[i]) > 1:
                    arr[date_index, lat_grid[i], [lon_grid[i][0], lon_grid[i][1]]] += val
                else:
                    arr[date_index, [lat_grid[i][0], lat_grid[i][1]], lon_grid[0]] += val
            else:
                arr[date_index, [lat_grid[i][0], lat_grid[i][1]], [lon_grid[i][0], lon_grid[i][1]]] += val
        else: 
            # add a value of 1 to the grid_cell
            arr[date_index[i], lat_grid[i], lon_grid[i]] += 1
    return arr

#%% specs
endyear = 2021          # last year to be processed
nyear = endyear-2000    # number of years to be processed (starts with beginning of MODIS era)

#%% read yearly burned area for 2000 - 2020

filenames = glob(wd + '2_pipeline/01_modis/regrid25/*week.tif')
# read grid cell area: count of modis pixels times sixe of a modis pixel
grid_area = np.load(sdir + '00_Total_pixels_025d.npy')[:160,:] * 463.312716528**2

# read burned area
BA_boreal = []
for file in filenames:
    BA_boreal_yr = []
    numweeks = 30
    for week in range(1, numweeks + 1):
        temp = readTif(file, week, return_ds = False)
        BA_boreal_yr.append(temp)
    BA_boreal_yr = np.array(BA_boreal_yr)
    # flatten to yearly
    BA_boreal.append(np.sum(BA_boreal_yr, axis = 0))
    
BA_boreal = np.array(BA_boreal)


#%% read active fires (this takes about 25 minutes)
# loop through active fires and write number of fires to a 0.25 degree grid
# we only take active fires in the main season (between March 1 and 30 weeks after/somewhere in September)

# create an ouptput grid
AF_boreal = np.empty((7*30*nyear, 160, 1440), dtype=float) # daily, 160 lats, 1440 lons
lats = np.arange(90- 0.25 / 2, 50, -0.25)
lons = np.arange(-180 + 0.25 / 2, 180, 0.25)    

# read shapefile (C6.1 only) or txt (collection 6)
t0 = time.time()

filenames = glob(sdir + 'MCD14MLv006_txt/*.txt') # list of ML txt files
for file in filenames:
    t1 = time.time()
    year = int(str.split(file,'.')[-4][:4])
    mon = int(str.split(file,'.')[-4][4:])
    startday = int(datetime.datetime(year, 3, 2, 12).strftime('%j')) # start with March 1
    if (year < 2001) | (mon < 3) | (mon > 9): # skip 2000, and ONDJF
        continue
    
    # read the file
    df_data = read_MCD14ML_clip(file, ext = [-180, 50, 180, 90])
    
    # list of grid indices from lon and lat coordinates
    lon_grid = [np.where(abs(lons-lon) == min(abs(lons-lon))) for lon in df_data['lon']]
    lat_grid = [np.where(abs(lats-lat) == min(abs(lats-lat))) for lat in df_data['lat']]
    
    # list of doy-year-indices from year and doy
    year_index = year - 2001
    startday = int(datetime.datetime(year, 3, 2, 12).strftime('%j'))
    doy_ind = [year_index*7*30 + doy - startday for doy in df_data['doy']]
    
    # add to active fires array
    AF_boreal = addAF2grid(AF_boreal, lon_grid, lat_grid, doy_ind)
    
    print(file, ',', str(len(df_data)), 'active fires processed in', str(time.time()-t1))
    
print('Read in', round(time.time()-t0, 2), 's')
   
#%% create a linear model predicting the burned area from the number of fire pixels

# first we flatten both datasets to 1d using the same order
BA_train_flat = BA_boreal[:,80:,:].flatten(order = 'C')
# we use yearly data for the aggregation
lower = np.arange(0, 30*7*nyear, 30*7) # create sequence for looping
upper = np.arange(30*7, 30*7*(nyear+1), 30*7)
AF_train_year = np.array([np.sum(AF_boreal[lower[year]:upper[year],80:,:], axis = 0) for year in range(nyear)]).flatten(order = 'C')

# build linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(AF_train_year, BA_train_flat)

## predict daily burned area from active fires for north of 70N
AF_pred_flat = AF_boreal[:,:80,:].flatten(order = 'C')
BA_pred_flat = AF_pred_flat * slope + intercept
BA_pred_flat[AF_pred_flat == 0] = 0

# reshape back to time, lat, lon
BA_pred = BA_pred_flat.reshape((nyear*30*7, 80, len(lons)), order = 'C')

#%% merge with MCD64 data and write out

lower = np.arange(0, 30*7*nyear, 30*7) # create sequence for looping
upper = np.arange(30*7, 30*7*(nyear+1), 30*7)
for yr in range(nyear):
    # read daily burned area into np array
    file = wd + '2_pipeline/01_modis/regrid25/MCD64A1C6_' + str(yr + 2001) + '.tif'
    BA_old = []
    day1 = int(datetime.datetime(yr + 2001, 3, 2, 12).strftime('%j'))
    
    for yday in range(day1, 30*7 + day1):
        temp = readTif(file, yday, return_ds = False)
        BA_old.append(temp)
    BA_old = np.array(BA_old)[:,80:,:]

    # extract year from prediction dataset
    BA_new = BA_pred[lower[yr]:upper[yr],:,:]
    # attach the predicted files
    BA_out = np.hstack((BA_new, BA_old))

    # sum weeks
    # create sequence for looping
    lowerwk = np.arange(0, 30*7, 7)
    upperwk = np.arange(7, 31*7, 7)
    BA_out_week = np.array([np.sum(BA_out[lowerwk[week]:upperwk[week],:,:], axis = 0) for week in range(30)])

    # write to file
    dsGEO = (-180.0, 0.25, 0.0, 90.0, 0.0, -0.25)
    metadata = gdal.Open(file).GetMetadata()
    dtype = 'Float32'
    outname = wd + '2_pipeline/01_modis/regrid25/70N/MCD64A1C6_70n' + str(yr + 2001) + '.tif'
    outname_week = wd + '2_pipeline/01_modis/regrid25/70N/MCD64A1C6_70n_week' + str(yr + 2001) + '.tif'
    
    writeTif(outname, BA_out, 4326, dsGEO, dtype, metadata)
    writeTif(outname_week, BA_out_week, 4326, dsGEO, dtype, metadata)