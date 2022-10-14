# -*- coding: utf-8 -*-
"""
Aggregates MCD64A1 data from original hdf tiles to 0.25 degree grid
@author: Rebecca Scholten
"""

import numpy as np
import pyproj
import glob
import os
import re
import gdal, gdalconst, osr
from calendar import monthrange, isleap
import datetime
import time

sdir = 'D:/shared_data/' # directory with data shared by projects (modis)
wdir = 'D:/waves/' # wdir on workstation 


#%%
''' Supporting functions: '''

def world2Pixel(gt, Xgeo, Ygeo):
    ''' Uses a geomatrix (gdal.GetGeoTransform()) to calculate 
    the pixel location of a geospatial coordinate'''
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
    nodata value is based on data type
    required inputs: filename, array, projection in epsg code, geotransform, data type
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


def tilemap(sampleloc, mask_earth=True, mres=None):
    '''
    Works like the MODIS tile calculator: https://landweb.modaps.eosdis.nasa.gov/cgi-bin/developer/tilemap.cgi
    Based on function provided by D. van Wees
    Geolocation array is made using 'reverse' mapping and:  index_x, index_y = np.meshgrid(range(2400), range(2400)). (see also the function 'construct_geolocation').
    :param sampleloc: sample location [tile, index_y, index_x]
    :param mask_earth (default=True): mask outside earth as NaN [bool]. mask_earth=False only for debugging purposes.
    :param mres (default=500): for MODIS resolution other than 500m, e.g. 250m (mres=250) or 1000m (mres=1000).
    :return: [lat, lon] in degrees
    '''
    
    if mres is None: mres = 500
    ndim = (500 * 2400) / mres  # number of row and column pixels in MODIS tile.
    
    sphere_radius = 6371007.181
    proj4str = ("+proj=sinu +a=%f +b=%f +units=m" % (sphere_radius, sphere_radius))
    p_modis_grid = pyproj.Proj(proj4str)
    
    R0 = 6371007.181000  # Earth radius in [m]
    limit_left = -20015109.354  # left limit of MODIS grid in [m]
    realres = ((abs(limit_left) * 2) / 36) / ndim  # actual size of each MODIS tile  (alternative: cell_size = ((limit_top*2)/18) / 2400)
    T = ndim * realres  # size of MODIS tile in meters.
     
    tilen = sampleloc[0]
    index_x = sampleloc[2]
    index_y = sampleloc[1]
    
    hn = int(tilen[1:3])
    vn = int(tilen[4:])
    
    lon_frac = (index_x + 0.5) / ndim + hn      # +0.5 to get cell midcenter.
    lat_frac = (index_y + 0.5) / ndim + vn
    
    x = (lon_frac - 36/2) * T
    y = - (lat_frac - 18/2) * T
    
    lon, lat = p_modis_grid(x, y, inverse=True)
    
    
    if mask_earth == True:
        
        phi = np.deg2rad(lat)
        x_border = np.deg2rad(180.0) * R0 * np.cos(phi)
        
        outside_earth = np.abs(x) > x_border
        
        if (type(outside_earth) is np.bool_):
            if outside_earth == True:
                lat = np.nan
                lon = np.nan
                y = np.nan
                x = np.nan
            else: pass
        else:
            lat[outside_earth] = np.nan
            lon[outside_earth] = np.nan
            y[outside_earth] = np.nan
            x[outside_earth] = np.nan
    
    output = [lat, lon]
        
    return output

def construct_geolocation(tilen, mask_earth=True, mres=None):
    '''
    Construct MODIS tile geolocation arrays.
    Original function provided by D. van Wees
    :param tilen: MODIS tile [str]
    :param mask_earth: mask outside earth as NaN [bool]
    :param mres (default=500): for MODIS resolution other than 500m, e.g. 250m (mres=250) or 1000m (mres=1000).
    :return: lats_geo: latitude geolocation array
    :return: lons_geo: longitude geolocation array
    '''
    
    if mres is None: mres = 500
    ndim = int((500 * 2400) / mres)    # number of row and column pixels in MODIS tile.
    
    index_x, index_y = np.meshgrid(range(ndim), range(ndim))
    lats_geo, lons_geo = tilemap([tilen, index_y, index_x], mask_earth=mask_earth, mres=mres)
    
    # No need to save: loading Float64 saved .tif geoloc opens in 0.5 s, just as fast as calculating. Float64 precision is necessary.
    
    return lats_geo, lons_geo

''' End '''

#%%
resdeg = 0.25  # degree resolution to aggregate to.
dtype = 'Int16' # datatype to be save to.

total_pix = np.load(sdir + '00_Total_pixels_025d.npy') # number of MODIS pixels in every degree grid cell

for year in range(2021, 2022):  # years to be processed.
    print(year)
    outname = wdir + '2_pipeline/01_modis/regrid25/MCD64A1C6_' + str(year) + '.tif'
    outname_week = wdir + '2_pipeline/01_modis/regrid25/MCD64A1C6_' + str(year) + '_week.tif'
    
    filelist = glob.glob(os.path.join(sdir + 'MCD64A1v061/', '*.hdf'))  # list of files in folder.
    tilelist = [re.search(r'h..v..', f).group() for f in filelist]  # list of tiles in folder.
    tilelist.sort()
    seen = set()  # removes duplicates from tilelist, perserving order.
    tilelist = [x for x in tilelist if not (x in seen or seen.add(x))]
    
    # check if no tiles are missing
    if len(tilelist) < 41:
        print('tiles are missing!')
    
    tdim = 365 + isleap(year)*1  # time dimension size. Set according to time dimension in input file.
    
    summ = np.zeros((tdim, 180 * int(1 / resdeg), 360 * int(1 / resdeg)))   # 0.25 grid where aggregated result is stored in.
    
    
    for tilen in tilelist:
        t0 = time.time()
        print('Processing tile: ' + tilen)
        
        ''' Creating aggregation matrix '''
        lats_geo, lons_geo = construct_geolocation(tilen, mask_earth=True, mres=500) # load geolocation arrays.
        outside_earth = np.isnan(lons_geo)  # pixels outside the real globe (in case of h10v02 for example).
        lats_geo[outside_earth] = 999.999
        lons_geo[outside_earth] = 999.999
        lat_index = np.floor(np.abs(lats_geo - 90.0) * int(1 / float(resdeg))).astype(int)  # <- This is where the magic happens. Be sure to understand this part.
        lon_index = np.floor((lons_geo + 180.0) * int(1 / float(resdeg))).astype(int)
        ''' End '''
        
        filelist_tile = [file for file in filelist if tilen in file and 'A'+str(year) in file]
        
        # check if months are missing
        if year == 2021:
            if len(filelist_tile) < 8:
                print('months are missing!' + str(len(filelist_tile)))
        else:
            if len(filelist_tile) < 12:
                print('months are missing!')
        
        for filename in filelist_tile:
            # print(filename)
            ''' Load input data to be aggregated (e.g. MODIS tif or hdf) '''
            ds = gdal.Open('HDF4_EOS:EOS_GRID:' + filename + ':MOD_Grid_Monthly_500m_DB_BA:Burn Date', gdalconst.GA_ReadOnly)
            Data = ds.ReadAsArray().astype(float)
            NoDataValue = ds.GetRasterBand(1).GetNoDataValue()
            
            Data[Data == NoDataValue] = np.nan
            
            if Data.ndim == 2: Data = np.expand_dims(Data, axis=0)
            ''' End '''
            
            Data[:, outside_earth] = np.nan  # Overdone, but better to be sure.
            nanarr = (np.isnan(Data) & ~outside_earth)  # .astype(int)  # nan values inside earth boundaries.
            Data[nanarr] = 0  # IMPORTANT: in case of arr[:,x,y] = [np.nan, 0.35]. The case e.g. if biome mask changes over years.
            
            # calculate start and end doy of that month for this year
            datestring = filename.split('.')[1][1:] # fetch date from filename
            mon = datetime.datetime.strptime(datestring, '%Y%j').month # extract month from Julian day
            totaldays = monthrange(year, mon)[1] # this returns number of days for each month
            doy_start = int(datestring[-3:]) # this is the first Julian day of the month
            doy_end = doy_start + totaldays # this already includes the cut off day
            idx_doy_start = doy_start - 1
            idx_doy_end = doy_end - 1
            
            # spread the Data to daily fire yes/no
            dataday = np.array([(Data == day)*1 for day in range(doy_start, doy_end)]).squeeze()
            
            
            ''' Main aggregation algorithm! '''
            # if no burned area happened in the dataset, we can just keep the zeroes
            if np.sum(dataday) > 0:
                index = np.where(np.nansum(Data, axis=0) > 0)  # 500m pixels with zero data or nan can be skipped. Taking the sum to see if its zero for every time step.
                if dataday.shape[0] > 1:  # much faster for multiple time steps.
                    for x, y in zip(index[0], index[1]):
                        summ[idx_doy_start:idx_doy_end, lat_index[x, y], lon_index[x, y]] = summ[idx_doy_start:idx_doy_end, lat_index[x, y], lon_index[x, y]] + dataday[:, x, y]
            
                elif dataday.shape[0] == 1:  # Faster for 1 time step. Indexing with 0 instead of ':'.
                    for x, y in zip(index[0], index[1]):
                        summ[0, lat_index[x, y], lon_index[x, y]] = summ[0, lat_index[x, y], lon_index[x, y]] + dataday[0, x, y]
        print('Processed in:', time.time() - t0, 's')
    
    ''' End '''
    
    X = summ / total_pix.astype('float')
    X = np.squeeze(np.array(X))
    
    ## clip to >50N
    X = X[:, :160,:]
    
    ## create weekly dataset
    # extract startday (March 2)
    startday = int(datetime.datetime(year, 3, 2, 12).strftime('%j')) # start with March 2
    
    # create sequence for looping
    # we need to substract 1 from arrays since python starts counting at 0
    lower = np.arange(startday, startday + (30*7), 7)-1
    upper = np.arange(startday + 7, startday + (31*7), 7)-1
    
    # sum weeks
    X_week = np.array([np.sum(X[lower[week]:upper[week],:,:], axis = 0) for week in range(30)])

    ''' Save data '''
    
    dsGEO = (-180.0, 0.25, 0.0, 90.0, 0.0, -0.25)
    metadata = ds.GetMetadata()
    
    # daily
    writeTif(outname, X, 4326, dsGEO, dtype, metadata)
    
    # weekly
    writeTif(outname_week, X_week, 4326, dsGEO, dtype, metadata)