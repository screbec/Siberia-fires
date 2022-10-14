# -*- coding: utf-8 -*-
"""
Write data from Xu et al. 2022 to csv file with fire locations'
Created on Tue Apr 26 10:09:41 2022

@author: Rebecca Scholten
"""

import gdal
import numpy as np
import pandas as pd
import geopandas as gpd

path = 'D:/shared_data/xu_erl_2022/' # path to data from Xu et al 2022

def readTif(file):
    ds = gdal.Open(file)
    arr = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    return gt, arr

def pixel2World(gt, Xpixel, Ypixel):
    '''Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the
    geospatial coordinate of a pixel location'''
    
    Xgeo = gt[0] + Xpixel*gt[1] + Ypixel*gt[2]
    Ygeo = gt[3] + Xpixel*gt[4] + Ypixel*gt[5]
    
    return (Xgeo, Ygeo)

def arr2coordlist(gt, arr):
    mask_mon = arr > 0
    index = np.where(mask_mon)
    coords = []
    for x, y in zip(index[0], index[1]):
        coords.append((x, y, arr[x, y]))
    x_mon, y_mon, doy_mon = zip(*coords)
    
    # transform pixel coordinates to real-world (3576)
    coords = [] # reset coords list
    for px in range(len(x_mon)):
        coords.append(pixel2World(gt, x_mon[px], y_mon[px]))
    x_geo, y_geo = zip(*coords)
        
    return x_geo, y_geo, doy_mon

#%% main
x_out = []
y_out = []
doys_out = []
years = []
for year in range(2012,2021):
    
    file = path + 'Daily_burned_area_' + str(year) + '.tif'
    gt, arr = readTif(file)
    
    # read burned area
    x_geo, y_geo, doy_mon = arr2coordlist(gt, arr)
    
    x_out.extend(x_geo)
    y_out.extend(y_geo)
    doys_out.extend(list(doy_mon))
    years.extend([year]*len(x_geo))

    print('year processed')

df = pd.DataFrame({'x': x_out, 'y': y_out, 'doy': doys_out, 'year':years})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y), crs=3576)
gdf_4326 = gdf.to_crs(crs=4326)
df['lon'] = gdf_4326['geometry'].x
df['lat'] = gdf_4326['geometry'].y
df.to_csv('D:/waves/2_pipeline/01_modis/wenxuan_erl_ba.csv')
