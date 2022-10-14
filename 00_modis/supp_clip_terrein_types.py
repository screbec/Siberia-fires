# -*- coding: utf-8 -*-
"""
Convert terrain shapefile to 0.05 degree res raster

@author: Rebecca Scholten
"""

import sys
import numpy as np
import gdal, ogr
import pandas as pd
import geopandas as gpd

wd = 'D:/waves/' # wdir


#%% functions
def readDS(filename, drvstring = 'ESRI Shapefile'):
    try:
        driver = ogr.GetDriverByName(drvstring)
    except:
        print('Driver not found.')
        sys.exit(1)
    ds = driver.Open(filename)
    lyr = ds.GetLayer()
    return(ds, lyr)

def layer2raster(inlyr, res, outname, extent = None, inattribute = None):
    '''takes a shapefile layer and transforms it to a raster with a 
    specified resolution and output name
    optionally the layer extent can be set and one attribute can be saved as values in the raster'''
    if extent:
        minx, miny, maxx, maxy = extent
    else:
        # get ulx and uly from extent
        minx, maxx, miny, maxy = inlyr.GetExtent() 
    # calculate number of x and y pixels from extent and resolution
    xpx = int(np.ceil((maxx - minx)/res))
    ypx = int(np.ceil((maxy - miny)/res))
    # create the geotransform
    geotrans = (np.floor(minx), res, 0, np.ceil(maxy), 0, -res) 
    # create output datasource
    driver = gdal.GetDriverByName('GTiff')
    outraster = driver.Create(outname, xpx, ypx, 1, gdal.GDT_Byte)
    outraster.SetGeoTransform(geotrans)
    outraster.SetProjection(inlyr.GetSpatialRef().ExportToWkt())
    if inattribute:
        err = gdal.RasterizeLayer(outraster, [1], inlyr, None, options=["ATTRIBUTE=%s" % inattribute])
    else:
        err = gdal.RasterizeLayer(outraster, [1], inlyr, None)
    if err != 0:
        raise Exception("error rasterizing layer: %s" % err)
    return outraster

#%% change terrain types shapefile to have numeric column for burning
gdf = gpd.read_file(r'D:\shared_data\fedorov\MLK_2019\commondata\shp\tm_veg_union.shp')
gdf = gdf.assign(TM_ID_num = gdf.TM_ID)

# remove some entries
gdf = gdf[gdf.TM_ID != '0'] # water bodies
gdf = gdf[gdf.TM_ID != 'XVIII'] # Tukulan (sand dunes)
gdf = gdf[gdf.TM_ID != 'XVII'] # Glaciers
gdf = gdf[gdf.TM_ID != 'XI'] # Rock-ridge

gdf['TM_ID_num'].replace(list(set(gdf.TM_ID)), list(range(1,18)), inplace=True)
gdf_4326 = gdf.to_crs('EPSG:4326')
gdf_4326.to_file(wd + '2_pipeline/00_ecoregions/tm_veg_union.shp')

lookup = gdf.filter(regex='TM_ID', axis=1).drop_duplicates()
lookup.to_csv(wd + '2_pipeline/00_ecoregions/terrain_types_lookup.csv')

#%% create an image with the same resolution and extent as regridded MODIS data

outname = wd + '2_pipeline/00_ecoregions/terrain_types_5.tif'
# read shapefile
ds_boreal, lyr_boreal = readDS(wd + '2_pipeline/00_ecoregions/tm_veg_union.shp')

# rasterize first shapefile into numpy array    
baseraster = layer2raster(lyr_boreal, 0.05, outname, (100, 55, 165, 77), inattribute = 'TM_ID_num')
base = baseraster.GetRasterBand(1).ReadAsArray()
baseraster = None

# warp to 0.25 degree, (majority)
outfile = wd + '2_pipeline/00_ecoregions/terrain_types_25.tif'
warp = gdal.Warp(outfile,outname,dstSRS='EPSG:4326',
                 xRes=0.25, yRes=0.25, resampleAlg = 'mode')
dat = warp.ReadAsArray()
warp = None

# write to csv
lats = np.arange(77- 0.25 / 2, 55, -0.25)
lons = np.arange(100 + 0.25 / 2, 165, 0.25)

lats_flat = np.repeat(lats, len(lons))
lons_flat = np.tile(lons, int(len(lats_flat)/len(lons)))

flat = dat.flatten(order = 'C')
df = pd.DataFrame({'lat':lats_flat, 'lon':lons_flat, 'terrain':flat})
df.to_csv(wd + '2_pipeline/00_ecoregions/terrain_types_25.csv')

