# -*- coding: utf-8 -*-
"""
Clipping of fire data to the larch forests of Siberia
Computation of each grid point's straight-line distance from the treeline
Weekly summary of fire activity above and below the tree line and within the Arctic circle

@author: Rebecca Scholten
"""

import sys
import math
import numpy as np
from glob import glob
import pandas as pd
from scipy import ndimage
import gdal,ogr

wd = 'D:/waves/' # wdir

#%% functions

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

def world2Pixel(gt, Xgeo, Ygeo):
    ''' Uses a geomatrix (gdal.GetGeoTransform()) to calculate the pixel location of a geospatial coordinate'''
    gt = list(gt)
    Xpx = int((Xgeo - gt[0]) / gt[1])
    Ypx = int((Ygeo - gt[3]) / gt[5])
    return (Xpx, Ypx)

def geo_to_polar(lon_arr, lat_arr):
    '''transform lists of geographic lat lon coordinates to polar LAEA grid (projected)'''
    
    import numpy as np
    import pyproj
    
    proj4str = ("epsg:3571")
    p_modis_grid = pyproj.Proj(proj4str)
        
    x_arr, y_arr = p_modis_grid(lon_arr, lat_arr)
    x_arr = np.array(x_arr)
    y_arr = np.array(y_arr)
    
    return x_arr, y_arr

def get_coords(arr, lats, lons, proj=True):
    coords = []
    for i in range(len(arr[0])):
       lat = lats[arr[0][i]]
       lon = lons[arr[1][i]]
       coords.append((lon,lat))
    if proj == True:
        lons,lats = zip(*coords)
        x, y = geo_to_polar(lons,lats) # reproject point to projected crs
        coords = list(zip(x,y))
    
    return coords

def build_rtree(coords):
    '''Builds Rtree from a shapely multipolygon shape
    and optionally uses list of fids as identifier'''
    import rtree
    
    idx = rtree.index.Index() # create new index
    for ind, coord in enumerate(coords):
        # left, bottom, right, top
        idx.insert(ind, (coord[0],coord[1],coord[0],coord[1]), ind)
    
    return idx

#%% set resolutionm and lons/lats
lats = np.arange(90- 0.25 / 2, 50, -0.25)
lons = np.arange(-180 + 0.25 / 2, 180, 0.25)
lons_clip = lons[lons > 110]

lats_flat = np.repeat(lats, len(lons_clip))
lons_flat = np.tile(lons_clip, int(len(lats_flat)/len(lons_clip)))

#%% create an image with the same resolution and extent as gfed/regridded MODIS data
# where boreal forest = 1, arctic/tundra = 2, and everything else = 0
# this can be used for summing up regions

# as input we here take the larch ecotone for boreal siberia
# and the circumpolar vegetation map for the treeline
ds_boreal, lyr_boreal = readDS(wd + '2_pipeline/00_ecoregions/eastern_siberia.shp')
ds_polar, lyr_polar = readDS(wd + '0_data/ecoregions/cav/cav_wgs84.shp')

outname = wd + '2_pipeline/00_ecoregions/temp/base_larch.tif'

# rasterize first shapefile into numpy array    
baseraster = layer2raster(lyr_boreal, 0.25, outname, (-180, 50, 180, 90))
base = baseraster.GetRasterBand(1).ReadAsArray()  
gt = baseraster.GetGeoTransform()
# base_clean = morphology.remove_small_objects(base == 255 , 6)
baseraster = None

# rasterize second shapefile into numpy array and merge with first
baseraster = layer2raster(lyr_polar, 0.25, outname, (-180, 50, 180, 90))
base1 = baseraster.GetRasterBand(1).ReadAsArray()  
baseraster = None
base[(base1 != 0) & (base == 255) ] = 2

# just to be safe, areas north of northern tundra area assigned to tundra as well
also_tundra, maxlab = ndimage.label(base == 255)
also_tundra[also_tundra == 16] = 0
base[also_tundra > 0] = 2

# write out the new raster to csv using lat and lon files
lats_flat = np.repeat(lats[lats > 50], len(lons))
lons_flat = np.tile(lons, int(len(lats_flat)/len(lons)))

#%% compute straight-line distance between pixel centroids and treeline

# insert all tundra pixels into rtree index
tundra = np.where(base == 2)
coords_tundra = get_coords(tundra, lats, lons, proj=True)
idx = build_rtree(coords_tundra)

# loop through boreal pixel, find nearest 
boreal = np.where(base == 255)
coords_taiga = get_coords(boreal, lats, lons, proj=True)
coords_geo = get_coords(boreal, lats, lons, proj=False)

res = []
for ind,coord in enumerate(coords_taiga):
    neighbours = list(idx.nearest(coord, 1, objects='raw'))
    if len(neighbours) > 1:
        print('more than one neighbour')
    tl_coord = coords_tundra[neighbours[0]]
    dist = math.sqrt((coord[0]-tl_coord[0])**2 + (coord[1]-tl_coord[1])**2)
    res.append((coords_geo[ind][0], coords_geo[ind][1], dist))
df = pd.DataFrame(res, columns=['lon','lat','dist_tl'])
df.to_csv(wd + '2_pipeline/01_modis/distance_treeline_centroids.csv')

#%% write weekly burned area, larch ecotone boreal and tundra to tidy dataframe

filepath = '2_pipeline/01_modis/regrid25/70N/*week*.tif'
outname = '2_pipeline/01_modis/perc_ba_0.25deg_week_70N_larch.csv'

# collect filenames
filenames = glob(wd + filepath)#[-2:]

# read number of pixels and reshape 
total_pix = np.load(wd + '2_pipeline/01_modis/regrid25/total_pixels_25.npy')
if total_pix.shape[0] > 160:
    total_pix = total_pix[:160,:]
total_pix = total_pix.flatten(order = 'C')

# read burned area tiffs
year = 2001 # start year
dfs = []
for file in filenames:
    year_arr = np.repeat(year, len(lats_flat))
    for week in range(1, 31):
        # create array for week number
        week_arr = np.repeat(week, len(lats_flat))
        # read percentage burned dataset
        temp = readTif(file, week, return_ds = False)
        for val in [2,255]:
            # set all pixels outside of region to 0 and reshape
            mask = base == val
            temp0 = (temp*mask).flatten(order = 'C')
            
            # write to pd_dataframe
            df_temp = pd.DataFrame({'lat':lats_flat, 'lon':lons_flat, 
                                    'reg':np.repeat(val, len(lats_flat)),
                                    'perc_burned':temp0, 'n_pix': total_pix,
                                    'year':year_arr, 'week':week_arr})
            # filter for perc burned larger 0
            dfs.append(df_temp[df_temp['perc_burned'] > 0])
            
    year += 1

#concatenate dataframes
df = pd.concat(dfs)
df.to_csv(wd + outname)

#%% same sums for Arctic circle

filepath = '2_pipeline/01_modis/regrid25/70N/*week*.tif'
outname = '2_pipeline/01_modis/perc_ba_0.25deg_week_70N_AC.csv'

# read dataset
filenames = glob(wd + filepath)#[-2:]

# read number of pixels and reshape 
total_pix = np.load(wd + '2_pipeline/01_modis/regrid25/total_pixels_25.npy')
if total_pix.shape[0] > 160:
    total_pix = total_pix[:160,:]
total_pix = total_pix[lats > 66.5635, :].flatten(order = 'C')

# clip lats to Arctic Circle
lats_AC = lats[lats > 66.5635]
lats_AC_flat = np.repeat(lats_AC, len(lons))
lons_AC_flat = np.tile(lons, int(len(lats_AC_flat)/len(lons)))

# read burned area tiffs
year = 2001 # start year
dfs = []
for file in filenames:
    year_arr = np.repeat(year, len(lats_AC_flat))
    for week in range(1, 31):
        # create array for week number
        week_arr = np.repeat(week, len(lats_AC_flat))
        # read percentage burned lats_AC_flat
        temp = readTif(file, week, return_ds = False)
        
        # set all pixels outside of region to 0 and reshape
        temp0 = temp[lats > 66.5635, :].flatten(order = 'C')
            
        # write to pd_dataframe
        df_temp = pd.DataFrame({'lat':lats_AC_flat, 'lon':lons_AC_flat,
                                'perc_burned':temp0, 'n_pix': total_pix,
                                'year':year_arr, 'week':week_arr})
        # filter for perc burned larger 0
        dfs.append(df_temp[df_temp['perc_burned'] > 0])
            
    year += 1

#concatenate dataframes
df = pd.concat(dfs)
df.to_csv(wd + outname)
