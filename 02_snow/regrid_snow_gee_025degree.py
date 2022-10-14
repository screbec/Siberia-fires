# -*- coding: utf-8 -*-
"""
Regrid snowmelt data to 0.25 deg
@author: Reecca Scholten
"""
import gdal

path_to_file = 'D:/waves/0_data/snow/gee_ldos/'

infile = path_to_file + 'ldos_larch_15_2021_mosaic_clip.tif'
outfile = path_to_file + 'ldos_larch_15_2021_mosaic_clip_025deg.tif'

# take bounds from burned area raster
extentraster = r'D:\waves\2_pipeline\01_modis\regrid25\MCD64A1C6_2020_week.tif'
src = gdal.Open(extentraster)
minX, xres, xskew, maxY, yskew, yres  = src.GetGeoTransform()
maxX = minX + (src.RasterXSize * xres)
minY = maxY + (src.RasterYSize * yres)
src = None

# warp to polar LAEA
warp = gdal.Warp(outfile,infile,dstSRS='EPSG:4326',
                 xRes=0.25, yRes=0.25, resampleAlg = 'average',
                 outputBounds=[minX, minY, maxX, maxY], 
                 outputBoundsSRS='EPSG:4326')
warp = None