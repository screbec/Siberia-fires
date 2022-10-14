# Siberian fires analysis




Python and R code for the analysis of drivers of Siberian fire extremes

This code uses satellite and reanalysis data to identify weeks with extreme fire activity in the larch forests of Siberia and subsequently analyse drivers of these extreme fires. The results of this work are in review (Scholten et al., 2022). The datasets required for the code are listed in the data requirements and linked in the data availability statement of Scholten et al. (2022). 
  
The scripts include the preprocessing stepsand the analysis of the drivers as shown in the workflow below. The source code for the display items in the publication is also included. Preprocessing was done in python, whereas analyses and plots were coded in R. For easier readability functions used in each script are pasted directly in the corresponding script. If you have any comments or suggestions regarding the code, please share them with us. Also feel free to contact us if you have any questions about the code, data or the analysis in general.

### Workflow: 

#### 1. Preprocessing:
- Regrid MCD64A1 data (hdf tiles) to 0.25 degree grid (00_modis/001_process_mcd64.py) 
- Supplement the burned area data with active fires for regions north of 70N (00_modis /002_ba_from_af.py)
- Clip the gridded fire product to the Siberian larch forests and compute weekly fire summary statistics (00_modis/003_clip2nesib.py)
- Aggregate weather data to weekly averages, compute additional variables and anomalies  (01_weather/020_weekly_fwi.py and 020_weekly_weather.R)
- Prepare snow data for analysis:
    - Compute first snow-free day from NSIDC data and compute anomalies (02_snow/050_snow_nsidc.R)
    - Aggregate first snow-free days from MXD10 (processed in GEE using doi:10.5281/zenodo.596556) to 0.25 degree grid (02_snow/regrid_snow_gee_025degree.py)
    - Prepare MODIS snowmelt data for analysis and compute anomalies (02_snow/051_snow_gee.R)
- Preprocess fire and lightning datasets (includes filtering, computation of climatologies and anomalies) (03_analysis/990_siberia_load.R)
#### 2. Analysis of drivers:
- Code for fire climatology and study area plots in Fig.1 and Fig.S8 (03_analysis/994_siberia_fireclim.R)
- Code for Fig.2 and Fig. S5 (03_analysis/991_figure2.py)
- Analyse compound effect of snow melt and Arctic front jet on fire activity as well as the trends in the drivers (includes code for Fig.3, Fig.4, Fig.S4 and Fig.S7) (03_analysis/993_drivers.R)\
- Code for fire stats, Table S1, and Figs S1 and S9 (03_analysis/992_siberia_fire.R)
- Code for Fig. S2 (03_analysis/997_figureS3.py)
- Code for Fig. S6 (03_analysis/996_figureS6.R)
#### 3. Supplementary analyses
- Extract average FRP from MCD14 for the study area (00_modis/frp_modis.R)
- Extract burned area data from Xu et al. (2022) for comparison with the product used in this study (00_modis/supp_compare_xu2022.py)
- Aggregate terrain types from Fedorov et al. (2018) to assess whether the influence of drivers is modulated by terrain types (00_modis/supp_clip_terrain_types.py)

### Data requirements:
- MODIS burned area (MCD64A1), 2001-2021
- MODIS active fire locations (MCD14), 2001-2021
- MODIS snow cover product (MXD10), accessed via GEE
- ERA5 reanalysis data: 250hPa wind, 500hPa geopotential, fire weather index, 2m-temperature,1000hPa relative humidity, total precipitation, convective available potential energy
- lightning strikes (VAISALA network)

### Software requirements:
Python code tested with Python 3.7.8, required packages: gdal, pyproj, pandas, geopandas, scipy
R code tested with R 4.2.1
