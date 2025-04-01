# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
Title: Amazon TROPOMI SIF 2023 anomaly maps
Author: Gerbrand Koren
Edited by Carlijn Jalink
Data needed: SIF data, amazon mask, 
"""

#%%################################################################################################################################################
#                                                                                                                                                 #
# INITIALIZE                                                                                                                                      #
#                                                                                                                                                 #
###################################################################################################################################################

# -- IPython settings
#%reset -f
#%matplotlib inline
#!conda install netCDF4

# -- Import statements
import pandas as pd  # Zorg dat pandas is geïmporteerd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as image
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import sys
import os


# -- Output settings
save = False

#%%################################################################################################################################################
#                                                                                                                                                 #
# INPUT                                                                                                                                           #
#                                                                                                                                                 #
################################################################################################################################################### 

# -- Coordinates
lat_min_coord = -30 # degN
lat_max_coord = 20 # degN
lon_min_coord = -90 # degE
lon_max_coord = -30 # degE

# -- Dimensions
ntime_year = 6 # -
ntime_month = 12 # -
nlat = 2275 # -
nlon = 2731 # -

# -- Directory
path = r"C:\Users\Jalin005\OneDrive - Universiteit Utrecht\Documents\01 Python Scripts\Data\SIF_test_data\\"
pathEC = r"C:\Users\Jalin005\OneDrive - Universiteit Utrecht\Documents\01 Python Scripts\Data\EC_data\\"


#%%################################################################################################################################################
#                                                                                                                                                 #
# DATA IMPORT                                                                                                                                     #
#                                                                                                                                                 #
###################################################################################################################################################
# TROPOMI SIF DATA

# -- Load coordinates from SIF data
h = nc.Dataset(path + r"\2018\s5p-l3grd-sif-001-month-20180501-20240325.nc")
lat = h.variables['latitude'][:]
lon = h.variables['longitude'][:]
h.close()

# -- Determine Index bounds
lat_min_index = np.min(np.arange(len(lat))[np.where(lat > lat_min_coord)])
lat_max_index = np.max(np.arange(len(lat))[np.where(lat < lat_max_coord)])
lon_min_index = np.min(np.arange(len(lon))[np.where(lon > lon_min_coord)])
lon_max_index = np.max(np.arange(len(lon))[np.where(lon < lon_max_coord)])

# -- Extract data within bounds
lat_sam = lat[lat_min_index:(lat_max_index + 1)]
lon_sam = lon[lon_min_index:(lon_max_index + 1)]

print("SIF latitude and longitude bounds determined.")

# -- Define time axes
SIF_years = np.arange(ntime_year) + 2018
SIF_months = np.arange(ntime_month) + 1

# -- Pre-allocation for SIF data
SIF = np.zeros([ntime_year, ntime_month, len(lat_sam), len(lon_sam)])  # Adjusted dimensions

# -- Loop over years and months to load SIF data
for ii in np.arange(ntime_year):
    for jj in np.arange(ntime_month):
        if SIF_years[ii] == 2018 and SIF_months[jj] < 5:
            SIF[ii, jj, :, :] = np.nan
            print("if statement correct")  # Checkpoint
        else:
            filename = glob.glob(path + str(SIF_years[ii]) + 
                                 '/s5p-l3grd-sif-001-month-' + 
                                 str(SIF_years[ii]) + 
                                 str(SIF_months[jj]).zfill(2) + 
                                 '01-*')[0]
            h = nc.Dataset(filename)
            SIF[ii, jj, :, :] = h.variables['solar_induced_fluorescence'][
                0, lat_min_index:(lat_max_index + 1), lon_min_index:(lon_max_index + 1)]
            print("else statement correct")  # Checkpoint
            h.close()

# -- Clear loop indices
del ii, jj

print("SIF data loaded.")
#%%
# AMAZON REGIONS MASK

# -- Load Amazon Mask Data using SIF bounds
mask_path = r"C:\Users\Jalin005\OneDrive - Universiteit Utrecht\Documents\01 Python Scripts\Data\Amazon_Ecoregions.nc"  # Path to mask file
mask_data = nc.Dataset(mask_path)

mask_lat = mask_data.variables['lat'][:]
mask_lon = mask_data.variables['lon'][:]

mask_values = mask_data.variables['mask_BiGeoGraLi'][:]
mask_AmazonFlatPlains = mask_data.variables['mask_AmazonFlatPlains'][:]
mask_CerradoCaat = mask_data.variables['mask_CerradoCaat'][:]
mask_AmazonAndesPiedmont = mask_data.variables['mask_AmazonAndesPiedmont'][:]
mask_GuianaShield = mask_data.variables['mask_GuianaShield'][:]
mask_BrazilianShield = mask_data.variables['mask_BrazilianShield'][:]

# Replace 0 values with NaN
mask_values[mask_values == 0] = np.nan
mask_AmazonFlatPlains[mask_AmazonFlatPlains == 0] = np.nan
mask_CerradoCaat[mask_CerradoCaat == 0] = np.nan
mask_AmazonAndesPiedmont[mask_AmazonAndesPiedmont == 0] = np.nan
mask_GuianaShield[mask_GuianaShield == 0] = np.nan
mask_BrazilianShield[mask_BrazilianShield == 0] = np.nan

# Subset mask using SIF index bounds
mask_subset = mask_values[
    np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][0]:np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][-1] + 1,
    np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][0]:np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][-1] + 1
]

mask_AmazonFlatPlains = mask_AmazonFlatPlains[
    np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][0]:np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][-1] + 1,
    np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][0]:np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][-1] + 1
]

mask_CerradoCaat = mask_CerradoCaat[
    np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][0]:np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][-1] + 1,
    np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][0]:np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][-1] + 1
]

mask_AmazonAndesPiedmont = mask_AmazonAndesPiedmont[
    np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][0]:np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][-1] + 1,
    np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][0]:np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][-1] + 1
]

mask_GuianaShield = mask_GuianaShield[
    np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][0]:np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][-1] + 1,
    np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][0]:np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][-1] + 1
]

mask_BrazilianShield = mask_BrazilianShield[
    np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][0]:np.where((mask_lat >= lat_sam.min()) & (mask_lat <= lat_sam.max()))[0][-1] + 1,
    np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][0]:np.where((mask_lon >= lon_sam.min()) & (mask_lon <= lon_sam.max()))[0][-1] + 1
]

print('amazon masks loaded and subset')
#%%
# ATTO EDDY COVARIANCE DATA

# locations EC-fluxtowers
lat_k34 = -2.6 # deg N
lon_k34 = -60.2 # deg E
lat_sa1, lon_sa1 = -2.8567, -54.9589  # SA1
lat_guyaflux, lon_guyaflux = 5.278, -52.92  # GUYAFlux
lat_atto, lon_atto = -2.1483, -59.0003  # Atto Tower

# Data EC-fluxtower
ec_file = "ATTO_monthly_eddyfluxes_deltaCO2_2014_2023.csv"
EC = pd.read_csv(pathEC + ec_file)
EC['Time'] = pd.to_datetime(EC['Time']) # 'Time'-column to datetime

# Filter the EC data for the wished for time frame
start_date = pd.to_datetime('2018-01-01')
end_date = pd.to_datetime('2023-12-31')
EC = EC[(EC['Time'] >= start_date) & (EC['Time'] <= end_date)]
EC.reset_index(drop=True, inplace=True)
EC_GPP_data = EC['GPP_gCper_m2_month']
EC['Month'] = EC['Time'].dt.month  # take the month out of time

# Calculate average per month over the years
EC_monthly_mean = EC.groupby('Month')['GPP_gCper_m2_month'].mean()
EC_monthly_std = EC.groupby('Month')['GPP_gCper_m2_month'].std() # Calculate the standard deviation per month over the years


# Print monthly average EC
print(EC_monthly_mean)


print("EC-data loaded and subset:", EC.shape)

#%%################################################################################################################################################
#                                                                                                                                                 #
# CALCULATION                                                                                                                                     #
#                                                                                                                                                 #
###################################################################################################################################################
"""
# Cut out the Amazon of the SIF data
SIF_southAmerica = SIF
SIF = SIF * mask_subset[np.newaxis, np.newaxis, :, :]  # Broadcast over time dimensions
SIF_AmazonFlatPlains = SIF * mask_AmazonFlatPlains[np.newaxis, np.newaxis, :, :]
SIF_CerradoCaat = SIF * mask_CerradoCaat[np.newaxis, np.newaxis, :, :]
SIF_AmazonAndesPiedmont = SIF * mask_AmazonAndesPiedmont[np.newaxis, np.newaxis, :, :]
SIF_GuianaShield = SIF * mask_GuianaShield[np.newaxis, np.newaxis, :, :]
SIF_BrazilianShield = SIF * mask_BrazilianShield[np.newaxis, np.newaxis, :, :]

# Function to calculate SIF statistics
def calculate_sif_statistics(SIF_input):
    SIF_avgYears = np.nanmean(SIF_input[:,:,:,:], axis=0)
    SIF_avg = np.nanmean(SIF_input[:,:,:,:], axis=(0, 1))
    SIF_2018 = np.nanmean(SIF_input[0,:,:,:], axis=(1,2))
    SIF_2019 = np.nanmean(SIF_input[1,:,:,:], axis=(1,2))
    SIF_2020 = np.nanmean(SIF_input[2,:,:,:], axis=(1,2))
    SIF_2021 = np.nanmean(SIF_input[3,:,:,:], axis=(1,2))
    SIF_2022 = np.nanmean(SIF_input[4,:,:,:], axis=(1,2))
    SIF_2023 = np.nanmean(SIF_input[5,:,:,:], axis=(1,2))
    SIF_comb = np.concatenate((SIF_2018, SIF_2019, SIF_2020, SIF_2021, SIF_2022, SIF_2023)) # Combine the spatially averaged SIF time series
    SIF_meanyear = np.nanmean(SIF_input, axis=(0, 2, 3))
    SIF_stdyear = np.nanstd(SIF_input, axis=(0, 2, 3))
    SIF_monthly_mean = np.nanmean([SIF_2018, SIF_2019, SIF_2020, SIF_2021, SIF_2022, SIF_2023], axis=0)
    SIF_5years_monthly_mean = np.concatenate((SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean))
    SIF_std_monthly = np.nanstd([SIF_2018, SIF_2019, SIF_2020, SIF_2021, SIF_2022, SIF_2023], axis=0)
    return (SIF_avgYears, SIF_avg, SIF_comb, SIF_meanyear, SIF_stdyear, SIF_monthly_mean, SIF_5years_monthly_mean, SIF_std_monthly)

# Compute statistics for each SIF region
(SIF_avgYears, SIF_avg, SIF_comb, SIF_meanyear, SIF_stdyear, SIF_monthly_mean, SIF_5years_monthly_mean, SIF_std_monthly) = calculate_sif_statistics(SIF)
(SIF_AmazonFlatPlains_avgYears, SIF_AmazonFlatPlains_avg, SIF_AmazonFlatPlains_comb, SIF_AmazonFlatPlains_meanyear, SIF_AmazonFlatPlains_stdyear, SIF_AmazonFlatPlains_monthly_mean, SIF_AmazonFlatPlains_5years_monthly_mean, SIF_AmazonFlatPlains_std_monthly) = calculate_sif_statistics(SIF_AmazonFlatPlains)
(SIF_CerradoCaat_avgYears, SIF_CerradoCaat_avg, SIF_CerradoCaat_comb, SIF_CerradoCaat_meanyear, SIF_CerradoCaat_stdyear, SIF_CerradoCaat_monthly_mean, SIF_CerradoCaat_5years_monthly_mean, SIF_CerradoCaat_std_monthly) = calculate_sif_statistics(SIF_CerradoCaat)
(SIF_AmazonAndesPiedmont_avgYears, SIF_AmazonAndesPiedmont_avg, SIF_AmazonAndesPiedmont_comb, SIF_AmazonAndesPiedmont_meanyear, SIF_AmazonAndesPiedmont_stdyear, SIF_AmazonAndesPiedmont_monthly_mean, SIF_AmazonAndesPiedmont_5years_monthly_mean, SIF_AmazonAndesPiedmont_std_monthly) = calculate_sif_statistics(SIF_AmazonAndesPiedmont)
(SIF_GuianaShield_avgYears, SIF_GuianaShield_avg, SIF_GuianaShield_comb, SIF_GuianaShield_meanyear, SIF_GuianaShield_stdyear, SIF_GuianaShield_monthly_mean, SIF_GuianaShield_5years_monthly_mean, SIF_GuianaShield_std_monthly) = calculate_sif_statistics(SIF_GuianaShield)
(SIF_BrazilianShield_avgYears, SIF_BrazilianShield_avg, SIF_BrazilianShield_comb, SIF_BrazilianShield_meanyear, SIF_BrazilianShield_stdyear, SIF_BrazilianShield_monthly_mean, SIF_BrazilianShield_5years_monthly_mean, SIF_BrazilianShield_std_monthly) = calculate_sif_statistics(SIF_BrazilianShield)

"""

#%% old version, that is quicker for the pythoncourse on the 1st of April

# -- Determine climatology

# temporal average
SIF_avgYears = np.nanmean(SIF[:,:,:,:], axis=(0))  # Average over years, ignoring NaNs

SIF_avg = np.nanmean(SIF[:,:,:,:], axis=(0, 1))  # Average over years and months, ignoring NaNs, only spatial information left


# spatially average SIF
SIF_2018 = np.nanmean(SIF[0,:,:,:],axis=(1,2)) #spatially average
SIF_2019 = np.nanmean(SIF[1,:,:,:],axis=(1,2)) #spatially average
SIF_2020 = np.nanmean(SIF[2,:,:,:],axis=(1,2)) #spatially average
SIF_2021 = np.nanmean(SIF[3,:,:,:],axis=(1,2)) #spatially average
SIF_2022 = np.nanmean(SIF[4,:,:,:],axis=(1,2)) #spatially average
SIF_2023 = np.nanmean(SIF[5,:,:,:],axis=(1,2)) #spatially average


# Combine the spatially averaged SIF time series
SIF_comb = np.concatenate((SIF_2018, SIF_2019, SIF_2020, SIF_2021, SIF_2022, SIF_2023))

SIF_meanyear = np.nanmean(SIF, axis=(0, 2, 3))
SIF_stdyear = np.nanstd(SIF, axis=(0, 2, 3))

# Bereken het maandgemiddelde over de jaren (per maand)
SIF_monthly_mean = np.nanmean([SIF_2018, SIF_2019, SIF_2020, SIF_2021, SIF_2022, SIF_2023], axis=0)

SIF_5years_monthly_mean = np.concatenate((SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean, SIF_monthly_mean))

# Bereken de standaarddeviatie over de jaren voor elke maand
SIF_std_monthly = np.nanstd([SIF_2018, SIF_2019, SIF_2020, SIF_2021, SIF_2022, SIF_2023], axis=0)





#%% ATTO calculations

### select SIF for the pixel that overlaps the ATTO tower

# find the closest pixel
lat_index_atto = np.abs(lat_sam - lat_atto).argmin()
lon_index_atto = np.abs(lon_sam - lon_atto).argmin()

# get the timeseries for the Atto-pixel out of the SIF dataset
SIF_atto_subset = SIF[:,:, lat_index_atto, lon_index_atto]
atto_pixel_time_series_comb = np.concatenate(SIF_atto_subset, axis=0)

# spatially average SIF
atto_SIF_2018 = SIF_atto_subset[0,:]  # select 2018
atto_SIF_2019 = SIF_atto_subset[1,:]  # select 2019
atto_SIF_2020 = SIF_atto_subset[2,:]  # select 2020
atto_SIF_2021 = SIF_atto_subset[3,:]  # select 2021
atto_SIF_2022 = SIF_atto_subset[4,:]  # select 2022
atto_SIF_2023 = SIF_atto_subset[5,:]  # select 2023

# Combine the spatially averaged SIF time series
atto_SIF_comb = np.concatenate((atto_SIF_2018, atto_SIF_2019, atto_SIF_2020, atto_SIF_2021, atto_SIF_2022, atto_SIF_2023))

# calculate the monthly average over the years 
atto_SIF_monthly_mean = np.nanmean([atto_SIF_2018, atto_SIF_2019, atto_SIF_2020, atto_SIF_2021, atto_SIF_2022, atto_SIF_2023], axis=0)


# calculate the monthly std over the years
atto_SIF_std_monthly = np.nanstd([atto_SIF_2018, atto_SIF_2019, atto_SIF_2020, atto_SIF_2021, atto_SIF_2022, atto_SIF_2023], axis=0)

### select SIF for the pixels around the ATTO tower

# Define a 3x3 grid around the ATTO pixel
lat_range = (lat_index_atto - 1, lat_index_atto + 2)
lon_range = (lon_index_atto - 1, lon_index_atto + 2)

# Ensure the indices stay within valid bounds of the dataset
lat_range = np.clip(lat_range, 0, len(lat_sam) - 1)
lon_range = np.clip(lon_range, 0, len(lon_sam) - 1)

# Extract SIF data for the surrounding pixels
SIF_around_atto_subset = SIF[:, :, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]


# spatially average SIF
around_atto_SIF_2018 = np.nanmean(SIF_around_atto_subset[0,:,:,:],axis=(1,2)) #spatially average
around_atto_SIF_2019 = np.nanmean(SIF_around_atto_subset[1,:,:,:],axis=(1,2)) #spatially average
around_atto_SIF_2020 = np.nanmean(SIF_around_atto_subset[2,:,:,:],axis=(1,2)) #spatially average
around_atto_SIF_2021 = np.nanmean(SIF_around_atto_subset[3,:,:,:],axis=(1,2)) #spatially average
around_atto_SIF_2022 = np.nanmean(SIF_around_atto_subset[4,:,:,:],axis=(1,2)) #spatially average
around_atto_SIF_2023 = np.nanmean(SIF_around_atto_subset[5,:,:,:],axis=(1,2)) #spatially average


# Combine the spatially averaged SIF time series
around_atto_SIF_comb = np.concatenate((around_atto_SIF_2018,around_atto_SIF_2019, around_atto_SIF_2020, around_atto_SIF_2021, around_atto_SIF_2022, around_atto_SIF_2023))

# Bereken het maandgemiddelde over de jaren (per maand)
around_atto_SIF_monthly_mean = np.nanmean([around_atto_SIF_2018, around_atto_SIF_2019, around_atto_SIF_2020, around_atto_SIF_2021, around_atto_SIF_2022, around_atto_SIF_2023], axis=0)


# Bereken de standaarddeviatie over de jaren voor elke maand
around_atto_SIF_std_monthly = np.nanstd([around_atto_SIF_2018, around_atto_SIF_2019, around_atto_SIF_2020, around_atto_SIF_2021, around_atto_SIF_2022, around_atto_SIF_2023], axis=0)


#%%################################################################################################################################################
#                                                                                                                                                 #
# RESULTS                                                                                                                                         #
#                                                                                                                                                 #
###################################################################################################################################################

#%%############################################################################
# PLOT 1) SPATIAL MAP OF AVERAGE SIF INCLUDING EC-TOWERS
###############################################################################

plt.rcParams['mathtext.default'] = 'regular'
plt.figure(figsize=(10,5),facecolor='w')
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-90,-30,-30,20],crs=ccrs.PlateCarree())
ax.coastlines(color='k')
ax.add_feature(cfeature.BORDERS,linestyle='-')
x = lon_sam
y = lat_sam
xs, ys = np.meshgrid(x,y)
zs = SIF_avg

mesh = ax.pcolormesh(xs,ys,zs,transform=ccrs.PlateCarree(),vmin=0.7,vmax=1.5)       # change y-axis range!
plt.plot(lon_k34, lat_k34,'ro',label='K34')
plt.plot(lon_sa1, lat_sa1, 'bo', label='SA1')  
plt.plot(lon_guyaflux, lat_guyaflux, 'mo', label='GUYAFlux')  

# Label to indicate EC towers
plt.text(lon_k34 + 1, lat_k34, 'K34', color='black', fontsize=10, transform=ccrs.PlateCarree())
plt.text(lon_sa1 + 1, lat_sa1, 'SA1', color='black', fontsize=10, transform=ccrs.PlateCarree())
plt.text(lon_guyaflux + 1, lat_guyaflux, 'GF-Guy', color='black', fontsize=10, transform=ccrs.PlateCarree())
# Voeg ATTO-toren toe aan de plot
plt.plot(lon_atto, lat_atto, 'go', label='ATTO')  # groen punt voor ATTO
plt.text(lon_atto + 1, lat_atto, 'ATTO', color='black', fontsize=10, transform=ccrs.PlateCarree())


gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                  linewidth=1,color='grey',alpha=0.5,linestyle='--')
gl.top_labels = False
gl.right_labels = False
plt.colorbar(mesh,label='SIF (mW/m$^2$/sr/nm)',extend='max')
plt.title('SIF_avg (2018-2023)')
if save == True:
    plt.savefig('SIF_avg.png',dpi=300,bbox_inches='tight')
else:
    plt.show()
  
#%%############################################################################
# PLOT 2) SEASONALITY OF SIF FOR DIFFERENT YEARS
###############################################################################

# Define months for the x-axis
months = np.arange(1, 13)  # From January (1) to December (12)

# Plot the data for each year, plotting SIF for each month
plt.figure()
plt.plot(months, SIF_2018, label='2018')
plt.plot(months, SIF_2019, label='2019')
plt.plot(months, SIF_2020, label='2020')
plt.plot(months, SIF_2021, label='2021')
plt.plot(months, SIF_2022, label='2022')
plt.plot(months, SIF_2023, label='2023')

# Labeling the axes and title
plt.xlabel('Month')
plt.ylabel('SIF (mW/m$^2$/sr/nm)')
plt.title('Monthly SIF for 2018-2023')
plt.legend()

# Adjusting the x-axis to show months properly
plt.xticks(months, ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

plt.show()

#%%############################################################################
# PLOT 3) SEASONALITY OF SIF FOR DIFFERENT YEARS, WITH AVERAGE AND Std SHADING
###############################################################################

months = np.arange(1, 13)
plt.figure(figsize=(10, 6))

# Plot the monthly mean
plt.plot(months, SIF_monthly_mean, color='black', linestyle='--', label='Monthly Average')
plt.fill_between(months, SIF_monthly_mean - SIF_std_monthly, SIF_monthly_mean + SIF_std_monthly, color='gray', alpha=0.2, label='Standard Deviation') # shading

# Plot the data for each year, plotting SIF for each month
plt.plot(months, SIF_2018, label='2018')
plt.plot(months, SIF_2019, label='2019')
plt.plot(months, SIF_2020, label='2020')
plt.plot(months, SIF_2021, label='2021')
plt.plot(months, SIF_2022, label='2022')
plt.plot(months, SIF_2023, label='2023')

# Labeling the axes and title
plt.xlabel('Month')
plt.ylabel('SIF (mW/m$^2$/sr/nm)')
plt.title('Monthly SIF for 2018-2023 with Std Dev Shading')
plt.legend()

# Adjusting the x-axis to show months properly
plt.xticks(months, ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

plt.show()

#%%############################################################################
# PLOT 4) AVERAGE SEASONALITY OF SIF WITH 1 STD DEV
###############################################################################
# Plotting
plt.figure(figsize=(10,6))
plt.plot(SIF_meanyear, label='Mean SIF')
plt.fill_between(np.arange(12), SIF_meanyear + SIF_stdyear, SIF_meanyear - SIF_stdyear, alpha=0.5, label='±1 Std Dev')
plt.ylim(0, 2)
plt.ylabel('SIF (unit)')
plt.xticks(np.arange(12), ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
plt.legend()
plt.title('Yearly SIF with Standard Deviation')
plt.show()


#%%############################################################################
# PLOT 5) TIMESERIES OF SPATIALLY AVERAGED SIF
###############################################################################
# Assuming SIF_comb is a flattened array with SIF data for each month from Jan 2018 to Dec 2023 (72 months in total)
plt.figure(figsize=(12, 6))  # Wider figure for better visibility

# Plot the combined SIF data
plt.plot(SIF_comb)

# Set the x-axis limits to reflect the range of months from Jan 2018 to Dec 2023
plt.xlim(0, 71)  # 72 months in total, from 0 to 71

# Set the x-ticks and labels for each half year from Jan 2018 to Dec 2023
plt.xticks(
    [0, 5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71],  # Every 6 months
    ['Jan 2018', 'Jun 2018', 'Dec 2018', 'Jun 2019', 'Dec 2019', 'Jun 2020', 'Dec 2020', 
     'Jun 2021', 'Dec 2021', 'Jun 2022', 'Dec 2022', 'Jun 2023', 'Dec 2023']
)

# Labeling the axes
plt.xlabel('Month')
plt.ylabel('SIF (mW/m$^2$/sr/nm)')
plt.title('Monthly Combined SIF from Jan 2018 to Dec 2023')

# Adjust the margins to add whitespace on the right
plt.subplots_adjust(right=0.9)  # Default is 1.0; reducing it adds space

# Show the plot
plt.show()


#%%############################################################################
# PLOT 6) TIMESERIES OF SPATIALLY AVERAGED SIF, WITH YEARLY AVERAGE
###############################################################################
import matplotlib.pyplot as plt
import numpy as np

# Assuming SIF_comb and SIF_5years_monthly_mean are already defined
plt.figure(figsize=(12, 6))  # Wider figure for better visibility

# Plot the combined SIF data
plt.plot(SIF_comb, label="Combined SIF Data", color='blue')

# Plot the 5-year monthly mean as a light blue dashed line
plt.plot(SIF_5years_monthly_mean, label="5-Year Monthly Mean", linestyle='--', color='cornflowerblue')

# Set the x-axis limits to reflect the range of months from Jan 2018 to Dec 2023
plt.xlim(0, 71)  # 72 months in total, from 0 to 71

# Set the x-ticks and labels for each half year from Jan 2018 to Dec 2023
plt.xticks(
    [0, 5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71],  # Every 6 months
    ['Jan 2018', 'Jun 2018', 'Dec 2018', 'Jun 2019', 'Dec 2019', 'Jun 2020', 'Dec 2020', 
     'Jun 2021', 'Dec 2021', 'Jun 2022', 'Dec 2022', 'Jun 2023', 'Dec 2023']
)

# Labeling the axes
plt.xlabel('Month')
plt.ylabel('SIF (mW/m$^2$/sr/nm)')
plt.title('Monthly Combined SIF from Jan 2018 to Dec 2023')

# Add a legend to differentiate the two lines
plt.legend()

# Adjust the margins to add whitespace on the right
plt.subplots_adjust(right=0.9)  # Default is 1.0; reducing it adds space

# Show the plot
plt.show()

#%%##########################################################################################
# PLOT 7) TIMESERIES OF SPATIALLY AVERAGED SIF, AND ATTO EC-FLUX DATA
#############################################################################################

# Aanpassen van de plot met eenvoudiger layout
plt.figure(figsize=(12, 6))  # Breder figuur voor betere zichtbaarheid

# Eerste as voor SIF
fig, ax1 = plt.subplots(figsize=(12, 6))

# Plot alleen de SIF data (geen 5-year average)
ax1.plot(SIF_comb, label="SIF Amazon", color='blue')

# Stel de x-as in
ax1.set_xlim(0, 71)  # 72 maanden, van 0 tot 71
ax1.set_xticks([0, 5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71])  # Elke 6 maanden
ax1.set_xticklabels(
    ['Jan 2018', 'Jun 2018', 'Dec 2018', 'Jun 2019', 'Dec 2019', 'Jun 2020', 'Dec 2020',
     'Jun 2021', 'Dec 2021', 'Jun 2022', 'Dec 2022', 'Jun 2023', 'Dec 2023']
)
ax1.set_xlabel('Month')
ax1.set_ylabel('SIF (mW/m$^2$/sr/nm)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')

# Tweede as voor EC GPP-data
ax2 = ax1.twinx()
ax2.plot(EC_GPP_data, label="GPP Data", color='green', linestyle='--')  # Stippellijn voor GPP
ax2.set_ylabel('GPP (gC m$^{-2}$ month$^{-1}$)', color='green')
ax2.tick_params(axis='y', labelcolor='green')

# Titel en legenda
fig.suptitle('Monthly SIF and GPP Data (Jan 2018 - Dec 2023)', fontsize=14)

# Combineer de legende
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
fig.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right", bbox_to_anchor=(0.9, 0.88))

# Toon de plot
plt.show()

#%%##########################################################################################
# PLOT 8) TIMESERIES OF ATTO ESA SIF, AND ATTO EC-FLUX DATA
#############################################################################################
#for just the atto pixel, use atto_pixel_time_series_comb 
# for the 9 pixels around the atto pixel, use around_atto_SIF_comb

# Aanpassen van de plot met eenvoudiger layout
plt.figure(figsize=(12, 6))  # Breder figuur voor betere zichtbaarheid

# Eerste as voor SIF
fig, ax1 = plt.subplots(figsize=(12, 6))

# Plot alleen de SIF data (geen 5-year average)
ax1.plot(atto_pixel_time_series_comb, label="SIF Data ATTO pixel", color='blue')

# Stel de x-as in
ax1.set_xlim(0, 71)  # 72 maanden, van 0 tot 71
ax1.set_xticks([0, 5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 71])  # Elke 6 maanden
ax1.set_xticklabels(
    ['Jan 2018', 'Jun 2018', 'Dec 2018', 'Jun 2019', 'Dec 2019', 'Jun 2020', 'Dec 2020',
     'Jun 2021', 'Dec 2021', 'Jun 2022', 'Dec 2022', 'Jun 2023', 'Dec 2023']
)
ax1.set_xlabel('Month')
ax1.set_ylabel('SIF (mW/m$^2$/sr/nm)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')

# Tweede as voor EC GPP-data
ax2 = ax1.twinx()
ax2.plot(EC_GPP_data, label="GPP Data", color='green', linestyle='--')  # Stippellijn voor GPP
ax2.set_ylabel('GPP (gC m$^{-2}$ month$^{-1}$)', color='green')
ax2.tick_params(axis='y', labelcolor='green')

# Titel en legenda
fig.suptitle('Monthly SIF and GPP Data (Jan 2018 - Dec 2023)', fontsize=14)

# Combineer de legende
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
fig.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right", bbox_to_anchor=(0.9, 0.88))

# Toon de plot
plt.show()

#%%##########################################################################################
# PLOT 8b) TIMESERIES OF ATTO ESA SIF, AND ATTO EC-FLUX DATA DIFFERENTLY
#############################################################################################
# Voor de plot gebruiken we de maandgemiddelden van SIF (rond ATTO en bij ATTO) en de EC-gegevens
months = np.arange(1, 13)

# Maak een figuur en eerste as voor SIF-gegevens
fig, ax1 = plt.subplots(figsize=(12, 6))

# Plot de maandgemiddelden en standaarddeviaties voor SIF rond ATTO
ax1.plot(months, around_atto_SIF_monthly_mean, label='Around ATTO Monthly Mean', color='blue', linestyle='--')
ax1.fill_between(months, around_atto_SIF_monthly_mean - around_atto_SIF_std_monthly, around_atto_SIF_monthly_mean + around_atto_SIF_std_monthly, color='blue', alpha=0.2, label='Around ATTO Std Dev')

# Plot de maandgemiddelden en standaarddeviaties voor SIF bij ATTO
ax1.plot(months, atto_SIF_monthly_mean, label='ATTO Monthly Mean', color='red', linestyle='--')
ax1.fill_between(months, atto_SIF_monthly_mean - atto_SIF_std_monthly, atto_SIF_monthly_mean + atto_SIF_std_monthly, color='red', alpha=0.2, label='ATTO Std Dev')

# Stel de x-as en labels voor de eerste y-as (SIF)
ax1.set_xlim(1, 12)
ax1.set_xticks(months)
ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax1.set_xlabel('Month')
ax1.set_ylabel('SIF (mW/m$^2$/sr/nm)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')

# Maak een tweede y-as voor de EC GPP-gegevens
ax2 = ax1.twinx()
ax2.plot(months, EC_monthly_mean, label='EC Monthly Mean', color='green', linestyle='-')
ax2.fill_between(months, EC_monthly_mean - EC_monthly_std, EC_monthly_mean + EC_monthly_std, color='green', alpha=0.2, label='EC Std Dev')

# Stel de labels en y-as voor de tweede as (EC)
ax2.set_ylabel('GPP (gC m$^{-2}$ month$^{-1}$)', color='green')
ax2.tick_params(axis='y', labelcolor='green')

# Titel en gecombineerde legenda
fig.suptitle('Monthly SIF and EC (GPP) Data (2018-2023)', fontsize=14)

# Combineer de legendes van beide assen
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
fig.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right", bbox_to_anchor=(0.6, 0.88))

# Toon de plot
plt.show()


