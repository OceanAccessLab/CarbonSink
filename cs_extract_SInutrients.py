'''
Data for Tina Schaefer
'''

# -*- coding: utf-8 -*-
"""

To extract bottom data from btl master file for Nina Schaefer

    
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import numpy as np
import h5py
import cmocean
import cmocean.cm as cmo
import cartopy. crs as ccrs
import cartopy.feature as cpf
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cc_variable_list2 as vl

# Read the entire AZMP dataset
df = pd.read_excel('/home/cyrf0006/github/AZMP-NL/datasets/water_samples/AZMP_Bottle_Data.xlsx', header=1)
df['GMT'] = df['GMT'].astype(str)
# correct some errors...
df['GMT'] = df['GMT'].replace('1902-10-26 00:00:00','00:00:00')
df['GMT'] = df['GMT'].replace('1900-01-13 02:38:24','00:00:00')
df['GMT'] = df['GMT'].replace('1900-01-21 13:55:12','00:00:00')
df['GMT'] = df['GMT'].replace('1900-01-16 04:19:12','00:00:00')
for idx, i in enumerate(df['GMT']):
    tmp = i.split(':')
    if len(tmp)==2:
         df['GMT'][idx] = tmp[0] + ':' +  tmp[1] + ':00'

df['timestamp'] = pd.to_datetime(df.Date.astype(str)) + pd.to_timedelta(df.GMT.astype(str)) 

# Find coastal SI
df = df[df.Station.str.contains('SI')]
df = df.drop(df[df.Station == 'SI-12'].index)
df = df.drop(df[df.Station == 'SI-13'].index)
df = df.drop(df[df.Station == 'SI-14'].index)


# Set depth as float (had some problem with some data set)
#df = df.astype({'Nominal_Depth':'float'})  
# find average 50-150m
df = df.loc[(df['Nominal_Depth'] >=50) & (df['Nominal_Depth'] <=150)]
#df = df.loc[df.groupby('Station')['Nominal_Depth'].mean()] #group by station then pull "min or max depth"

# Set index
df.set_index('timestamp', inplace=True)

# monhtly mean
df_monthly = df.resample('M').mean()    

df_monthly.Nitrate.dropna().plot(marker='.', linestyle=' ')    
df_monthly.Silicate.dropna().plot(marker='.', linestyle=' ')
df_monthly.Nitrate.dropna().rolling(6).mean().plot(color='steelblue', linewidth=3)
df_monthly.Silicate.dropna().rolling(6).mean().plot(color='orange', linewidth=3) 
plt.legend(['NO3', 'SiO', 'NO3 smooth', 'SiO smooth'])
plt.grid()
plt.ylabel(r'Concentration ($mmol\,m^{-3}$)')
plt.xlabel(r' ')


R = df_monthly.Nitrate.dropna() / df_monthly.Phosphate.dropna()
R.plot(marker='.', linestyle=' ')
R.rolling(5).mean().plot(color='steelblue', linewidth=2) 
plt.grid()
plt.ylabel(r'NO3 / PO4')
plt.xlabel(r' ')


df_monthly.Nitrate.dropna().plot(marker='.', linestyle=' ')    
df_monthly.Silicate.dropna().plot(marker='.', linestyle=' ')
df_monthly.Nitrate.dropna().rolling(5).mean().plot(color='steelblue', linewidth=2)
df_monthly.Silicate.dropna().rolling(5).mean().plot(color='orange', linewidth=2) 
plt.legend(['NO3', 'SiO', 'NO3 5yr', 'SiO 5yr'])
plt.grid()
plt.ylabel(r'Concentration ($mmol\,m^{-3}$)')
plt.xlabel(r' ')

df_monthly.to_pickle('monthly_nutrients.pkl')

# Drop some columns
#df = df.drop(columns=['Fluorescence_uncalibrated', 'CDOM','Beam_Attenuation_Coefficient', 'PI'])

# save
df.to_csv('carbon_bottom_schaefer.csv', float_format='%.4f', index=True)


