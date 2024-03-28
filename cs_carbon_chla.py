'''
Data for B. Zemskova.
Consist of adding chl-a data to NL carbon.

Frederic.Cyr@dfo-mpo.gc.ca
'''

# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import numpy as np
import h5py
import cmocean
import cmocean.cm as cmo
import cc_variable_list2 as vl
import gsw
from PyCO2SYS import CO2SYS 
from PyCO2SYS.meta import version

## ---- Read bottle Masterfile ---- ##
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

# Select time period
df.index = pd.to_datetime(df.timestamp)
df = df[df.index.year>=2014]
# Drop wrong IDs and set to int.
df = df.dropna(subset=['ID'])
df.ID = df.ID.astype('int64')
# Keep only chl-a
df = df[['ID', 'Chla', 'Fluorescence_calibrated']]
df = df.dropna(subset=['Fluorescence_calibrated', 'Chla'], how='all')
# remove duplicates in sample ID
df = df.drop_duplicates(subset=['ID'])
# rename ID to match carbon data
df = df.rename(columns={'ID' : 'Sample_ID'})


## ---- Read AZMP carbon data ---- ##
df_NL = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP-NL_carbon_data_ID.csv')
df_NL.set_index('Timestamp', inplace=True)


## ---- Merge data Chla into carbonate data ---- ##
index = df_NL.index
df.reset_index(inplace=True)
df_NL = df_NL.merge(df, on='Sample_ID', how='left')
df_NL.index = index
del df


## ---- some cleaning and saving ---- ##
df_NL.drop(columns=['timestamp', 'Sample_ID'], inplace=True)
df_NL = df_NL.rename(columns={'Chla' : 'Chla_concentration_(mg/m3)'})
df_NL = df_NL.rename(columns={'Fluorescence_calibrated' : 'Fluorescence_(mg/m3)'})

# Sort the dataset
df_NL.sort_index(inplace=True)        

# Save final dataset
df_NL.to_csv('AZMP-NL_carbon_chla.csv', float_format='%.4f', index=True)

