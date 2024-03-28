'''
Data for Tina Schaefer
'''


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

DEPTH='range'

DEPTH='bottom'


# Read the entire AZMP dataset
df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv', delimiter=',')

# to datetime
df['Timestamp'] = pd.to_datetime(df['Timestamp'])

# Only NL data
df = df[df.Region=='NL']


# Round to the neareast hour (or 2)
df['Timestamp'] = df['Timestamp'].dt.floor('2h')

# Set depth as float (had some problem with some data set)
df = df.astype({'Depth_(dbar)':'float'})  

if DEPTH == 'surface':
    df = df.loc[df.groupby('Timestamp')['Depth_(dbar)'].idxmin()] #group by station then pull "min or max depth"
    df = df.loc[df['Depth_(dbar)'] <20] #take all depths >10m
if DEPTH == 'bottom':
    df = df.loc[df.groupby('Timestamp')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"
    df = df.loc[df['Depth_(dbar)'] >50] #take all depths >10m (for bottom) to eliminate lone surface samples
if DEPTH == 'range':
    df = df.loc[(df['Timestamp'] >=50) & (df.depth <=150)]
    df = df.loc[df.groupby('Station_Name')['Depth_(dbar)'].idxmax()] #group by station then pull "min or max depth"

# Set index
df.set_index('Timestamp', inplace=True)

# Drop some columns
#df = df.drop(columns=['Fluorescence_uncalibrated', 'CDOM','Beam_Attenuation_Coefficient', 'PI'])

# save
df.to_csv('carbon_bottom_schaefer.csv', float_format='%.4f', index=True)


