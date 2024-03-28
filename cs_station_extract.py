'''
Extracting Stations for Canada OA White paper.

Frederic.Cyr@dfo-mpo.gc.ca - Sept. 2023
'''

import os
import pandas as pd
#matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import from_levels_and_colors
import numpy as np
from math import radians, cos, sin, asin, sqrt
import azmp_sections_tools as azst
import cmocean.cm as cmo


## ----- Parameters to edit ---- ##
REGION = 'NL'
SEASON = 'fall'
YEAR = 2014# plotting parameters
# File to load
my_file = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv'


# Get the data
df = pd.read_csv(my_file)
# set index
df.index = pd.to_datetime(df.Timestamp)
df.drop(columns='Timestamp', inplace=True)
# Extract Region
#df = df[df.Region==REGION]

# Extract year
#df = df[df.index.year == YEAR]

# Extract season
if SEASON == 'spring':
    df = df[df.index.month <= 5]
elif SEASON == 'summer':
    df = df[(df.index.month>=6) & (df.index.month<=8)]
elif SEASON == 'fall':
    df = df[df.index.month >= 9]


# keep some columns
df = df[['Region', 'Station_Name', 'Latitude_(degNorth)', 'Longitude_(degEast)']]

# drop duplicates
df.drop_duplicates('Station_Name', inplace=True)

# save file
df.to_csv('Regular_stations_all.csv')
