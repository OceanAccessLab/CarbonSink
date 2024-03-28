'''
WORK IN PROGRESS

script to plot annual OA sections

Frederic.Cyr@dfo-mpo.gc.ca - March 2021
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
#import cc_variable_list_NL as vl
# For NL work:
import cc_variable_list_NL as vl
import gsw

## ----- Parameters to edit ---- ##
# Variable you want to plot
VAR = 'Omega_Aragonite_(unitless)'
#VAR = 'pH_Total_(total_scale)'
#VAR = 'Oxygen_Saturation_(%)'
#VAR = 'Dissolved_Oxygen_(mL/L)'
VAR = 'Nitrate_Concentration_(mmol/m3)'
#VAR = 'Silicate_Concentration_(mmol/m3)'
#VAR = 'Phosphate_Concentration_(mmol/m3)'
#VAR = 'Salinity_(psu)'
#VAR = 'Total_Alkalinity_(umol/kg)'
#VAR = 'Inorganic_Carbon_(umol/kg)'
#VAR = 'AOU'
REGION = 'NL'

# File to load
my_file = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data_with2021.csv'
# For colorbar:
vv = vl.variable_parameters(VAR)
num_levels = vv[0]
vmin = vv[1]
vmax = vv[2]
midpoint = vv[3]
colors = vv[4]
ticks = vv[5]
axis_label = vv[6]
extent = vv[7]
v_anom = np.linspace(vv[8], vv[9], vv[10])

# Text for figure:
VAR_text = VAR.split('(')[0][0:-1] 

# adjust the colorbar
levels = np.linspace(vmin, vmax, num_levels)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
colors = colors(vals)
colors=np.concatenate([[colors[0,:]], colors, [colors[-1,:]]],0)
cmap, norm = from_levels_and_colors(levels, colors, extend=extent)

# Get the data
df = pd.read_csv(my_file)
# set index
df.index = pd.to_datetime(df.Timestamp)
df.drop(columns='Timestamp', inplace=True)
# Extract NL data
df = df[df.Region==REGION]

# Calculate Apparent Oxygen Utilisation (AOU) if needed
if VAR == 'AOU':
    SA = gsw.SA_from_SP(df['Salinity_(psu)'], df['Depth_(dbar)'], df['Longitude_(degEast)'], df['Latitude_(degNorth)'])
    CT = gsw.CT_from_t(SA, df['Temperature_(degC)'], df['Depth_(dbar)'])
    O2sol = gsw.O2sol(SA, CT, df['Depth_(dbar)'],  df['Longitude_(degEast)'], df['Latitude_(degNorth)']) # in umol/kg
    O2sol = O2sol/43.570 # in ml/l 
    df['AOU'] = O2sol - df['Dissolved_Oxygen_(mL/L)']

        
## ---- Compute density ---- ##
Z = df['Depth_(dbar)'].values
PT = df['Temperature_(degC)'].values
SP = df['Salinity_(psu)'].values
SA = gsw.SA_from_SP(SP, Z, -50, 47)
CT = gsw.CT_from_pt(SA, PT)
SIG0 = gsw.sigma0(SA, CT)
df['sigma_0'] = SIG0

# Chose shelf only
df = df.loc[df['sigma_0'] <= 27.5]

# Group by time and depth average
df = df.groupby('Timestamp').mean()
df = df.resample('D').mean()

fig = plt.figure()
plt.scatter(df.index, df['Latitude_(degNorth)'], c=df[VAR], s=30)
plt.plot([df.index.min(), df.index.max()],[47, 47], '--k')
plt.colorbar()
plt.ylabel('Latitude (degN)')
plt.title(VAR)
fig.set_size_inches(w=6.5,h=4)
fig_name1 = 'hovmoller_' + VAR_text + '.png' 
fig.savefig(fig_name1, dpi=200)
os.system('convert -trim ' + fig_name1 + ' ' + fig_name1)


# plot integrated (keep only north of 47N)
df = df[df['Latitude_(degNorth)']>=47]
fig = plt.figure()
df[VAR].plot()
plt.grid()
plt.ylabel(VAR)
plt.xlabel(' ')
fig.set_size_inches(w=5.15,h=4)
fig_name2 = 'integrated_hovmoller_' + VAR_text + '.png' 
fig.savefig(fig_name2, dpi=200)
os.system('convert -trim ' + fig_name2 + ' ' + fig_name2)

os.system('montage ' +  fig_name1 + ' ' + fig_name2 + ' -tile 1x2 -geometry +1+1  -background white ' + fig_name1)
os.system('rm ' +  fig_name2)

