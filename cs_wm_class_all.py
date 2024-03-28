# -*- coding: utf-8 -*-
"""
Script to plot T-S diagram of water mass class derived by B. Zemskova

Frederic.Cyr@dfo-mpo.gc.ca
July 2021
Update Fall 2023 with new data

"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import water_masses as wm
import gsw

# Adjust fontsize/weight
font = {'weight' : 'bold',
        'size'   : 15}
plt.rc('font', **font)


## Input parameters
TSDIC = False
class_type = '6_classes'


var_y = 'Temperature_(degC)'
var_x = 'Salinity_(psu)'
var_col = 'Depth_(dbar)' #kml or region_number or Region
#var_col = 'Omega_Aragonite_(--)' #kml or region_number or Region
var_col_sort = False
var_sym = 'Region' #kml or Region
if TSDIC:
    fig_name = 'azmp_T-S-DIC_carbon_wm_' + class_type + '.png'
else:
    fig_name = 'azmp_T-S_carbon_wm_' + class_type + '.png'
# Default colors cycle
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

## Load data and some checks
# 2023 version:
if TSDIC:
    df = pd.read_csv('./Testing_classification/all_data_CT_SA_DIC_classes.csv', delimiter=',')
else:
    df = pd.read_csv('./Testing_classification/all_data_CT_SA_classes.csv', delimiter=',')
# 2021 version:
#df = pd.read_csv('./Testing_classification/AZMP_df_class_TSDIC.csv', delimiter=',')
#df = pd.read_csv('./Testing_classification/AZMP_df_class_TSonly.csv', delimiter=',')
df.set_index('Timestamp', inplace=True)
df.index = pd.to_datetime(df.index)

## Derive variables using TEOS-10
SP = df['Salinity_(psu)']
T = df['Temperature_(degC)']
p = df['Depth_(dbar)']
lat = df['Latitude_(degNorth)']
lon = df['Longitude_(degEast)']
SA = gsw.SA_from_SP(SP,p,lon, lat)
CT = gsw.CT_from_t(SA,T,p)
SIG0 = gsw.sigma0(SA, CT)

fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_figheight(6)
fig.set_figwidth(8)


wm_class = np.sort(df[class_type].unique())
wm_class_text = wm_class+1
co = np.full(len(wm_class),'nan')
for i, iclass in enumerate(wm_class):
    SA_tmp = SA[df[class_type]==iclass]
    CT_tmp = CT[df[class_type]==iclass]
    s = plt.scatter(SA_tmp, CT_tmp, s=15, color=colors[i], lw=0.2, alpha=0.9)
    #s = plt.plot(SA_tmp, CT_tmp, '.', alpha=0.9)
    #co[i] = s[0].get_color()

plt.legend(wm_class_text.astype(str))
plt.ylabel(r'$\Theta_T$', fontsize=15, fontweight='normal')
plt.xlabel(r'$S_A$', fontsize=15, fontweight='normal')
ax.tick_params(axis='y', direction='out')
ax.tick_params(axis='x', direction='out')


## add isopycnals    
# Figure out boudaries (mins and maxs)
smin = np.nanmin(np.ceil(SA)) - (0.01 * np.nanmin(np.ceil(SA)))
smax = np.nanmax(np.floor(SA)) + (0.01 * np.nanmax(np.floor(SA)))
tmin = np.nanmin(np.ceil(CT)) - (0.01 * np.nanmin(np.ceil(CT)))
tmax = np.nanmax(np.floor(CT)) + (0.01 * np.nanmax(np.floor(CT)))
smin = np.nanmin(np.ceil(SA))
smax = np.nanmax(np.floor(SA))
tmin = np.nanmin(np.ceil(CT))
tmax = np.nanmax(np.floor(CT))
# Calculate how many gridcells we need in the x and y dimensions
xdim = int(round((smax-smin)/0.05+1,0))
ydim = int(round((tmax-tmin)/0.05+1,0))
# Create empty grid of zeros
sig0 = np.full((ydim,xdim), 0)
# Create temp and salt vectors of appropiate dimensions
ti = np.linspace(1,ydim-1,ydim)*0.05+tmin
si = np.linspace(1,xdim-1,xdim)*0.05+smin
# Loop to fill in grid with densities
for j in range(0,int(ydim)):
    for i in range(0, int(xdim)):
        sig0[j,i]=gsw.sigma0(si[i],ti[j])
# add to plot
CS = plt.contour(si,ti,sig0, linestyles=':', colors='silver')
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level

## Add freezing seawater line 
Tf = CT_freezing = gsw.CT_freezing(si,0,0)
plt.plot(si, Tf, '--k')

## Add water masses definitions
wm_def = wm.water_masses_def_petrie()
for i, w in enumerate(wm_def.keys()):
    shape = np.array(wm_def.get(w))
    p = plt.plot(shape[:,0], shape[:,1], color='k' )
    co = p[0].get_color() 
    plt.text(np.mean(shape[:,0]), np.mean(shape[:,1]), w, color=co)

## Add Fratantoni def
mS = [34.3, 35, 35, 34.3, 34.3]
mT = [4, 4, 8, 8, 4]
p = plt.plot(mS, mT, color='k', linewidth=2 )
plt.text(np.mean(mS), np.mean(mT), 'LSlW', color='k')
mS = [34.7, 35.5, 35.5, 34.7, 34.7]
mT = [8, 8, 12, 12, 8]
p = plt.plot(mS, mT, color='k', linewidth=2 )
plt.text(np.mean(mS), np.mean(mT), 'WSW', color='k')

## Add Jutras Def.
## # CIL
## mS = [31.3, 33.5, 33.5, 31.3, 31.3]
## mT = [-1.7, -1.7, 2.6, 2.6, -1.7]
## p = plt.plot(mS, mT, color='m', linewidth=2 )
## plt.text(np.mean(mS), np.mean(mT), 'CIL', color='m')
## # NACW
## mS = [33.8, 36.2, 36.4, 37.7, 37.7, 35.3, 34, 33.8, 33.8]
## mT = [4.2, 4.2, 7.8, 17.5, 17.9, 17.9, 8.2, 4.6, 4.2] 
## p = plt.plot(mS, mT, color='m', linewidth=2 )
## plt.text(np.mean(mS), np.mean(mT), 'NACW', color='m')
## # LCW
## mS = [32.9, 35.5, 35.5, 32.9, 32.9]
## mT = [-.9, -.9, 3.4, 3.4, -.9]
## p = plt.plot(mS, mT, color='m', linewidth=2 )
## plt.text(np.mean(mS), np.mean(mT), 'LCW', color='m')
   
#### ---- Save Figure ---- ####
fig.set_size_inches(w=12, h=10)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#plt.show()

