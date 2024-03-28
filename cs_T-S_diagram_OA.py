# -*- coding: utf-8 -*-
"""
Script to plot T-S diagram of Ocean Carbon parameters.
(ESSD OA paper)

Script originally written by O. Gibb (OA_plots.py)

Modified by: Frederic.Cyr@dfo-mpo.gc.ca
April 2021
"""

from matplotlib.colors import from_levels_and_colors
import gsw
import cc_variable_list2 as vl
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import water_masses as wm

## Input parameters
var_y = 'Temperature_(degC)'
var_x = 'Salinity_(psu)'
var_col = 'Depth_(dbar)' #kml or region_number or Region
#var_col = 'Omega_Aragonite_(--)' #kml or region_number or Region
var_col_sort = False
var_sym = 'Region' #kml or Region
water_masses = True

## Load data and some checks
df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv', delimiter=',')
df.set_index('Timestamp', inplace=True)
df.index = pd.to_datetime(df.index)
df = df[df.index.year<2020] # limit to 2014-2019

## Derive variables using TEOS-10
SP = df['Salinity_(psu)']
T = df['Temperature_(degC)']
p = df['Depth_(dbar)']
lat = df['Latitude_(degNorth)']
lon = df['Longitude_(degEast)']
SA = gsw.SA_from_SP(SP,p,lon, lat)
CT = gsw.CT_from_t(SA,T,p)

## Select Variables
plt_var = [var_y, var_x, var_col]
axis = ['y','x','z']

conditions = [
     df['Region'] == 'GSL',
     df['Region'] == 'NL',
     df['Region'] == 'SS']
choices = ['s','^','o']
df['region_marker'] = np.select(conditions, choices, default=0)
markers = {'SS': 's','NL' : 'o', 'GSL' : '^'}
mark_style = markers
M=['SS', 'NL', 'GSL']

## Plot display
for q, p in zip(plt_var, axis):
    
    v = vl.variable_parameters(q)
    num_levels = v[0]
    vmin = v[1]
    vmax = v[2]
    midpoint = v[3]
    colors = v[4]
    ticks = v[5]
    axis_label = v[6]
    extent = v[7] 
        
    locals()['vmin_'+(p)] = vmin
    locals()['vmax_'+(p)] = vmax
    locals()['axis_label_'+(p)] = axis_label
    locals()['ticks_'+(p)] = ticks

levels = np.linspace(vmin, vmax, num_levels)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
colors = colors(vals)
colors=np.concatenate([[colors[0,:]], colors, [colors[-1,:]]],0)
cmap, norm = from_levels_and_colors(levels, colors, extend='both')

## plot
df = df.sort_values([var_col], ascending=[var_col_sort])

fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_figheight(6)
fig.set_figwidth(8)
for kind in mark_style:
    d = df[df[var_sym]==kind]
    p_tmp = d['Depth_(dbar)']
    SA_tmp = gsw.SA_from_SP(d['Salinity_(psu)'], p_tmp, d['Longitude_(degEast)'], d['Latitude_(degNorth)'])
    CT_tmp = gsw.CT_from_t(SA_tmp,d['Temperature_(degC)'],p_tmp)
    s = plt.scatter(SA_tmp, CT_tmp, marker=mark_style[kind], s=15, lw=0.2, alpha=0.9, edgecolor='black', c=d[var_col].values, cmap=cmap,  vmin=vmin, vmax=vmax, norm=norm)
plt.xlim([vmin_x, vmax_x]); #need to use with water masses figure
plt.ylim([vmin_y, vmax_y]); #need to use with water masses figure
plt.ylabel(axis_label_y, fontsize=11, fontweight='normal')
plt.xlabel(r'$S_A$', fontsize=11, fontweight='normal')
plt.xticks(ticks_x)
plt.yticks(ticks_y)
ax.tick_params(axis='y', direction='out')
ax.tick_params(axis='x', direction='out')
cb = plt.colorbar(s, extend=extent, ticks=ticks)
cb.set_label(axis_label, fontsize=11, fontweight='normal') # <--- plain text

## add isopycnals    
# Figure out boudaries (mins and maxs)
smin = np.nanmin(SA) - (0.01 * np.nanmin(SA))
smax = np.nanmax(SA) + (0.01 * np.nanmax(SA))
tmin = np.nanmin(CT) - (0.01 * np.nanmin(CT))
tmax = np.nanmax(CT) + (0.01 * np.nanmax(CT))
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
if water_masses:
    wm_def = wm.water_masses_def_petrie()
    for w in wm_def.keys():
        shape = np.array(wm_def.get(w))
        p = plt.plot(shape[:,0], shape[:,1])
        co = p[0].get_color() 
        plt.text(np.mean(shape[:,0]), np.mean(shape[:,1]), w, color=co)



## Save
fig.savefig('./AZMP_OA_'+var_y+'_'+var_x+'_'+var_col+'_'+var_sym+'.svg', format='svg', dpi=1000, bbox_inches='tight')
