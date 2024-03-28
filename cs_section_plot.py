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
import gsw

## ----- Parameters to edit ---- ##
# Variable you want to plot
#VAR = 'Omega_Aragonite_(unitless)'
#VAR = 'pH_Total_(total_scale)'
VAR = 'Oxygen_Saturation_(%)'
#VAR = 'Dissolved_Oxygen_(mL/L)'
#VAR = 'Phosphate_Concentration_(mmol/m3)'
#VAR = 'Temperature_(degC)'
#VAR = 'pH_tot'
#VAR = 'Omega_Aragonite_(--)'
#VAR = 'Salinity_(psu)'
#VAR = 'Total_Alkalinity_(umol/kg)'
VAR = 'Inorganic_Carbon_(umol/kg)'
#VAR = 'AOU'
REGION = 'NL'
SECTION = 'FC'
SEASON = 'spring'
YEAR = 2016# plotting parameters
#v = 10
#v_anom=10
CMAP = cmo.cm.seismic
ZMAX = 1000
# File to load
my_file = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv'
# For colorbar:
if REGION == 'NL':
    import cc_variable_list_NL as vl
else:
    import cc_variable_list2 as vl
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
section_EN = ['BI', 'MB', 'SI', 'WB', 'BB', 'S27', 'FC', 'SEGB', 'SWSPB']
section_FR = ['IB', 'MB', 'IS', 'WB', 'BB', 'S27', 'BF', 'GBSE', 'SWSPB']
SECTION_FR = section_FR[section_EN.index(SECTION)]
season_EN = ['spring', 'summer', 'fall']
season_FR = ['printemps', 'été', 'automne']
SEASON_FR = season_FR[season_EN.index(SEASON)]
VAR_text = VAR.split('(')[0][0:-1] 
if VAR == 'Oxygen_Saturation_(%)':
    title = 'Oxygen Saturation for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
    title_FR = 'Saturation en oxygène pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
elif VAR == 'pH_Total_(total_scale)':
    title = 'pH Total for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
    title_FR = 'pH Total pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
elif VAR == 'Omega_Aragonite_(unitless)':
    title = 'Aragonite saturation state for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
    title_FR = 'Saturation en aragonite pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
elif VAR == 'AOU':
    VAR_text = 'AOU'
    title = 'AOU for section ' + SECTION + ' - ' + SEASON + ' ' + str(YEAR)
    title_FR = 'UAO pour la section ' + SECTION_FR + ' - ' + SEASON_FR + ' ' + str(YEAR)
else:
    title=VAR_text
    title_FR=VAR_text
    
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
# Extract section
#df = df.loc[(df.Station_Name.str.contains(SECTION + '-')) | (df.Station_Name.str.contains(SECTION + '_')) ]
df = df.loc[(df.Station_Name.str.contains(SECTION))]
# Extract season
#df = df[df.index.year == YEAR]
if SEASON == 'spring':
    df = df[df.index.month <= 5]
elif SEASON == 'summer':
    df = df[(df.index.month>=6) & (df.index.month<=8)]
elif SEASON == 'fall':
    df = df[df.index.month >= 9]

# Calculate Apparent Oxygen Utilisation (AOU) if needed
if VAR == 'AOU':
    SA = gsw.SA_from_SP(df['Salinity_(psu)'], df['Depth_(dbar)'], df['Longitude_(degEast)'], df['Latitude_(degNorth)'])
    CT = gsw.CT_from_t(SA, df['Temperature_(degC)'], df['Depth_(dbar)'])
    O2sol = gsw.O2sol(SA, CT, df['Depth_(dbar)'],  df['Longitude_(degEast)'], df['Latitude_(degNorth)']) # in umol/kg
    O2sol = O2sol/43.570 # in ml/l 
    df['AOU'] = O2sol - df['Dissolved_Oxygen_(mL/L)'] 
    
# Build climatology
df_list = []
for i in df.index.year.unique():
    df_tmp = df[df.index.year == i]
    # Extract variable
    df_tmp = df_tmp[['Depth_(dbar)', 'Station_Name', VAR]]
    df_tmp = df_tmp.pivot(index='Depth_(dbar)', columns='Station_Name') #<--- this is cool!
    df_tmp = df_tmp[VAR]
    # So I re-define a constant vertical axis every 5m.
    depth_range = np.arange(2.5, 2000, 5) # range to look for data
    reg_depth = (depth_range[1:] + depth_range[:-1]) / 2 # mid point of the range
    df_tmp = df_tmp.groupby(pd.cut(df_tmp.index, depth_range)).mean() # <--- This is cool!
    df_tmp.index = reg_depth # replace range by mean depth
    # interpolate vertically and horisontally where possible
    df_tmp.interpolate(axis=0, limit_area='inside', inplace=True)
    df_tmp.interpolate(axis=1, limit_area='inside', inplace=True)
    # Drop depth where there is no values associated to a certain depth (because my axis is every 5m...)
    #df_tmp = df_tmp.dropna(how='all')
    df_list.append(df_tmp)
    #del df_tmp
    

# Create multi-index
df_all = pd.concat(df_list, keys=df.index.year.unique(), sort=True)

# extract year current year
df_year = df_all.xs((YEAR), level=('Timestamp'))

# compute climatology
df_clim = df_all.groupby(level=1).apply(lambda x: x.mean())

# vertically fill NaNs
df_year.interpolate(axis=0, limit_area='inside', inplace=True)
df_clim.interpolate(axis=0, limit_area='inside', inplace=True)
# horizontally fill NaNs
df_year.interpolate(axis=1, limit_area='inside', inplace=True)
df_clim.interpolate(axis=1, limit_area='inside', inplace=True)

# calculate anomaly
df_anom = df_year-df_clim

## ---- Load station lat/lon ---- ##
df_stn = pd.read_excel('/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx')
df_stn = df_stn.drop(['SECTION', 'LONG'], axis=1)
df_stn = df_stn.rename(columns={'LONG.1': 'LON'})
df_stn = df_stn.dropna()
df_stn = df_stn[df_stn.STATION.str.contains(SECTION)]
df_stn = df_stn.reset_index(drop=True)
# remove hyphen in stn names (normally should keep it everywhere, but simpler here)
#df_stn['STATION'] = df_stn.STATION.str.replace('-', '')

## ---- Compute distance vector ---- ##
distance = np.full((df_clim.keys().shape), np.nan)
lat0 = df_stn[df_stn.index==0]['LAT']
lon0 = df_stn[df_stn.index==0]['LON']
for i, stn in enumerate(df_clim.keys().values):
    stn = stn.replace('_','-') # replace in case underscore is used
    lat_stn = df_stn[df_stn.STATION==str(stn)]['LAT']
    lon_stn = df_stn[df_stn.STATION==str(stn)]['LON']      
    distance[i] = azst.haversine(lon0, lat0, lon_stn, lat_stn)
XLIM =distance.max()
# Distance for current year
## distance_year = np.full((df_year.index.shape), np.nan)
## for i, stn in enumerate(df_year.index):
##     lat_stn = df_stn[df_stn.STATION==stn]['LAT']
##     lon_stn = df_stn[df_stn.STATION==stn]['LON']      
##     distance_year[i] = haversine(lon0, lat0, lon_stn, lat_stn)
## # Distance for sigt (rare but sometimes not same as other)
## distance_sigt = np.full((df_sigt_year.index.shape), np.nan)
## for i, stn in enumerate(df_sigt_year.index):
##     lat_stn = df_stn[df_stn.STATION==stn]['LAT']
##     lon_stn = df_stn[df_stn.STATION==stn]['LON']      
##     distance_sigt[i] = haversine(lon0, lat0, lon_stn, lat_stn)



## ---- Retrieve bathymetry using function ---- ##
if REGION == 'NL':
    bathymetry = azst.section_bathymetry(SECTION)
else:
    bathymetry = []
## ---- plot Figure ---- ##
#XLIM = df_section_itp.index[-1][1]
fig = plt.figure()
# ax1
ax = plt.subplot2grid((3, 1), (0, 0))
c = plt.contourf(distance, df_year.index, df_year, levels, cmap=cmap, extend='both')
#c_sig1 = plt.contour(distance_sigt, df_sigt_year.columns, df_sigt_year, v_sig, colors='gray', linewidths=1)
for i in distance:
    plt.plot([i, i], [0, ZMAX], '--k', alpha=.5)
ax.set_ylim([0, ZMAX])
ax.set_xlim([0, XLIM])
#plt.clabel(c_sig1, inline=1, fontsize=10, colors='gray', fmt='%1.1f')
ax.set_ylabel('Depth (m)', fontweight = 'bold')
ax.invert_yaxis()
if REGION == 'NL':
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
    ax.add_patch(Bgon)
cb = plt.colorbar(c)
cb.ax.tick_params(labelsize=8)
cb.set_label(axis_label, fontsize=12, fontweight='normal')
ax.xaxis.label.set_visible(False)
ax.tick_params(labelbottom='off')
ax.set_title(title)

# ax2
ax2 = plt.subplot2grid((3, 1), (1, 0))
c = plt.contourf(distance, df_clim.index, df_clim, levels, cmap=cmap, extend='both')
#c_sig2 = plt.contour(distance, df_sigt_clim.columns, df_sigt_clim, v_sig, colors='gray', linewidths=1)
for i in distance:
    plt.plot([i, i], [0, ZMAX], '--k', alpha=.5)
ax2.set_ylim([0, ZMAX])
ax2.set_xlim([0,  XLIM])
#plt.clabel(c_sig2, inline=1, fontsize=10, colors='gray', fmt='%1.1f')
ax2.set_ylabel('Depth (m)', fontweight = 'bold')
ax2.invert_yaxis()
if REGION == 'NL':
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
    ax2.add_patch(Bgon)
cb = plt.colorbar(c)
cb.ax.tick_params(labelsize=8)
cb.set_label(axis_label, fontsize=12, fontweight='normal')
ax2.xaxis.label.set_visible(False)
ax2.tick_params(labelbottom='off')
ax2.set_title('Climatology (2014-2020)')

# ax3
ax3 = plt.subplot2grid((3, 1), (2, 0))
#c = plt.contourf(distance, df_anom.index, df_anom, v_anom, cmap=cmo.cm.seismic, extend='both')
c = plt.contourf(distance, df_anom.index, df_anom, v_anom, cmap=plt.cm.RdBu_r, extend='both')
ax3.set_ylim([0, ZMAX])
ax3.set_xlim([0,  XLIM])
ax3.set_ylabel('Depth (m)', fontweight = 'bold')
ax3.set_xlabel('Distance (km)', fontweight = 'bold')
ax3.invert_yaxis()
if REGION == 'NL':
    Bgon = plt.Polygon(bathymetry,color=np.multiply([1,.9333,.6667],.4), alpha=1, zorder=10)
    ax3.add_patch(Bgon)
cb = plt.colorbar(c)
cb.ax.tick_params(labelsize=8)
cb.set_label(axis_label, fontsize=12, fontweight='normal')
ax3.set_title(r'Anomaly')

fig.set_size_inches(w=8,h=12)
fig_name = 'btl_' + VAR_text + '_' + SECTION + '_' + SEASON + '_' + str(YEAR) + '.png' 
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)


## French figure
ax.set_title(title_FR)
ax2.set_title('Climatologie (2014-2020)')
ax3.set_title(r'Anomalie')
ax.set_ylabel('Profondeur (m)', fontweight = 'bold')
ax2.set_ylabel('Profondeur (m)', fontweight = 'bold')
ax3.set_ylabel('Profondeur (m)', fontweight = 'bold')
fig.set_size_inches(w=8,h=12)
fig_name = 'btl_' + VAR_text + '_' + SECTION + '_' + SEASON + '_' + str(YEAR) + '_FR.png' 
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)

# For NL ResDoc
 ## montage btl_Omega_Aragonite_SI_summer_2019.png btl_Omega_Aragonite_SI_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_SI_2019-2020_summer.png
 ## montage btl_Omega_Aragonite_BB_summer_2019.png btl_Omega_Aragonite_BB_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_BB_2019-2020_summer.png
 ## montage btl_Omega_Aragonite_FC_summer_2019.png btl_Omega_Aragonite_FC_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_FC_2019-2020_summer.png

 ## montage btl_Oxygen_Saturation_SI_summer_2019.png btl_Oxygen_Saturation_SI_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_SI_2019-2020_summer.png
 ## montage btl_Oxygen_Saturation_BB_summer_2019.png btl_Oxygen_Saturation_BB_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_BB_2019-2020_summer.png
 ## montage btl_Oxygen_Saturation_FC_summer_2019.png btl_Oxygen_Saturation_FC_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_FC_2019-2020_summer.png

 ## montage btl_pH_Total_SI_summer_2019.png btl_pH_Total_SI_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_SI_2019-2020_summer.png
 ## montage btl_pH_Total_BB_summer_2019.png btl_pH_Total_BB_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_BB_2019-2020_summer.png
 ## montage btl_pH_Total_FC_summer_2019.png btl_pH_Total_FC_summer_2020.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_FC_2019-2020_summer.png

## In French:
## montage btl_Omega_Aragonite_SI_summer_2019_FR.png btl_Omega_Aragonite_SI_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_SI_2019-2020_summer_FR.png
## montage btl_Omega_Aragonite_BB_summer_2019_FR.png btl_Omega_Aragonite_BB_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_BB_2019-2020_summer_FR.png
## montage btl_Omega_Aragonite_FC_summer_2019_FR.png btl_Omega_Aragonite_FC_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Omega_Aragonite_FC_2019-2020_summer_FR.png

## montage btl_Oxygen_Saturation_SI_summer_2019_FR.png btl_Oxygen_Saturation_SI_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_SI_2019-2020_summer_FR.png
## montage btl_Oxygen_Saturation_BB_summer_2019_FR.png btl_Oxygen_Saturation_BB_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_BB_2019-2020_summer_FR.png
## montage btl_Oxygen_Saturation_FC_summer_2019_FR.png btl_Oxygen_Saturation_FC_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_Oxygen_Saturation_FC_2019-2020_summer_FR.png

## montage btl_pH_Total_SI_summer_2019_FR.png btl_pH_Total_SI_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_SI_2019-2020_summer_FR.png
## montage btl_pH_Total_BB_summer_2019_FR.png btl_pH_Total_BB_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_BB_2019-2020_summer_FR.png
## montage btl_pH_Total_FC_summer_2019_FR.png btl_pH_Total_FC_summer_2020_FR.png -tile 2x1 -geometry +10+10  -background white NL_OA_pH_Total_FC_2019-2020_summer_FR.png
