'''
Map for OA White paper.

Run in: 
/home/cyrf0006/research/carbonSink/OA_WhitePaper

Frederic.Cyr@dfo-mpo.gc.ca
Sept. 2023

'''


import matplotlib.pyplot as plt
import geopandas as gpd
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import os
import pandas as pd

import h5py
import cmocean
import cmocean.cm as cmo
#import cartopy. crs as ccrs
#import cartopy.feature as cpf
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.lines import Line2D
#import cc_variable_list2 as vl
import netCDF4

# Box North America
NAlon = [-170, -20]
#NAlat = [-20, 80]
NAlat = [20, 80]
# Box Pacific
Plon = [-150, -119]
Plat = [45, 60]
# Box Arctic
#Alon = [-100, -45]
#Alat = [50, 75]
Alon = [-110, -41]
Alat = [55, 78]
# Box Atlanttic
Atlon = [-73, -41]
Atlat = [40, 61.5]

# NEXT UPDATE, TRY THIS:
#world = gpd.read_file(gpd.datasets.get_path('/home/cyrf0006/data/country_shapefiles/110m_cultural/ne_110m_admin_0_boundary_lines_land.shp'))

# Deprecated
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world = world.to_crs({'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})
world = world.to_crs({'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})
# Provinces of Canada
province = gpd.read_file('/home/cyrf0006/data/country_shapefiles/lpr_000b16a_e.shp')
tmpWGS84 = province.to_crs({'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})

# drop Canada with low res
world_noCan = world[world.name!='Canada']

# Some plot parameters
lightdeep = cmocean.tools.lighten(cmo.deep, 0.5)
vlightdeep = cmocean.tools.lighten(cmo.deep, 0.75)
ll = np.linspace(0, 2000, 11)
ls = np.linspace(0, 5500, 20)
lss = np.linspace(0, 5000, 11)


### ----- Load station data ----- ##
xls = pd.ExcelFile('./OA_data_resources_v2_AZMP.xlsx')
Obs_point = pd.read_excel(xls, 'Obs_point')
Obs_region = pd.read_excel(xls, 'Obs_region')
Biological_loc = pd.read_excel(xls, 'Biological_loc')
Biological_spec = pd.read_excel(xls, 'Biological_spec')
Obs_AZMP = pd.read_excel(xls, 'Obs_AZMP')
# type str
#Obs_AZMP['Station Name'] = Obs_AZMP['Station Name'].astype('str')

Obs_map = pd.read_excel('OA_data_resources_for_map.xlsx', header=1)
Obs_map = Obs_map.dropna(how='all')
#some dirty cleaning
A =  Obs_map['Station / Asset Name'].str.contains('NOAA')
Obs_map = Obs_map.loc[A.dropna().index]

# Some groups
NOAA = Obs_map.loc[Obs_map['Station / Asset Name'].str.contains('NOAA')]
NOAA = NOAA.loc[NOAA['Unnamed: 0']=='Observational/Monitoring']
NOAA.drop(index=NOAA.index[0], axis=0, inplace=True) # drop Station P
LineP = Obs_map.loc[Obs_map['Station / Asset Name'].str.contains('Line P')]
ONC = Obs_map.loc[Obs_map['Station / Asset Name'].str.contains('ONC')]
#Hakai = Obs_map.loc[Obs_map['Station / Asset Name'].str.contains('Hakai')]
Hakai = Obs_map.loc[Obs_map['Station / Asset Name'].str.contains('KC Buoy')]



###  ----- Bathymetry | North America ----- ###
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
print('Load and grid bathymetry - North America')
decim_scale = 4
# h5 file
h5_outputfile = 'oa_NAmerica_bathymetry.h5'
if os.path.isfile(h5_outputfile):
     print([h5_outputfile + ' exists! Reading directly'])
     h5f = h5py.File(h5_outputfile,'r')
     NAlons = h5f['lon'][:]
     NAlats = h5f['lat'][:]
     NAZ = h5f['Z'][:]
     h5f.close()
else:
    lonLims = [NAlon[0], NAlon[1]]  
    latLims = [NAlat[0], NAlat[1]]
    # Extract variables
    dataset = netCDF4.Dataset(dataFile)
    x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
    y = [-89-59.75/60, 89+59.75/60]
    spacing = dataset.variables['spacing']
    # Compute Lat/Lon
    nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
    ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir
    lon = np.linspace(x[0],x[-1],nx)
    lat = np.linspace(y[0],y[-1],ny)
    ## interpolate data on regular grid (temperature grid)
    ## Reshape data
    zz = dataset.variables['z']
    Z = zz[:].reshape(ny, nx)
    Z = np.flipud(Z) # <------------ important!!!
    # Reduce data according to Region params
    idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
    idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
    Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
    lon = lon[idx_lon[0]]
    lat = lat[idx_lat[0]]
    # Reduce data according to Decim_scale
    lon = lon[::decim_scale]
    lat = lat[::decim_scale]
    Z = Z[::decim_scale, ::decim_scale]
# Save data for later use
    h5f = h5py.File(h5_outputfile, 'w')
    h5f.create_dataset('lon', data=lon)
    h5f.create_dataset('lat', data=lat)
    h5f.create_dataset('Z', data=Z)
    h5f.close()
    NAlons = lon.copy()
    NAlats = lat.copy()
    NAZ = Z.copy()
    print(' -> Done!')


###  ----- Bathymetry | Atlantic ----- ###
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
print('Load and grid bathymetry - Atlantic')
# h5 file
h5_outputfile = 'oa_Atlantic_bathymetry.h5'
if os.path.isfile(h5_outputfile):
     print([h5_outputfile + ' exists! Reading directly'])
     h5f = h5py.File(h5_outputfile,'r')
     Atlons = h5f['lon'][:]
     Atlats = h5f['lat'][:]
     AtZ = h5f['Z'][:]
     h5f.close()
else:
    lonLims = [Atlon[0], Atlon[1]]  
    latLims = [Atlat[0], Atlat[1]]
    # Extract variables
    dataset = netCDF4.Dataset(dataFile)
    x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
    y = [-89-59.75/60, 89+59.75/60]
    spacing = dataset.variables['spacing']
    # Compute Lat/Lon
    nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
    ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir
    lon = np.linspace(x[0],x[-1],nx)
    lat = np.linspace(y[0],y[-1],ny)
    ## interpolate data on regular grid (temperature grid)
    ## Reshape data
    zz = dataset.variables['z']
    Z = zz[:].reshape(ny, nx)
    Z = np.flipud(Z) # <------------ important!!!
    # Reduce data according to Region params
    idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
    idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
    Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
    lon = lon[idx_lon[0]]
    lat = lat[idx_lat[0]]
    # Save data for later use
    h5f = h5py.File(h5_outputfile, 'w')
    h5f.create_dataset('lon', data=lon)
    h5f.create_dataset('lat', data=lat)
    h5f.create_dataset('Z', data=Z)
    h5f.close()
    Atlons = lon.copy()
    Atlats = lat.copy()
    AtZ = Z.copy()
    print(' -> Done!')

###  ----- Bathymetry | Pacific ----- ###
print('Load and grid bathymetry - Pacific')
# h5 file
h5_outputfile = 'oa_Pacific_bathymetry.h5'
if os.path.isfile(h5_outputfile):
     print([h5_outputfile + ' exists! Reading directly'])
     h5f = h5py.File(h5_outputfile,'r')
     Plons = h5f['lon'][:]
     Plats = h5f['lat'][:]
     PZ = h5f['Z'][:]
     h5f.close()
else:
    lonLims = [Plon[0], Plon[1]]  
    latLims = [Plat[0], Plat[1]]
    # Extract variables
    dataset = netCDF4.Dataset(dataFile)
    x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
    y = [-89-59.75/60, 89+59.75/60]
    spacing = dataset.variables['spacing']
    # Compute Lat/Lon
    nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
    ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir
    lon = np.linspace(x[0],x[-1],nx)
    lat = np.linspace(y[0],y[-1],ny)
    ## interpolate data on regular grid (temperature grid)
    ## Reshape data
    zz = dataset.variables['z']
    Z = zz[:].reshape(ny, nx)
    Z = np.flipud(Z) # <------------ important!!!
    # Reduce data according to Region params
    idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
    idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
    Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
    lon = lon[idx_lon[0]]
    lat = lat[idx_lat[0]]
    # Save data for later use
    h5f = h5py.File(h5_outputfile, 'w')
    h5f.create_dataset('lon', data=lon)
    h5f.create_dataset('lat', data=lat)
    h5f.create_dataset('Z', data=Z)
    h5f.close()
    Plons = lon.copy()
    Plats = lat.copy()
    PZ = Z.copy()
    print(' -> Done!')


###  ----- Bathymetry | Arctic ----- ###
print('Load and grid bathymetry - Arctic')
# h5 file
h5_outputfile = 'oa_Arctic_bathymetry.h5'
if os.path.isfile(h5_outputfile):
     print([h5_outputfile + ' exists! Reading directly'])
     h5f = h5py.File(h5_outputfile,'r')
     Alons = h5f['lon'][:]
     Alats = h5f['lat'][:]
     AZ = h5f['Z'][:]
     h5f.close()
else:
    lonLims = [Alon[0], Alon[1]]  
    latLims = [Alat[0], Alat[1]]
    # Extract variables
    dataset = netCDF4.Dataset(dataFile)
    x = [-179-59.75/60, 179+59.75/60] # to correct bug in 30'' dataset?
    y = [-89-59.75/60, 89+59.75/60]
    spacing = dataset.variables['spacing']
    # Compute Lat/Lon
    nx = int((x[-1]-x[0])/spacing[0]) + 1  # num pts in x-dir
    ny = int((y[-1]-y[0])/spacing[1]) + 1  # num pts in y-dir
    lon = np.linspace(x[0],x[-1],nx)
    lat = np.linspace(y[0],y[-1],ny)
    ## interpolate data on regular grid (temperature grid)
    ## Reshape data
    zz = dataset.variables['z']
    Z = zz[:].reshape(ny, nx)
    Z = np.flipud(Z) # <------------ important!!!
    # Reduce data according to Region params
    idx_lon = np.where((lon>=lonLims[0]) & (lon<=lonLims[1]))
    idx_lat = np.where((lat>=latLims[0]) & (lat<=latLims[1]))
    Z = Z[idx_lat[0][0]:idx_lat[0][-1]+1, idx_lon[0][0]:idx_lon[0][-1]+1]
    lon = lon[idx_lon[0]]
    lat = lat[idx_lat[0]]
    # Save data for later use
    h5f = h5py.File(h5_outputfile, 'w')
    h5f.create_dataset('lon', data=lon)
    h5f.create_dataset('lat', data=lat)
    h5f.create_dataset('Z', data=Z)
    h5f.close()
    Alons = lon.copy()
    Alats = lat.copy()
    AZ = Z.copy()
    print(' -> Done!')
    
### --- Main Plot --- ###
fig, ax = plt.subplots(nrows=1, ncols=1)
c = ax.contourf(NAlons, NAlats, -NAZ, lss, cmap=lightdeep, extend='max', zorder=5)
world.plot(ax=ax, color='tan')

#ax.set_xlim([-75, 38])b
ax.set_xlim([NAlon[0], NAlon[1]])
ax.set_ylim([NAlat[0], NAlat[1]])
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax.xaxis.set_label_position('top')

plt.xlabel(r'Longitude ($\rm ^{\circ}E$)', fontsize=13)
plt.ylabel(r'Latitude ($\rm ^{\circ}N$)', fontsize=13)
cax = fig.add_axes([0.91, .4, 0.015, 0.2])
cb = plt.colorbar(c, cax=cax, orientation='vertical')
cb.set_label('Depth (m)', fontsize=12, fontweight='normal')

## --- Add insets --- ##
# Define news axes
#axins = zoomed_inset_axes(ax, zoom=2, loc=3)
#axins = zoomed_inset_axes(ax, zoom=2.4, loc=3, bbox_to_anchor=(55,135))
axins = zoomed_inset_axes(ax, zoom=3, loc=3, bbox_to_anchor=(-30,50))
#axins2 = zoomed_inset_axes(ax, zoom=2, loc=8)
axins2 = zoomed_inset_axes(ax, zoom=2.3, loc=8, bbox_to_anchor=(450,50))
#axins3 = zoomed_inset_axes(ax, zoom=2, loc=4)
#axins3 = zoomed_inset_axes(ax, zoom=2.3, loc=4, bbox_to_anchor=(1190,135))
axins3 = zoomed_inset_axes(ax, zoom=2.7, loc=4, bbox_to_anchor=(930,50))
#bbox_to_anchor=(600, 0, 400, 0)

# Box Pacific
axins.set_xlim(Plon[0], Plon[1])
axins.set_ylim(Plat[0], Plat[1])
# Box Arctic
axins2.set_xlim(Alon[0], Alon[1])
axins2.set_ylim(Alat[0], Alat[1])
# Box Atlanttic
axins3.set_xlim(Atlon[0], Atlon[1])
axins3.set_ylim(Atlat[0], Atlat[1])

mark_inset(ax, axins, loc1=1, loc2=2, fc="none", ec="0", zorder=100)
mark_inset(ax, axins2, loc1=1, loc2=2, fc="none", ec="0", zorder=100)
mark_inset(ax, axins3, loc1=1, loc2=2, fc="none", ec="0", zorder=100)
axins.set_zorder(100)
axins2.set_zorder(100)
axins3.set_zorder(100)
axins.set_title('Pacific', color='k')
axins2.set_title('Arctic', color='k')
axins3.set_title('Atlantic', color='k')

# Manual Legend
handles, labels = ax.get_legend_handles_labels()
line1 = Line2D([0], [0], label='Core AZMP section', color='salmon')
line2 = Line2D([0], [0], label='Other DFO section', color='r')
point1 = Line2D([0], [0], label='DFO HR Station', marker='*', markersize=10, markeredgecolor='m', markerfacecolor='m', linestyle='')
point2 = Line2D([0], [0], label='NOAA Station', marker='*', markersize=10, markeredgecolor='tab:blue', markerfacecolor='tab:blue', linestyle='')
point3 = Line2D([0], [0], label='ONC Station', marker='*', markersize=10, markeredgecolor='tab:green', markerfacecolor='tab:green', linestyle='')
point4 = Line2D([0], [0], label='Hakai Station', marker='*', markersize=10, markeredgecolor='tab:orange', markerfacecolor='tab:orange', linestyle='')
handles.extend([line1, line2, point1, point2, point3, point4])
ax.legend(handles=handles, loc='lower center')

## ---- Plot zooms ---- ##
## 1. Pacific
c = axins.contourf(Plons, Plats, -PZ, ls, cmap=lightdeep, extend='max', zorder=5)
#cc = axins.contour(lon, lat, -Z, [100, 500, 1000, 2000, 3000, 4000, 5000], colors='silver', linewidths=.5, zorder=10)
world_noCan.plot(color='tan',ax=axins, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'British Columbia / Colombie-Britannique'].plot(color='tan',ax=axins)
tmpWGS84[tmpWGS84.PRNAME == 'Alberta'].plot(color='tan',ax=axins)
#world.plot(edgecolor='black', color='tan',ax=axins, alpha=1)


## 2. Arctic
c = axins2.contourf(Alons, Alats, -AZ, ls, cmap=lightdeep, extend='max', zorder=5)
#cc = axins2.contour(lon, lat, -Z, [100, 500, 1000, 2000, 3000, 4000, 5000], colors='silver', linewidths=.5, zorder=10)
world_noCan.plot(color='tan',ax=axins2, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'Newfoundland and Labrador / Terre-Neuve-et-Labrador'].plot(color='tan',ax=axins2, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'Nova Scotia / Nouvelle-Écosse'].plot(color='tan',ax=axins2, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'New Brunswick / Nouveau-Brunswick'].plot(color='tan',ax=axins2, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'Prince Edward Island / Île-du-Prince-Édouard'].plot(color='tan',ax=axins2)
tmpWGS84[tmpWGS84.PRNAME == 'Quebec / Québec'].plot(color='tan',ax=axins2, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'Ontario'].plot(color='tan',ax=axins2)
tmpWGS84[tmpWGS84.PRNAME == 'Manitoba'].plot(color='tan',ax=axins2)
tmpWGS84[tmpWGS84.PRNAME == 'Saskatchewan'].plot(color='tan',ax=axins2)
tmpWGS84[tmpWGS84.PRNAME == 'Alberta'].plot(color='tan',ax=axins2)
tmpWGS84[tmpWGS84.PRNAME == 'Yukon'].plot(color='tan',ax=axins2)
tmpWGS84[tmpWGS84.PRNAME == 'Northwest Territories / Territoires du Nord-Ouest'].plot(color='tan',ax=axins2)
tmpWGS84[tmpWGS84.PRNAME == 'Nunavut'].plot(color='tan',ax=axins2)
#world.plot(edgecolor='black', color='tan',ax=axins2, alpha=1)


## 3. Atlantic
c = axins3.contourf(Atlons, Atlats, -AtZ, ls, cmap=lightdeep, extend='max', zorder=5)
#cc = axins3.contour(lon, lat, -Z, [100, 500, 1000, 2000, 3000, 4000, 5000], colors='silver', linewidths=.5, zorder=10)
world_noCan.plot(color='tan',ax=axins3, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'Newfoundland and Labrador / Terre-Neuve-et-Labrador'].plot(color='tan',ax=axins3, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'Nova Scotia / Nouvelle-Écosse'].plot(color='tan',ax=axins3, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'New Brunswick / Nouveau-Brunswick'].plot(color='tan',ax=axins3, alpha=1)
tmpWGS84[tmpWGS84.PRNAME == 'Prince Edward Island / Île-du-Prince-Édouard'].plot(color='tan',ax=axins3)
tmpWGS84[tmpWGS84.PRNAME == 'Quebec / Québec'].plot(color='tan',ax=axins3, alpha=1)
#world.plot(edgecolor='black', color='tan',ax=axins3, alpha=1)

## 4. Add stations
# AZMP
AZMP_sections = Obs_AZMP.loc[Obs_AZMP.type=='Section']
BANQ = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('BANQ')]
BB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('BB-')]
BBL = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('BBL_')]
BP = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('BP_')]
CMO = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('CMO')]
CSL = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('CSL')]
FC = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('FC-')]
HL = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('HL_')]
LCM = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('LCM')] #LC Mouth
LL = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('LL')]
NEC = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('NEC')]
PL = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('PL')] #GoM
PS = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('PS')]
SAG = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('SAG')]
SEGB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('SEGB-')]
SESPB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('SESPB-')]
SWSPB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('SWSPB-')]
SI = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('SI-')]
TCEN = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('TCEN')]
TDC = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('TDC')] #Cabot St.
TIDM = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('TIDM')]
TSI = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('TSI')] #SeptIles
YL = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('YL')] #GoM
MB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('MB-')]
BI = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('BI-')]
WB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('WB-')]
STAB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('STAB')]
TASO = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('TASO')]
TBB = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('TBB')]
# Those with wrong order (ignore anyway, not aligned)
CH = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('CH')]
CH = CH.iloc[[9,12],]
IF = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('IF')]
IF = IF.iloc[[0,10],]
# Those not in straight line
LC = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('LC-')]
XGSL = AZMP_sections.loc[AZMP_sections['Station Name'].str.contains('XGSL-')]


sections_toplot = ['BB', 'BBL', 'BP', 'CSL', 'FC', 'HL', 'LCM', 'LL', 'NEC', 'PL', 'SAG', 'SEGB', 'SESPB', 'SWSPB', 'SI', 'TDC', 'TIDM', 'TSI', 'YL', 'BI', 'MB', 'WB', 'STAB', 'TASO', 'TBB']
AZMP_sections = ['BB', 'BBL', 'CSL', 'FC', 'HL', 'LL', 'SEGB', 'SESPB', 'SWSPB', 'SI', 'TIDM', 'TSI', 'BI', 'MB', 'WB', 'TASO', 'TBB']
sections_toignore = ['PS', 'BANQ', 'TCEN', 'CMO', 'CH', 'IF']

for sec in sections_toplot:
    command = "axins3.plot(" + sec + ".iloc[[0,-1],].Longitude.values, " + sec + ".iloc[[0,-1],].Latitude.values, 'r', zorder=10)"
    eval(command)
for sec in AZMP_sections:
    command = "axins3.plot(" + sec + ".iloc[[0,-1],].Longitude.values, " + sec + ".iloc[[0,-1],].Latitude.values, 'salmon', zorder=10)"
    eval(command)


# Non-straight transects
axins3.plot(LC.Longitude.values, LC.Latitude.values, 'r', zorder=10)
axins3.plot(XGSL.Longitude.values, XGSL.Latitude.values, 'r', zorder=10)

# HR Stations
AZMP_HR = Obs_AZMP.loc[Obs_AZMP.type=='Hrstation']
axins3.plot(AZMP_HR[AZMP_HR['Station Name']=='HL_02'].Longitude.values, AZMP_HR[AZMP_HR['Station Name']=='HL_02'].Latitude.values, marker='*', color='m', zorder=10)
axins3.plot(AZMP_HR[AZMP_HR['Station Name']=='Rimouski'].Longitude.values, AZMP_HR[AZMP_HR['Station Name']=='Rimouski'].Latitude.values, marker='*', color='m', zorder=10)
axins3.plot(AZMP_HR[AZMP_HR['Station Name']=='STN27'].Longitude.values, AZMP_HR[AZMP_HR['Station Name']=='STN27'].Latitude.values, marker='*', color='m', zorder=10)
axins3.annotate(' Halifax 2', xy=[AZMP_HR[AZMP_HR['Station Name']=='HL_02'].Longitude.values, AZMP_HR[AZMP_HR['Station Name']=='HL_02'].Latitude.values], horizontalalignment='left', verticalalignment='top', color='black', zorder=10)
axins3.annotate(' Rimouski', xy=[AZMP_HR[AZMP_HR['Station Name']=='Rimouski'].Longitude.values, AZMP_HR[AZMP_HR['Station Name']=='Rimouski'].Latitude.values-1], horizontalalignment='center', verticalalignment='top', color='black', zorder=10)
axins3.annotate(' Station 27', xy=[AZMP_HR[AZMP_HR['Station Name']=='STN27'].Longitude.values, AZMP_HR[AZMP_HR['Station Name']=='STN27'].Latitude.values], horizontalalignment='left', verticalalignment='bottom',  color='black', zorder=10)


# b) Line P - Pacific
axins.plot(-LineP['Longitude (°W)'].values, LineP['Latitude (°N)'].values, color='r', linewidth=1, zorder=10)
axins.plot(-LineP['Longitude (°W)'].iloc[-1], LineP['Latitude (°N)'].iloc[-1], color='m', marker='*', linewidth=0, zorder=10)
axins.annotate('Station P', xy=[-LineP['Longitude (°W)'].iloc[-1], LineP['Latitude (°N)'].iloc[-1]], horizontalalignment='center', verticalalignment='bottom', color='black', zorder=10)

# c) NOAA Atlantic - Pacific
axins3.plot(-NOAA['Longitude (°W)'].values, NOAA['Latitude (°N)'].values, color='tab:blue', marker='*', linewidth=0, zorder=10)
axins.plot(-NOAA['Longitude (°W)'].values, NOAA['Latitude (°N)'].values, color='tab:blue', marker='*', linewidth=0, zorder=10)

# d) ONC - Arctic Pacific
axins.plot(-ONC['Longitude (°W)'].values, ONC['Latitude (°N)'].values, color='tab:green', marker='*', linewidth=0, zorder=10)
axins2.plot(-ONC['Longitude (°W)'].values, ONC['Latitude (°N)'].values, color='tab:green', marker='*', linewidth=0, zorder=10)

# e) Hakai - Pacific
axins.plot(Hakai['Longitude (°W)'].values, Hakai['Latitude (°N)'].values, color='tab:orange', marker='*', linewidth=0, zorder=10)

# f) DFO Arctic Transect (manual from map)
axins3.plot([-55.660, -48.330], [53.570, 60.570], 'r', zorder=10) #Ar7W
axins2.plot([-61.260, -53.950], [66.650, 67.310], 'r', zorder=10) # Davis St.
axins2.plot([-90.440, -90.140], [74.090, 74.540], 'r', zorder=10) # Barrow St.
#axins3.plot([], [], 'r', zorder=10)


plt.setp(axins.get_xticklabels(), visible=False)
plt.setp(axins.get_yticklabels(), visible=False)
plt.setp(axins2.get_xticklabels(), visible=False)
plt.setp(axins2.get_yticklabels(), visible=False)
plt.setp(axins3.get_xticklabels(), visible=False)
plt.setp(axins3.get_yticklabels(), visible=False)


fig.set_size_inches(w=8.5, h=9.5)
#fig.tight_layout() 
fig.set_dpi(300)
figname = 'map_Oa_WhitePaper.png'
fig.savefig(figname)
os.system('convert -trim ' + figname + ' ' + figname)

