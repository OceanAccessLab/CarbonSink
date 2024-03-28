'''
This script generates the maps for AZMP carbon project with B. Zemskova.

Run in ~/research/carbonSink/Zemskova

Frederic.Cyr@dfo-mpo.gc.ca
March 2021
Update Fall 2023 with new data

TODO: redo this figure after we decide the final classification
One figure per season.


'''


import os
import h5py                                                                
import netCDF4
import numpy as  np
import matplotlib.pyplot as plt
import openpyxl, pprint
import pandas as pd
import cmocean as cmo
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from cartopy.crs import EqualEarth, PlateCarree
import cartopy.crs as ccrs
import cartopy.feature as cpf
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

## ---- Region parameters ---- ##
dataFile = '/home/cyrf0006/data/GEBCO/GEBCO_2014_1D.nc'
lon_0 = -50
lat_0 = 50
lonLims = [-65, -39.9]
latLims = [40, 60]
proj = 'merc'
decim_scale = 4
stationFile = '/home/cyrf0006/github/AZMP-NL/data/STANDARD_SECTIONS.xlsx'
lightdeep = cmo.tools.lighten(cmo.cm.deep, .7)
v = np.linspace(0, 3500, 36)
#v = np.linspace(-4000, 0, 9)
#v = np.linspace(0, 5400, 55) # Olivia's (good one, but lose channels on the shelf)
v = np.linspace(0, 4000, 41) # t
zlims = [0, 25]
#zlims = [75, 125]
#zlims = [150, 300]
zlims = [300, 500]
#zlims = [500, 1000]
#zlims = [1000, 5000]

class_type = '6_classes'
TSDIC = False
#season = 'spring'
rand_range = [-.4, .4]
# Default colors cycle
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
if TSDIC:
     file_name = './Testing_classification/NLonly_data_CT_SA_DIC_classes.csv'
     fig_name = 'map_carbon_wm_' + class_type + '_' + str(zlims[0]) + '-' + str(zlims[1]) + 'm.png'
else:
     file_name = './Testing_classification/NLonly_data_CT_SA_classes.csv'
     fig_name = 'map_carbon_wmTS_' + class_type + '_' + str(zlims[0]) + '-' + str(zlims[1]) + 'm.png'
     
## ---- Carbon data ---- ##
df_wm = pd.read_csv(file_name, delimiter=',')
df_wm.set_index('Timestamp', inplace=True)
df_wm.index = pd.to_datetime(df_wm.index)
df_wm = df_wm[(df_wm['Depth_(dbar)']>=zlims[0]) & (df_wm['Depth_(dbar)']<=zlims[1])]
wm_lat = df_wm['Latitude_(degNorth)']
wm_lon = df_wm['Longitude_(degEast)']
# add random perturbation to latitude (for better display on map)
#wm_lat = wm_lat + np.random.random(len(wm_lat))*(rand_range[1]-rand_range[0]) + rand_range[0] + wm_lon.index.month/10
wm_lat = wm_lat + np.random.random(len(wm_lat))*(rand_range[1]-rand_range[0])/5 + wm_lon.index.month/5 - .6
#wm_lat = wm_lat + wm_lon.index.month/10

## ---- Bathymetry ---- ##
print('Load and grid bathymetry')
# h5 file
#h5_outputfile = '/home/cyrf0006/AZMP/utils/nl_stacfis_bathymetry.h5'
h5_outputfile = '/home/cyrf0006/AZMP/utils/nl_climate_bathymetry.h5'

if os.path.isfile(h5_outputfile):
     print([h5_outputfile + ' exists! Reading directly'])
     h5f = h5py.File(h5_outputfile,'r')
     lon = h5f['lon'][:]
     lat = h5f['lat'][:]
     Z = h5f['Z'][:]
     h5f.close()

else:
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
    # Reshape data
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
    print(' -> Done!')


## ---- Station info ---- ##
import pandas as pd
df = pd.read_excel(stationFile)
sections = df['SECTION'].values
stations = df['STATION'].values
stationLat = df['LAT'].values
stationLon = df['LONG.1'].values

index_BI = df.SECTION[df.SECTION=="BEACH ISLAND"].index.tolist()
index_MB = df.SECTION[df.SECTION=="MAKKOVIK BANK"].index.tolist()
index_SI = df.SECTION[df.SECTION=="SEAL ISLAND"].index.tolist()
index_WB = df.SECTION[df.SECTION=="WHITE BAY"].index.tolist()
index_BB = df.SECTION[df.SECTION=="BONAVISTA"].index.tolist()
index_FC = df.SECTION[df.SECTION=="FLEMISH CAP"].index.tolist()
index_SEGB = df.SECTION[df.SECTION=="SOUTHEAST GRAND BANK"].index.tolist()
index_SESPB = df.SECTION[df.SECTION=="SOUTHEAST ST PIERRE BANK"].index.tolist()
index_SWSPB = df.SECTION[df.SECTION=="SOUTHWEST ST PIERRE BANK"].index.tolist()

index_S27 = df.SECTION[df.SECTION=="STATION 27"].index.tolist()


## ---- plot ---- ##
x,y = lon[0::4],lat[0::4]
Z = Z[0::4,0::4]
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111, projection=ccrs.Mercator())
ax.set_extent([lonLims[0], lonLims[1], latLims[0], latLims[1]], crs=ccrs.PlateCarree())
ax.add_feature(cpf.NaturalEarthFeature('physical', 'coastline', '50m', edgecolor='k', alpha=0.7, linewidth=0.6, facecolor='tan'), zorder=15)#cpf.COLORS['land']))
m=ax.gridlines(linewidth=0.5, color='black', draw_labels=True, alpha=0.5)
m.top_labels = False
m.right_labels = False
m.xlocator = mticker.FixedLocator([-75, -70, -60, -50, -40])
m.ylocator = mticker.FixedLocator([40, 45, 50, 55, 60, 65])
m.xformatter = LONGITUDE_FORMATTER
m.yformatter = LATITUDE_FORMATTER
m.ylabel_style = {'size': 7, 'color': 'black', 'weight':'bold'}
m.xlabel_style = {'size': 7, 'color': 'black', 'weight':'bold'}
lightdeep = cmo.tools.lighten(cmo.cm.deep, 0.5)
c = plt.contourf(x, y, -Z, v, transform=ccrs.PlateCarree(), cmap=lightdeep, extend='max', zorder=1)
cc = plt.contour(x, y, -Z, [100, 500, 1000, 2000, 3000, 4000], transform=ccrs.PlateCarree(), colors='lightgrey', linewidths=.5, zorder=1)
#ccc = plt.contour(x, y, -Z, [1000], colors='darkgrey', linewidths=2, zorder=1)
plt.clabel(cc, inline=1, fontsize=10, colors='lightgrey', fmt='%d')


#c = plt.contourf(lon, lat, -Z, ls, transform=ccrs.PlateCarree(), cmap=lightdeep, extend='max', zorder=5)
#cc = plt.contour(lon, lat, -Z, [100, 500, 1000, 2000, 3000, 4000, 5000], colors='silver', linewidths=.5, transform=ccrs.PlateCarree(), zorder=10)
#plt.clabel(cc, inline=True, fontsize=7, fmt='%i')

# adjust the colorbar
#levels = np.linspace(vmin, vmax, num_levels)
#midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
#vals = np.interp(midp, [vmin, midpoint, vmax], [0, 0.5, 1])
#colors = colors(vals)
#colors=np.concatenate([[colors[0,:]], colors, [colors[-1,:]]],0)
#cmap, norm = from_levels_and_colors(levels, colors, extend=extent)

# plot WM classes
wm_class = np.sort(df_wm[class_type].unique())
wm_class_text = wm_class+1
co = np.full(len(wm_class),'nan')
class_idx = np.array([]).astype(int)
for i, iclass in enumerate(wm_class):
    lat_tmp = wm_lat[df_wm[class_type]==iclass].values
    lon_tmp = wm_lon[df_wm[class_type]==iclass].values
    if len(lon_tmp)>0:
        plt.scatter(lon_tmp, lat_tmp, 10, transform=ccrs.PlateCarree(), color=colors[iclass], alpha=.5, zorder=100)
        class_idx = np.append(class_idx, i)

plt.legend(wm_class_text[class_idx].astype(str)).set_zorder(102)



# Plot contours
#x,y = m(*np.meshgrid(lon,lat))
#c = plt.contourf(x, y, -Z, v, cmap=lightdeep, extend='max', zorder=1)
#cc = plt.contour(x, y, -Z, [100, 500, 1000, 2000, 3000, 4000], colors='lightgrey', linewidths=.5, zorder=1)
#ccc = plt.contour(x, y, -Z, [1000], colors='darkgrey', linewidths=2, zorder=1)
#plt.clabel(cc, inline=1, fontsize=10, colors='lightgrey', fmt='%d')

#m.fillcontinents(color='tan', zorder=50);
#m.drawparallels(np.arange(10,70,5), labels=[1,0,0,0], fontsize=12, fontweight='bold');
#m.drawmeridians(np.arange(-80, 5, 5), labels=[0,0,0,1], fontsize=12, fontweight='bold');

# plot AZMP stations
x, y = stationLon[index_BI],stationLat[index_BI]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Beachy Island', transform=ccrs.PlateCarree(), horizontalalignment='left', verticalalignment='center', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_MB],stationLat[index_MB]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Makkovik', transform=ccrs.PlateCarree(), horizontalalignment='left', verticalalignment='center', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_SI],stationLat[index_SI]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Seal Island', transform=ccrs.PlateCarree(), horizontalalignment='left', verticalalignment='center', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_WB],stationLat[index_WB]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' White Bay', transform=ccrs.PlateCarree(), horizontalalignment='left', verticalalignment='top', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_BB],stationLat[index_BB]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' Bonavista', transform=ccrs.PlateCarree(), horizontalalignment='left', verticalalignment='center', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_FC],stationLat[index_FC]
#m.scatter(x,y,10,marker='o',color='r', zorder=100)
plt.text(x[-6], y[-1], ' Flemish Cap', transform=ccrs.PlateCarree(), horizontalalignment='center', verticalalignment='top', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_SEGB],stationLat[index_SEGB]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' SE Grand Bank', transform=ccrs.PlateCarree(), horizontalalignment='left', verticalalignment='center', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_SESPB],stationLat[index_SESPB]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' SE St. Pierre Bank', transform=ccrs.PlateCarree(), horizontalalignment='right', verticalalignment='top', fontsize=12, color='r', fontweight='bold', zorder=100)

x, y = stationLon[index_SWSPB],stationLat[index_SWSPB]
#m.scatter(x,y,10,marker='o',color='r', zorder=10)
plt.text(x[-1], y[-1], ' SW St. Pierre Bank', transform=ccrs.PlateCarree(), horizontalalignment='right', verticalalignment='top', fontsize=12, color='r', fontweight='bold', zorder=100)

# Add high-res stations
#x, y = stationLon[index_S27],stationLat[index_S27]
#m.scatter(x[0],y[0],100,marker='*',color='k', zorder=20)
#plt.text(x[0], y[0]+50000, 'S27', horizontalalignment='right', verticalalignment='top', fontsize=12, color='r', fontweight='bold', zorder=100)

   
#### ---- Save Figure ---- ####
fig.set_size_inches(w=10, h=12)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
#plt.show()

