import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import datetime
import gsw
import cmocean as cmo

df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv', delimiter=',')

# Select NL only
df = df[df.Region=='NL']

# set index
df.index = pd.to_datetime(df.Timestamp)
df.drop(columns='Timestamp', inplace=True)
df.drop(columns='Region', inplace=True)

# Calculate AOU
SA = gsw.SA_from_SP(df['Salinity_(psu)'], df['Depth_(dbar)'], df['Longitude_(degEast)'], df['Latitude_(degNorth)'])
CT = gsw.CT_from_t(SA, df['Temperature_(degC)'], df['Depth_(dbar)'])
O2sol = gsw.O2sol(SA, CT, df['Depth_(dbar)'],  df['Longitude_(degEast)'], df['Latitude_(degNorth)']) # in umol/kg
O2sol = O2sol/43.570 # in ml/l 
df['AOU'] = O2sol - df['Dissolved_Oxygen_(mL/L)'] # calculate Apparent Oxygen Utilisation

# Select 20515 only
df15 = df[df.index.year==2015]

cmap = plt.cm.get_cmap('viridis', 9)
cmap = plt.cm.get_cmap('Paired', 9)
#cmap = plt.cm.Paired(5)
plt.scatter(df['Salinity_(psu)'], df['Temperature_(degC)'], c=df.index.year, cmap=cmap, alpha=.7)
plt.colorbar()
plt.scatter(df15['Salinity_(psu)'], df15['Temperature_(degC)'], color='r')

salt = df['Salinity_(psu)']
temp = df['Temperature_(degC)']
# Figure out boudaries (mins and maxs)
smin = salt.min() - (0.01 * salt.min())
smax = salt.max() + (0.01 * salt.max())
tmin = temp.min() - (0.01 * temp.max())
tmax = temp.max() + (0.01 * temp.max())
 
# Calculate how many gridcells we need in the x and y dimensions
xdim = int(round((smax-smin)/0.05+1,0))
ydim = int(round((tmax-tmin)/0.05+1,0))

# Create empty grid of zeros
dens = np.full((ydim,xdim), 0)
 
# Create temp and salt vectors of appropiate dimensions
ti = np.linspace(1,ydim-1,ydim)*0.05+tmin
si = np.linspace(1,xdim-1,xdim)*0.05+smin
 
# Loop to fill in grid with densities
for j in range(0,int(ydim)):
    for i in range(0, int(xdim)):
        dens[j,i]=gsw.rho(si[i],ti[j],0)
 
# Substract 1000 to convert to sigma-t
dens = dens - 1000
 
# Plot data ***********************************************
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
CS = plt.contour(si,ti,dens, linestyles=':', colors='silver')
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.0f') # Label every second level
 
c = ax1.scatter(salt, temp, c=df.index.year, cmap=cmap, alpha=.7)

ax1.set_xlabel('Salinity')
ax1.set_ylabel(r'Temperature ($^{\circ}$C)')
cb = plt.colorbar(c)



### Playing around ###
# TA vs DIC
plt.scatter(df['Total_Alkalinity_(umol/kg)'], df['Inorganic_Carbon_(umol/kg)'], c=df.index.year, cmap=cmap, alpha=.7)
plt.colorbar()
plt.show()

# S vs DO
plt.scatter(df['Salinity_(psu)'], df['Dissolved_Oxygen_(mL/L)'], c=df.index.year, cmap=cmap, alpha=.7)
plt.colorbar()
plt.show()

# DIC vs AOU
plt.scatter(df['AOU'], df['Inorganic_Carbon_(umol/kg)'], c=df.index.year, cmap=cmap, alpha=.7)
plt.colorbar()
plt.show()
