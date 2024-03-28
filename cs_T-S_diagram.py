import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import datetime
import gsw

df = pd.read_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/AZMP-NL_carbon.csv', delimiter=',')

# set index
df.index = pd.to_datetime(df.Timestamp)
df.drop(columns='Timestamp', inplace=True)
df.drop(columns='Region', inplace=True)

cmap = plt.cm.get_cmap('viridis', 6)
cmap = plt.cm.get_cmap('Paired', 6)
#cmap = plt.cm.Paired(5)
plt.scatter(df['Salinity_(psu)'], df['Temperature_(degC)'], c=df.index.year, cmap=cmap, alpha=.7)
plt.colorbar()


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
