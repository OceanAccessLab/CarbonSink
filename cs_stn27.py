import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import datetime
from scipy import stats
import seawater as swx
from PyCO2SYS import CO2SYS 
from PyCO2SYS.meta import version
import gsw

# parameters
ztop = 15
zbot = 50
my_var = 'AOU'

if my_var == 'pH':
    variable = 'pH_tot'
elif my_var == 'NO3':
    variable = 'Nitrate_Concentration_(mmol/m3)'
elif my_var == 'PO4':
    variable = 'Phosphate_Concentration_(mmol/m3)'
elif my_var == 'SiO':
    variable = 'Silicate_Concentration_(mmol/m3)'
elif my_var == 'O2':
    variable = 'Dissolved_Oxygen_(mL/L)'
elif my_var == 'O2s':
    variable = 'Oxygen_Saturation_(%)'
elif my_var == 'omegaC':
    variable = 'Omega_Calcite_(--)'
elif my_var == 'pco2':
    variable = 'pCO2_(uatm)'
elif my_var == 'TA':
    variable = 'Total_Alkalinity_(umol/kg)'
elif my_var == 'DIC':
    variable = 'Inorganic_Carbon_(umol/kg)'
elif my_var == 'S':
    variable = 'Salinity_(psu)'
elif my_var == 'T':
    variable = 'Temperature_(degC)'
else:
    variable = my_var

# Load all data
df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv', delimiter=',')
df.set_index('Timestamp', inplace=True)
df.index = pd.to_datetime(df.index)   

# Calculate AOU
SA = gsw.SA_from_SP(df['Salinity_(psu)'], df['Depth_(dbar)'], df['Longitude_(degEast)'], df['Latitude_(degNorth)'])
CT = gsw.CT_from_t(SA, df['Temperature_(degC)'], df['Depth_(dbar)'])
O2sol = gsw.O2sol(SA, CT, df['Depth_(dbar)'],  df['Longitude_(degEast)'], df['Latitude_(degNorth)']) # in umol/kg
O2sol = O2sol/43.570 # in ml/l 
df['AOU'] = O2sol - df['Dissolved_Oxygen_(mL/L)'] # calculate Apparent Oxygen Utilisation

# keep only Stn 27 and Riki
df_s27 = df[df.Station_Name=='STN27']

# drop non-numeric columns
df_s27 = df_s27.drop(columns=['Region', 'Trip_Name', 'Station_Name', 'Latitude_(degNorth)', 'Longitude_(degEast)'])

# isolate surface (top 10m) and bottom 50m
df_s27_top = df_s27[df_s27['Depth_(dbar)']<=ztop].resample('M').mean()
df_s27_top_std = df_s27[df_s27['Depth_(dbar)']<=ztop].resample('M').std()
df_s27_bot = df_s27[df_s27['Depth_(dbar)']>=(df_s27['Depth_(dbar)'].max() - zbot)].resample('M').mean()
df_s27_bot_std = df_s27[df_s27['Depth_(dbar)']>=(df_s27['Depth_(dbar)'].max() - zbot)].resample('M').std()
s27_top = df_s27_top[variable].dropna()
s27_top_std = df_s27_top_std[variable].dropna()
s27_bot = df_s27_bot[variable].dropna()
s27_bot_std = df_s27_bot_std[variable].dropna()

plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(s27_top.index, s27_top.values, linewidth=2, marker='o')
ax.plot(s27_bot.index, s27_bot.values, linewidth=2, marker='o')
ax.legend(['top 15m', 'bottom 50m'])
ax.fill_between(s27_top.index, s27_top+s27_top_std/2, s27_top-s27_top_std/2, facecolor='gray', alpha=.2)
ax.fill_between(s27_bot.index, s27_bot+s27_bot_std/2, s27_bot-s27_bot_std/2, facecolor='gray', alpha=.2)
ax.grid()
plt.xlim([pd.to_datetime('2014-07-01'), pd.to_datetime('2020-03-01')])
YLIM = ax.get_ylim()
Y1 = [YLIM[0], YLIM[0]]
Y2 = [YLIM[1], YLIM[1]]
X = [pd.to_datetime('2015-11-01'), pd.to_datetime('2016-01-30')]
ax.fill_between(X, Y1, Y2, facecolor='red', alpha=.15)
plt.ylim(YLIM)
# ylabel
if my_var == 'pH':
    plt.ylabel(r'$pH_{is,T}$')
elif my_var == 'NO3':
    plt.ylabel(r'$[NO_3]~(mmol\,m^{-3})$')
elif my_var == 'PO4':
    plt.ylabel(r'$[PO_4]~(mmol\,m^{-3})$')
elif my_var == 'SiO':
    plt.ylabel(r'$[Si]~(mmol\,m^{-3})$')
elif my_var == 'O2s':
    plt.ylabel(r'$O_{2_{sat}}~(\%)$')
elif my_var == 'AOU':
    plt.ylabel(r'$AOU~(mL\,L^{-1})$')
elif my_var == 'O2':
    plt.ylabel(r'$[O_2]~(mL\,L^{-1})$')
elif my_var == 'omegaC':
    plt.ylabel(r'$\Omega_C$')
elif my_var == 'pco2':
    plt.ylabel(r'$pCO_2~(\mu atm)$')
elif my_var == 'TA':
    plt.ylabel(r'$TA~(\mu mol\,kg^{-1})$')    
elif my_var == 'DIC':
    plt.ylabel(r'$DIC~(\mu mol\,kg^{-1})$')
elif my_var == 'S':
    plt.ylabel(r'S')
elif my_var == 'T':
    plt.ylabel(r'T')
else:
    variable = ''

fig.set_size_inches(w=12,h=4)
fig_name2 = 'Stn27_' + my_var + '.png'
fig.savefig(fig_name2, dpi=200)
os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name2 + ' ' + fig_name2)

# montage in linux
# montage Stn27_TA.png Stn27_DIC.png Stn27_pH.png Stn27_omegaC.png Stn27_pco2.png -tile 1x5 -geometry +10+10  -background white  oa15_stn27_pH.png

#montage Stn27_T.png Stn27_S.png Stn27_O2.png  Stn27_NO3.png Stn27_PO4.png -tile 1x5 -geometry +10+10  -background white  oa15_stn27_nut.png

