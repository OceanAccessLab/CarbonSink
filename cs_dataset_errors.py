'''
Script to estimate the sensitivity of parametrization during CO2sys application.
This is for Table 2 of AZMP carbonate dataset paper.

To run in /home/cyrf0006/AZMP/oa/Dataset_paper

Files need for this script were generated manually by editing cc_merge_azmp.py and copy AZMP_carbon_data.csv to the files below.

Frederic.Cyr@dfo-mpo.gc.ca - October 2022
'''

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PyCO2SYS import CO2SYS 


## Different files to test (generated manually...):
my_file0 = '/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv'
my_file1 = 'AZMP_carbon_data_K1K2_4.csv'
my_file2 = 'AZMP_carbon_data_K1K2_4_noNuts.csv'
my_file3 = 'AZMP_carbon_data_K1K2_4_TA-2sigma.csv'
my_file4 = 'AZMP_carbon_data_K1K2_4_TA+2sigma.csv'
my_file5 = 'AZMP_carbon_data_K1K2_9.csv'

## Load and arange
df1 = pd.read_csv(my_file1)
df2 = pd.read_csv(my_file2)
df3 = pd.read_csv(my_file3)
df4 = pd.read_csv(my_file4)
df5 = pd.read_csv(my_file5)
df1.set_index('Timestamp', inplace=True)
df2.set_index('Timestamp', inplace=True)
df3.set_index('Timestamp', inplace=True)
df4.set_index('Timestamp', inplace=True)
df5.set_index('Timestamp', inplace=True)
df1.index = pd.to_datetime(df1.index)
df2.index = pd.to_datetime(df2.index)
df3.index = pd.to_datetime(df3.index)
df4.index = pd.to_datetime(df4.index)
df5.index = pd.to_datetime(df5.index)
# Save copy for b)
df_tmp1 = df1.copy()
df_tmp2 = df2.copy()
# find indices for Gulf 2015
idx_TA = df1.loc[((df1.index.year==2015) & (df1.index.month==6)) & (df1['Region'].str.contains('GSL')), 'Total_Alkalinity_(umol/kg)'].index
# Keep only relevant variables
df1 = df1[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]
df2 = df2[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]
df3 = df3[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]
df4 = df4[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]
df5 = df5[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]


## a) Test dissociation constant choices

# Check % of S<20
print('Percentage of S<20:')
tmp = df1[df1['Salinity_(psu)']<20].shape[0] / df1.shape[0] * 100
print(str(tmp))

# Check error associated with choice of constant
tmp1 = df1[df1['Salinity_(psu)']<20]
tmp5 = df5[df5['Salinity_(psu)']<20]
diff = (tmp1-tmp5)/tmp1*100
diff.mean()
del diff, tmp1, tmp5

tmp1 = df1[df1['Salinity_(psu)']>=20]
tmp5 = df5[df5['Salinity_(psu)']>=20]
diff = (tmp1-tmp5)/tmp1*100
diff.mean()
del diff, tmp1, tmp5


## b) Test the effect of using nutrients or not
# Drop index with 3 parameters
df_TA_TIC_pH = df_tmp1.dropna(subset=['Total_Alkalinity_Measured_(umol/kg)', 'Inorganic_Carbon_Measured_(umol/kg)', 'pH_lab_(seawater_scale)'], thresh=3, axis=0)
df_tmp1.drop(df_TA_TIC_pH.index, inplace=True)
df_tmp2.drop(df_TA_TIC_pH.index, inplace=True)
# Isolate index with TA/TIC
df_TA_TIC1 = df_tmp1.dropna(subset=['Total_Alkalinity_Measured_(umol/kg)', 'Inorganic_Carbon_Measured_(umol/kg)'], thresh=2, axis=0)
df_TA_TIC2 = df_tmp2.dropna(subset=['Total_Alkalinity_Measured_(umol/kg)', 'Inorganic_Carbon_Measured_(umol/kg)'], thresh=2, axis=0)
# Isolate index with TA/pH
df_TA_pH1 = df_tmp1.dropna(subset=['Inorganic_Carbon_Measured_(umol/kg)', 'pH_lab_(seawater_scale)'], thresh=2, axis=0)
df_TA_pH2 = df_tmp2.dropna(subset=['Inorganic_Carbon_Measured_(umol/kg)', 'pH_lab_(seawater_scale)'], thresh=2, axis=0)
# Keep only some variables
df_TA_TIC1 = df_TA_TIC1[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]
df_TA_TIC2 = df_TA_TIC2[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]
df_TA_pH1 = df_TA_pH1[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]
df_TA_pH2 = df_TA_pH2[['Salinity_(psu)', 'Total_Alkalinity_(umol/kg)', 'Inorganic_Carbon_(umol/kg)', 'pH_Total_(total_scale)', 'Omega_Aragonite_(unitless)', 'Omega_Calcite_(unitless)', 'pCO2_(uatm)']]

# Differences
diff = (df_TA_TIC1-df_TA_TIC2)/df_TA_TIC1*100
diff.mean()
del diff
diff = (df_TA_pH1- df_TA_pH2)/df_TA_pH1*100
diff.mean()
del diff


## c) Effect of replacing missing TA
diff = (df1.loc[idx_TA]-df3.loc[idx_TA])/df1.loc[idx_TA]*100
diff.mean()
del diff
diff = (df1.loc[idx_TA]-df4.loc[idx_TA])/df1.loc[idx_TA]*100
diff.mean()


## d) Effect of calculating Omega, pCO2, etc. with TA-DIC or TA-pH couples (only possible where we have all 3 params)

# Isolate index with 3 parameters
df_tmp = pd.read_csv(my_file0)
df_TA_TIC_pH = df_tmp.dropna(subset=['Total_Alkalinity_Measured_(umol/kg)', 'Inorganic_Carbon_Measured_(umol/kg)', 'pH_lab_(total_scale)'], thresh=3, axis=0)
# Select Gulf only (N=2224)
df_TA_TIC_pH = df_TA_TIC_pH[df_TA_TIC_pH.Region=='GSL']

# Only consider S>20
df_TA_TIC_pH = df_TA_TIC_pH[df_TA_TIC_pH['Salinity_(psu)']>20]

# apply CO2sys (TA/TIC)
CO2dict = CO2SYS(df_TA_TIC_pH['Total_Alkalinity_Measured_(umol/kg)'], df_TA_TIC_pH['Inorganic_Carbon_Measured_(umol/kg)'], 1, 2, df_TA_TIC_pH['Salinity_(psu)'], df_TA_TIC_pH['pH_lab_temp_(degC)'], df_TA_TIC_pH['Temperature_(degC)'], 0, df_TA_TIC_pH['Depth_(dbar)'], df_TA_TIC_pH['Silicate_Concentration_(mmol/m3)'], df_TA_TIC_pH['Phosphate_Concentration_(mmol/m3)'], 1, 4, 1)
co2sys_TA_TIC = pd.DataFrame.from_dict(CO2dict)
co2sys_TA_TIC.index = df_TA_TIC_pH.index
del CO2dict
# apply CO2sys (TA/pH)
CO2dict = CO2SYS(df_TA_TIC_pH['Total_Alkalinity_Measured_(umol/kg)'], df_TA_TIC_pH['pH_lab_(total_scale)'], 1, 3, df_TA_TIC_pH['Salinity_(psu)'], df_TA_TIC_pH['pH_lab_temp_(degC)'], df_TA_TIC_pH['Temperature_(degC)'], 0, df_TA_TIC_pH['Depth_(dbar)'], df_TA_TIC_pH['Silicate_Concentration_(mmol/m3)'], df_TA_TIC_pH['Phosphate_Concentration_(mmol/m3)'], 1, 4, 1)
co2sys_TA_pH = pd.DataFrame.from_dict(CO2dict)
co2sys_TA_pH.index = df_TA_TIC_pH.index
del CO2dict
# apply CO2sys (TIC/pH)
CO2dict = CO2SYS(df_TA_TIC_pH['Inorganic_Carbon_Measured_(umol/kg)'], df_TA_TIC_pH['pH_lab_(total_scale)'], 2, 3, df_TA_TIC_pH['Salinity_(psu)'], df_TA_TIC_pH['pH_lab_temp_(degC)'], df_TA_TIC_pH['Temperature_(degC)'], 0, df_TA_TIC_pH['Depth_(dbar)'], df_TA_TIC_pH['Silicate_Concentration_(mmol/m3)'], df_TA_TIC_pH['Phosphate_Concentration_(mmol/m3)'], 1, 4, 1)
co2sys_TIC_pH = pd.DataFrame.from_dict(CO2dict)
co2sys_TIC_pH.index = df_TA_TIC_pH.index
del CO2dict


# Calculate differences
diff1 = (co2sys_TA_TIC-co2sys_TA_pH)/co2sys_TA_TIC*100
diff2 = (co2sys_TA_TIC-co2sys_TIC_pH)/co2sys_TA_TIC*100
diff3 = (co2sys_TIC_pH-co2sys_TA_pH)/co2sys_TA_TIC*100

diff1[['pHoutTOTAL', 'TAlk', 'TCO2', 'OmegaCAout', 'OmegaARout', 'pCO2out']].mean()
diff2[['pHoutTOTAL', 'TAlk', 'TCO2', 'OmegaCAout', 'OmegaARout', 'pCO2out']].mean()
diff3[['pHoutTOTAL', 'TAlk', 'TCO2', 'OmegaCAout', 'OmegaARout', 'pCO2out']].mean()



np.sort(diff1['pCO2out'])[round(diff1.shape[0]*.975)]
np.sort(diff1['pCO2out'])[round(diff1.shape[0]*.025)]
np.sort(diff1['pCO2out'])[round(diff2.shape[0]*.975)]
np.sort(diff1['pCO2out'])[round(diff2.shape[0]*.025)]

np.sort(diff1['OmegaARout'])[round(diff1.shape[0]*.975)]
np.sort(diff1['OmegaARout'])[round(diff1.shape[0]*.025)]
np.sort(diff1['OmegaARout'])[round(diff2.shape[0]*.975)]
np.sort(diff1['OmegaARout'])[round(diff2.shape[0]*.025)]
