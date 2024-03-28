'''
Script to estimate error with the different carbonate equilibria pairs 
Comparision exercise with M. Starr.

Frederic.Cyr@dfo-mpo.gc.ca - July 31 2023
'''

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PyCO2SYS import CO2SYS 


## Load file (manually added "x" in important columns
my_file0 = '/home/cyrf0006/AZMP/oa/Dataset_paper/Cyr_Comparaison.xlsx'



## d) Effect of calculating Omega, pCO2, etc. with TA-DIC or TA-pH couples (only possible where we have all 3 params)

# Isolate index with 3 parameters
df_tmp = pd.read_excel(my_file0, skiprows=[0,2])
df_TA_TIC_pH = df_tmp.dropna(subset=['xAt', 'xCIDT', 'xpH  labo'], thresh=3, axis=0)


### *** Need to check which scale to use for pH ***
#HERE!!!
# apply CO2sys (TA/TIC)
CO2dict = CO2SYS(df_TA_TIC_pH['xAt'], df_TA_TIC_pH['xCIDT'], 1, 2, df_TA_TIC_pH['PSAL'], df_TA_TIC_pH['LABT_01'], df_TA_TIC_pH['TE90'], 0, df_TA_TIC_pH['PRES'], df_TA_TIC_pH['xSi_03'], df_TA_TIC_pH['xPO4_03'], 1, 4, 1)
co2sys_TA_TIC = pd.DataFrame.from_dict(CO2dict)
co2sys_TA_TIC.index = df_TA_TIC_pH.index
del CO2dict
# apply CO2sys (TA/pH)
CO2dict = CO2SYS(df_TA_TIC_pH['xAt'], df_TA_TIC_pH['xpH  labo'], 1, 3, df_TA_TIC_pH['PSAL'], df_TA_TIC_pH['LABT_01'], df_TA_TIC_pH['TE90'], 0, df_TA_TIC_pH['PRES'], df_TA_TIC_pH['xSi_03'], df_TA_TIC_pH['xPO4_03'], 1, 4, 1)
co2sys_TA_pH = pd.DataFrame.from_dict(CO2dict)
co2sys_TA_pH.index = df_TA_TIC_pH.index
del CO2dict
# apply CO2sys (TIC/pH)
CO2dict = CO2SYS(df_TA_TIC_pH['xCIDT'], df_TA_TIC_pH['xpH  labo'], 2, 3, df_TA_TIC_pH['PSAL'], df_TA_TIC_pH['LABT_01'], df_TA_TIC_pH['TE90'], 0, df_TA_TIC_pH['PRES'], df_TA_TIC_pH['xSi_03'], df_TA_TIC_pH['xPO4_03'], 1, 4, 1)
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


#* Note on CO2sys usage:
#CO2dict = CO2SYS(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, TEMPIN, TEMPOUT, PRESIN, PRESOUT,
#    SI, PO4, pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS, NH3=0.0, H2S=0.0, KFCONSTANT=1)
    
