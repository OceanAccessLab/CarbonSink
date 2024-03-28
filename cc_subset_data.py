import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


df = pd.read_csv('/home/cyrf0006/github/AZMP-NL/datasets/carbonates/AZMP_carbon_data.csv', delimiter=',')


# Only NL data
df = df[df.Region=='NL']

# Add date columns (was for Cassie)
## df['Year'] = df.Timestamp.dt.year  
## df['Month'] = df.Timestamp.dt.month                            
## df['Day'] = df.Timestamp.dt.day

# Save NL
df.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/AZMP-NL_carbon.csv', float_format='%.4f', index=False)         

# Save per section
df_BI = df.loc[df.Station_Name.str.contains('MB-')]
df_MB = df.loc[df.Station_Name.str.contains('BI-')]
df_SI = df.loc[df.Station_Name.str.contains('SI-')]
df_BB = df.loc[df.Station_Name.str.contains('BB-')]
df_FC = df.loc[df.Station_Name.str.contains('FC-')]
df_SEGB = df.loc[df.Station_Name.str.contains('SEGB-')]
df_SESPB = df.loc[df.Station_Name.str.contains('SESPB-')]
df_SWSPB = df.loc[df.Station_Name.str.contains('SWSPB-')]
df_BI.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/BI_carbon.csv', float_format='%.4f', index=False)         
df_MB.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/MB_carbon.csv', float_format='%.4f', index=False)
df_SI.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/SI_carbon.csv', float_format='%.4f', index=False)         
df_BB.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/BB_carbon.csv', float_format='%.4f', index=False)         
df_FC.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/FC_carbon.csv', float_format='%.4f', index=False)         
df_SEGB.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/SEGB_carbon.csv', float_format='%.4f', index=False)         
df_SESPB.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/SESPB_carbon.csv', float_format='%.4f', index=False)         
df_SWSPB.to_csv('/home/cyrf0006/github/oceanaccess/CarbonSink/data/SWSPB_carbon.csv', float_format='%.4f', index=False)         


# 
