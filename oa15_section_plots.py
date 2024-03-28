# -*- coding: utf-8 -*-
"""
SCript to make multiple section plots
Same as in AZMP CC SAR figure

variable options are:
VARIABLE = 'Temperature_(degC)'
VARIABLE = 'Salinity_(psu)'
VARIABLE = 'Oxygen_Saturation_(%)'
VARIABLE = 'Total_Alkalinity_(umol/kg)'
VARIABLE = 'Inorganic_Carbon_(umol/kg)'
VARIABLE = 'pCO2_(uatm)',
VARIABLE = 'pH_Total_(total_scale)'

Frederic.Cyr@dfo-mpo.gc.ca
Sept. 2023
"""


import cc_tools as cc

years = [2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022]
seasons = ['spring','summer','fall']
#seasons = ['spring']
variables = ['Inorganic_Carbon_(umol/kg)', 'Total_Alkalinity_(umol/kg)', 'Oxygen_Saturation_(%)', 'AOU']
sections = ['SI','BB','FC']


for year in years:
    print(str(year))
    for season in seasons:
        print(' ' + season)
        for variable in variables:
            print('  ' + variable )
            for section in sections:
                print('   ' + section)
                
                cc.section_plot(variable, year, season, section, 1000)
            

