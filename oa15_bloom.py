import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gsw
import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from lmfit.models import LorentzianModel
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy import optimize
from datetime import timedelta


# Read Lab Sea bloom data
bloomt = pd.read_csv('/home/cyrf0006/research/carbonSink/OA15/CentralLabSea_SpringBloom.csv')
bloomt.index = bloomt.Year


stam = bloomt['t[start]'].mean()
stas = bloomt['t[start]'].std()/2
durm = bloomt['t[duration]'].mean()
durs = bloomt['t[duration]'].std()/2
ampm = bloomt['Amplitude[fit]'].mean()
amps = bloomt['Amplitude[fit]'].std()/2
magm = bloomt['Magnitude[fit]'].mean()
mags = bloomt['Magnitude[fit]'].std()/2

## ---- plot index ---- ##
width=.5
fig = plt.figure()
fig.clf()
ax = plt.subplot2grid((4, 1), (0, 0))
bloomt['t[start]'].plot(kind='bar', width = width, zorder=10)
plt.bar(12, bloomt['t[start]'].loc[2015], width = width, zorder=20, color='indianred')
ticks = ax.xaxis.get_ticklocs()
plt.plot([ticks[0]-1, ticks[-1]+1], [stam, stam], '--k')
plt.fill_between([ticks[0]-1, ticks[-1]+1], [stam-stas, stam-stas], [stam+stas, stam+stas], facecolor='gray', alpha=.2)
plt.ylabel('[DOY]', fontsize=14)
plt.xlabel(' ')
plt.grid()
#plt.ylim([-2.5,2.5])
plt.title('Bloom initiation', weight='bold', fontsize=14, loc='left')
ax.tick_params(labelbottom=False)

ax2 = plt.subplot2grid((4, 1), (1, 0))
bloomt['t[duration]'].plot(kind='bar', width = width, zorder=10)
plt.bar(12, bloomt['t[duration]'].loc[2015], width = width, zorder=20, color='indianred')
ticks = ax2.xaxis.get_ticklocs()
plt.plot([ticks[0]-1, ticks[-1]+1], [durm, durm], '--k')
plt.fill_between([ticks[0]-1, ticks[-1]+1], [durm-durs, durm-durs], [durm+durs, durm+durs], facecolor='gray', alpha=.2)
plt.ylabel('[days]', fontsize=14)
plt.xlabel(' ')
plt.grid()
#plt.ylim([-2.5,2.5])
plt.title('Bloom duration', weight='bold', fontsize=14, loc='left')
ax2.tick_params(labelbottom=False)
#ax2.tick_params(axis='x', rotation=0)

ax3 = plt.subplot2grid((4, 1), (2, 0))
bloomt['Amplitude[fit]'].plot(kind='bar', width = width, zorder=10)
plt.bar(12, bloomt['Amplitude[fit]'].loc[2015], width = width, zorder=20, color='indianred')
ticks = ax3.xaxis.get_ticklocs()
plt.plot([ticks[0]-1, ticks[-1]+1], [ampm, ampm], '--k')
plt.fill_between([ticks[0]-1, ticks[-1]+1], [ampm-amps, ampm-amps], [ampm+amps, ampm+amps], facecolor='gray', alpha=.2)
plt.ylabel(r'[$\rm mg\,m^{-3}$]', fontsize=14)
plt.xlabel(' ')
plt.grid()
#plt.ylim([-2.5,2.5])
plt.title('Bloom Amplitude', weight='bold', fontsize=14, loc='left')
ax3.tick_params(labelbottom=False)

ax4 = plt.subplot2grid((4, 1), (3, 0))
bloomt['Magnitude[fit]'].plot(kind='bar', width = width, zorder=10)
plt.bar(12, bloomt['Magnitude[fit]'].loc[2015], width = width, zorder=20, color='indianred')
ticks = ax4.xaxis.get_ticklocs()
plt.plot([ticks[0]-1, ticks[-1]+1], [magm, magm], '--k')
plt.fill_between([ticks[0]-1, ticks[-1]+1], [magm-mags, magm-mags], [magm+mags, magm+mags], facecolor='gray', alpha=.2)
plt.ylabel(r'[$\rm mg\,m^{-3}\,day$]', fontsize=14)
plt.xlabel(' ')
plt.grid()
#plt.ylim([-2.5,2.5])
plt.title('Bloom Magnitude', weight='bold', fontsize=14, loc='left')
ax4.tick_params(axis='x', rotation=90)


fig_name = 'Bloom_CLS.png'
fig.set_size_inches(w=9,h=15)
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)
