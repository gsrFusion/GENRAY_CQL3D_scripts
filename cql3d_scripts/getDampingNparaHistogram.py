"""
Plots the ray traces and the RF power deposition density
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.collections import LineCollection
from scipy.interpolate import interp2d
import matplotlib
import os, sys
from scipy.signal import find_peaks
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import getGfileDict
gfileDict = getGfileDict.getGfileDict()
import helperFunctions as helper
import getInputFileDictionary
genray_in = getInputFileDictionary.getInputFileDictionary('genray')
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d')

import netCDF4
import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
print(targetDir)
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)

nparas = cqlrf_nc.variables['wnpar'][:]
delpwrs= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point

nparasBinEdges = np.linspace(np.nanmin(nparas)-.1, np.nanmax(nparas) + .1, 100)
binCenters = (nparasBinEdges[1:] + nparasBinEdges[:-1])/2
numPerBin = np.zeros(len(binCenters))
"""
for ray in range(len(nparas)):
    npara = nparas[ray]
    nanIndex = np.where(np.isnan(nparas[ray]))[0]
    if len(nanIndex) == 0 or nanIndex > 5:
        if len(nanIndex) == 0:
            nanIndex = len(npara)
        else:
            nanIndex = nanIndex[0]
        npara = npara[:nanIndex]
        delpwr = delpwrs[ray][:nanIndex]

        avgNpara = npara[:-1]#(npara[1:]+npara[:-1])/2
        powerDiff = -1*np.diff(delpwr)
        indices = np.digitize(avgNpara, nparasBinEdges, right = False)
        numPerBin[indices] += powerDiff

fig,ax = plt.subplots()
ax.plot(binCenters, (numPerBin/np.max(numPerBin))*(1/binCenters**2), lw = 2)
ax.set_ylabel(f'deposited power (A.U.)')
ax.set_xlabel(r'$N_{||}$')
"""
fig,ax = plt.subplots()
for ray in range(len(nparas)):
    if ray >= 28:
        continue
    ax.plot(nparas[ray], (delpwrs[ray]))
ax.set_xlim([0,5])
#ax.set_ylim([0,-3e9])
ax.set_ylabel(r'delpwr')
ax.set_xlabel(r'$N_{||}$')
fig.tight_layout()
plt.show()