###
# Plots the Nperp of an ECE ray. May be useful for debugging
###


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.collections import LineCollection

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
import helperFunctions as helper
import getInputFileDictionary
import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
machine = getTargetInfo.getMachine()
genray_in = getInputFileDictionary.getInputFileDictionary('genray')
gfileDict = getGfileDict.getGfileDict()
genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')#netCDF4.Dataset(f'{targetDir}/genray.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 14)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

#returns the index of the array whose element is closest to value
def findNearestIndex(value, array):
    idx = (np.abs(array - value)).argmin()

    return idx

#adds the ray traces to ax
def plotRays():
    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    
    wr  = genray_ece_nc.variables['wr_em_nc'][:] #major radius of the ray at each point along the trace, in m
    wcnpar_em_nc = genray_ece_nc.variables['wcnpar_em_nc'][:]
    wcnper_em_nc = genray_ece_nc.variables['wcnper_em_nc'][:]

    norm = plt.Normalize(0, 1)

    fig,ax = plt.subplots()
    for i in range(len(wr)):
        if i > 0:
            continue
        ax.plot(wr[i][wr[i] > 1],wcnper_em_nc[i][wr[i] > 1],lw = 3)
        ax.scatter(wr[i][wr[i] > 1][:1],wcnper_em_nc[i][wr[i] > 1][:1],color = 'r',marker = '*',s=100)
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Nperp")

    fig.tight_layout()


def main():
    plotRays()
    plt.show()

main()


