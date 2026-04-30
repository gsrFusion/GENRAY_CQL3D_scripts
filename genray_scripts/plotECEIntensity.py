"""
Plots the ray traces and the RF power deposition density
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
import matplotlib

import os, sys
from matplotlib.path import Path
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

genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('figure', titlesize = 22)
plt.rc('legend', fontsize = 14)

#plots either the toroidal and/or poloidal ray trajectories
def plotECE():
    specIntensity = genray_ece_nc.variables['w_specific_intensity_nc'][:][:]
    freqs = genray_ece_nc.variables['wfreq_nc'][:]
    wr_em_nc = genray_ece_nc['wr_em_nc'][:]

    fig,ax = plt.subplots()

    for i in range(len(freqs)):
        ax.plot(wr_em_nc[i], specIntensity[i])

    ax.set_ylabel(r'Specific Intensity')
    ax.set_xlabel(r'major radius')
    ax.set_xlim([1.4,2.35])
    ax.set_ylim(bottom = 0)
    #ax.set_xlim([80,115])
    ax.set_title(shotNum)
    ax.legend(loc = 'best')
    fig.tight_layout()
    plt.show()

def main():
    plotECE()

main()
