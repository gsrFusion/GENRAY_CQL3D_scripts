###
# Plot the ray trajectories and damping as predicted by GENRAY
# since it's a linear code, the damping is very wrong for LH, but the ray trajectory can be useful
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
genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')#netCDF4.Dataset(f'{targetDir}/genray.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 14)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

#adds the ray traces to ax
def plotPowerDep():
    powden_e = genray_nc.variables["powden_e"][:]
    rho_bin_center = genray_nc.variables["rho_bin_center"][:]
    print(rho_bin_center)
    fig,ax = plt.subplots()
    plt.subplots_adjust(left=0.22,bottom = .1)
    ax.set_ylabel("power (erg/s cm^3)")
    ax.set_xlabel("rho_p")

    ax.plot(rho_bin_center, powden_e, lw = 2)
    fig.tight_layout()
    plt.show()


plotPowerDep()
