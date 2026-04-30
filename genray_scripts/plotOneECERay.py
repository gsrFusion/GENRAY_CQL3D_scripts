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
    
    wr  = genray_ece_nc.variables['wr'][:]/100 #major radius of the ray at each point along the trace, in m
    nparas = genray_ece_nc.variables['wnpar'][:]

    wz  = genray_ece_nc.variables['wz'][:]/100 #height of the ray at each point along the trace, in m
    delpwr= genray_ece_nc.variables["delpwr"][:] #power in the ray at each point

    norm = plt.Normalize(0, 1)

    fig,ax = plt.subplots(figsize = (4.25,7.1))
    plt.subplots_adjust(left=0.22,bottom = .1)
    ax.set_ylabel("Z (m)")
    ax.set_xlabel("R (m)")

    print(nparas.shape)
    for ray in range(len(wr)):
        #delpwr[ray,:] = delpwr[ray,:]/delpwr[ray,0] #normalize the ray power to that ray's starting power
        ax.plot(wr[ray], wz[ray],lw = 3)

    ax.plot(xlim, ylim, 'r', lw = 2)#plot wall
    ax.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS

    ax.set_aspect('equal')
    
    if machine == 'FENIX':
        ax.set_ylim(-.35,.35)
        ax.set_xlim(.45, 1)
    else:
        pass
        #ax.set_ylim(min(ylim)*1.05, max(ylim)*1.05)
        #ax.set_xlim(min(xlim)*.95, max(xlim)*1.05)

    helper.drawFluxSurfaces(ax)

def plotToroidalTrace():
    wr  = genray_ece_nc.variables["wr"][:]/100 #major radius of the ray at each point along the trace, in m
    toroidalAngle = genray_ece_nc.variables["wphi"][:] 

    delpwr= genray_ece_nc.variables["delpwr"][:] #power in the ray at each point

    norm = plt.Normalize(0, 1)

    fig,ax = plt.subplots(figsize = (6.5,5.5))

    ax.set_ylabel("Y (m)")
    ax.set_xlabel("X (m)")
    thetas = np.linspace(0,2*np.pi,100)
    ax.plot(np.cos(thetas), np.sin(thetas), color = 'k', linewidth = 2)
    ax.plot(np.cos(thetas)*2.355, np.sin(thetas)*2.355, color = 'k', linewidth = 2)
    ax.set_aspect('equal')
    for ray in range(len(wr)):
        rayX = (np.cos(toroidalAngle[ray])*wr[ray])
        rayY = (np.sin(toroidalAngle[ray])*wr[ray])
        points = np.array([rayX, rayY]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap=plt.cm.turbo)
        # Set the values used for colormapping
        lc.set_array(delpwr[ray]/delpwr[ray][0])
        lc.set_linewidth(2)
        ax.add_collection(lc)

def main():
    plotRays()
    #plotToroidalTrace()
    plt.show()

main()
