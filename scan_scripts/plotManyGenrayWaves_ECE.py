###
# Plot the ray trajectories and damping as predicted by GENRAY
# since it's a linear code, the damping is very wrong for LH, but the ray trajectory can be useful
###


import numpy as np
import matplotlib.pyplot as plt
import os, sys

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getGfileDict
import helperFunctions as helper
import netCDF4


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
def plotRays(targetDir, label, ax, color):
    gfileDict = getGfileDict.getGfileDict(targetDir = targetDir)
    genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    
    wr  = genray_ece_nc.variables['wr_em_nc'][:] #major radius of the ray at each point along the trace, in m
    wz  = genray_ece_nc.variables['wz_em_nc'][:] #height of the ray at each point along the trace, in m

    for ray in range(len(wr)):
    
        ax.plot(wr[ray][wr[ray] > 1], wz[ray][wr[ray] > 1], lw = 3, color = color)

    ax.plot(xlim, ylim, 'r', lw = 2)#plot wall
    ax.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS

    ax.set_ylim(min(ylim)*1.05, max(ylim)*1.05)
    ax.set_xlim(min(xlim)*.95, max(xlim)*1.05)

    helper.drawFluxSurfaces(ax, targetDir=targetDir)

def main():
    fig,ax = plt.subplots(figsize = (4.25,7.1))
    plt.subplots_adjust(left=0.22,bottom = .1)
    ax.set_ylabel("Z (m)")
    ax.set_xlabel("R (m)")
    ax.set_aspect('equal')

    stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_203619.04130/DIIID_203619.04130_expSpectrum_2Zeff/DIIID_203619.04130_expSpectrum_2Zeff_'
    targetDirs = [f'{stem}first',
                  f'{stem}second',]
    colors = ['tab:blue','tab:orange']
    labels = ['first','second']
    for i in range(len(targetDirs)):
        plotRays(targetDirs[i], labels[i], ax,colors[i])
    
    ax.legend()
    fig.tight_layout()
    plt.show()

main()
