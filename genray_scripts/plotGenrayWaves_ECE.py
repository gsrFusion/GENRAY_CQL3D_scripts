###
# Plot the ray trajectories and damping as predicted by GENRAY
# since it's a linear code, the damping is very wrong for LH, but the ray trajectory can be useful
###


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.collections import LineCollection
from scipy.interpolate import interp2d, interp1d, RectBivariateSpline

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
    gfileDict = getGfileDict.getGfileDict()

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    
    rgrid = gfileDict["rgrid"] #R points of the LCFS
    zgrid = gfileDict["zgrid"]
    B_zGrid = gfileDict["bzrz"]
    B_TGrid = gfileDict["btrz"]
    B_rGrid = gfileDict["brrz"]

    wr  = genray_ece_nc.variables['wr_em_nc'][:] #major radius of the ray at each point along the trace, in m
    freqs = genray_ece_nc.variables['wfreq_nc'][:]
    print(f'freqs: {freqs}')
    nparas = genray_ece_nc.variables['wnpar'][:]

    wz  = genray_ece_nc.variables['wz_em_nc'][:] #height of the ray at each point along the trace, in m
    
    B_totstrength = np.sqrt(np.square(B_zGrid) + np.square(B_rGrid) + np.square(B_TGrid))
    f2_ce = 27.9906*B_totstrength*2

    fig,ax = plt.subplots(figsize = (4.25,7))
    #plt.subplots_adjust(left=0.22,bottom = .1)
    ax.set_ylabel("Z (m)")
    ax.set_xlabel("R (m)")

    #ax.text(1.98, 1.28, r'$2\Omega_e = $77 GHz', rotation = 0, fontsize = 16, color = 'blue', ha='center', va='center')
    ax.text(2.3, .99, r'$2\Omega_e = $91 GHz', rotation = 90, fontsize = 17, color = 'blue', ha='center', va='center')

    ax.contour(rgrid, zgrid, f2_ce, levels = [94], colors= ['blue'], linewidths=3, linestyles = 'dotted', zorder = 10)


    print(nparas.shape)
    for ray in range(len(wr)):
        #if (freqs[ray] < 90 or freqs[ray] > 94) :
        #    continue
        #if not (int(freqs[ray]) in [91,]) :
        #    continue
        print(ray)
        #delpwr[ray,:] = delpwr[ray,:]/delpwr[ray,0] #normalize the ray power to that ray's starting power
        ax.plot(wr[ray][wr[ray] > 1], wz[ray][wr[ray] > 1],lw = 2, color = 'limegreen')

    ax.plot(xlim, ylim, color = 'grey', lw = 3,zorder = 6)#plot wall
    ax.scatter([2.5],[0.004], color ='r',marker='*',zorder = 10, s=150)
    ax.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS
    ax.set_ylim(min(ylim)-.05, max(ylim)+.05)
    ax.set_xlim(min(xlim)-.05, max(xlim)+.05)
    ax.set_aspect('equal')
    
    gfileDict = getGfileDict.getGfileDict(targetDir = '/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_203619.04130Temp')
    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    ax.plot(xlim, ylim, color = 'grey', lw = 3, linestyle = 'dashed',zorder = 5)#plot wall

    helper.drawFluxSurfaces(ax)
    fig.tight_layout()
    #plt.savefig('203912_ECErays_first.jpeg',dpi=300)



def plotToroidalTrace():
    wr  = genray_ece_nc.variables["wr_em_nc"][:] #major radius of the ray at each point along the trace, in m
    toroidalAngle = genray_ece_nc.variables["wphi_em_nc"][:] 

    norm = plt.Normalize(0, 1)

    fig,ax = plt.subplots(figsize = (6.5,5.5))

    ax.set_ylabel("Y (m)")
    ax.set_xlabel("X (m)")
    thetas = np.linspace(0,2*np.pi,100)
    ax.plot(np.cos(thetas), np.sin(thetas), color = 'k', linewidth = 2)
    ax.plot(np.cos(thetas)*2.355, np.sin(thetas)*2.355, color = 'k', linewidth = 2)
    ax.set_aspect('equal')
    for ray in range(len(wr)):
        rayX = (np.cos(toroidalAngle[ray][wr[ray] > 1])*wr[ray][wr[ray] > 1])
        print(rayX)
        rayY = (np.sin(toroidalAngle[ray][wr[ray] > 1])*wr[ray][wr[ray] > 1])
        
        ax.plot(rayX, rayY)
    

def main():
    plotRays()
    #plotToroidalTrace()
    plt.show()

main()
