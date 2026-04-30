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
import getInputFileDictionary
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
def plotRaysOfDirectory(ax, targetDir, maxDelPwrPlot):
    genray_in = getInputFileDictionary.getInputFileDictionary('genray',pathprefix=f'{parentdir}/')
    gfileDict = getGfileDict.getGfileDict(pathprefix=f'{parentdir}/')
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    
    wr  = genray_nc.variables["wr"][:] #major radius of the ray at each point along the trace

    wz  = genray_nc.variables["wz"][:] #height of the ray at each point along the trace
    delpwr= genray_nc.variables["delpwr"][:] #power in the ray at each point
    wr *= .01; wz*=.01 #convert to m from cm
    nparas = genray_nc.variables['wnpar'][:]
    

    norm = plt.Normalize(0, 1)

    #plt.subplots_adjust(left=0.22,bottom = .1)
    ax.set_ylabel("z (m)")
    ax.set_xlabel("R (m)")

    for ray in range(len(wr)):
        delpwr[ray,:] = delpwr[ray,:]/delpwr[ray,0] #normalize the ray power to that ray's starting power
        mostPowerDep = findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]) #find the index of the last ray point we want to plot
        #print(f' ray: {ray}. nparas[ray][0]: {nparas[ray][0]}')

        points = np.array([wr[ray][:mostPowerDep], wz[ray][:mostPowerDep]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap=plt.cm.jet)
        # Set the values used for colormapping
        lc.set_array(delpwr[ray][:mostPowerDep])
        lc.set_linewidth(1)
        ax.add_collection(lc)

    ax.plot(xlim, ylim, 'r', lw = 2)#plot wall
    ax.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS

    #ax.set_title(f"Plotting Rays until {(maxDelPwrPlot) * 100} %\n ray power deposition")
    ax.set_aspect('equal')
    
    ax.set_ylim(min(ylim)*1.05, max(ylim)*1.05)
    ax.set_xlim(min(xlim)*.95, max(xlim)*1.05)

    drawFluxSurfaces(ax, gfileDict)

#draw poloidal flux surfaces
def drawFluxSurfaces(ax, gfileDict):
    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    ## THIS NEEDS TO BE TOROIDAL RHO
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 200)
    zInterp = np.linspace(np.min(z), np.max(z), 200)
    psirzNormInterp = interp2d(r,z, psirzNorm, kind = 'cubic')(rInterp, zInterp)
    
    rhosToPlot = [.2,.4,.6,.8,1]#np.arange(.1,1.1,.1)

    ax.contour(rInterp, zInterp, psirzNormInterp, np.square(rhosToPlot), colors= 'k')

def main():
    maxDelPwrPlot = .9 #what portion of ray power must have been damped before we stop plotting that ray

    #in row major order
    fig, axs = plt.subplots(3,3, figsize = (13,7))

    stem = '/home/grantr/scratch/genray_batch/WEST_shots/WEST_56898.6000/WEST_56898.6000'
    NPara_forwards = [-1.8,-2,-2.2]
    NPara_revs = [-4,-5,-6]

    whatCode = 'both'
    totalNum = -1
    print(f'starting making scans')
    for i in range(len(NPara_forwards)):
        for j in range(len(NPara_revs)):
            totalNum +=1
            targetDir = f'{stem}_n{np.abs(NPara_forwards[i])}Npara_n{np.abs(NPara_revs[j])}Npara_thgrill140_Tscale3'

            plotRaysOfDirectory(axs.flat[totalNum], targetDir, maxDelPwrPlot)

    plt.show()

main()
