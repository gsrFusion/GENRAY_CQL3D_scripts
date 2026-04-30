import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.collections import LineCollection
import matplotlib
from scipy.interpolate import interp2d
import matplotlib
from matplotlib.path import Path

import os, sys
from scipy.signal import find_peaks
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import helperFunctions as helper
import getGfileDict
import getInputFileDictionary
import netCDF4
#"""
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)
#"""


def plotRays(ax_pol, targetDir, maxDelPwrPlot):
    genray_in = getInputFileDictionary.getInputFileDictionary('genray', targetDir = targetDir)
    gfileDict = getGfileDict.getGfileDict(targetDir = targetDir)
    cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    toroidalAngle = cqlrf_nc.variables["wphi"][:] 
    wr  = cqlrf_nc.variables["wr"][:] #major radius of the ray at each point along the trace
    wz  = cqlrf_nc.variables["wz"][:] #height of the ray at each point along the trace
    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    wr *= .01; wz*=.01 #convert to m from cm

    radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))#radial variable. I think it's rho_pol in my case. Doesn't really matter in this application
    nparas = cqlrf_nc.variables['wnpar'][:]
    norm = plt.Normalize(0, 1)

    
    #plot the ray using a LineCollection which allows the colormap to be applied to each ray
    for ray in range(len(wr)):

        #if not(genray_in['grill'][f'anmax({3})'] >= nparas[ray][0] >= genray_in['grill'][f'anmin({3})']):
        #    continue

        bounceIndex = helper.findBounceIndex(radialVariable[ray],bounceToFind = 1)
        mostPowerDep = helper.findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]/delpwr[ray][0]) #find the index of the last ray point we want to plot
        endingIndex = mostPowerDep#min(bounceIndex,mostPowerDep)
        #if rays are in the forwrad lobe, use them to calcualte first and second pass absorption
        points = np.array([wr[ray][:endingIndex], wz[ray][:endingIndex]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap='turbo')
        # Set the values used for colormapping
        lc.set_array(delpwr[ray][:endingIndex]/delpwr[ray][0])
        lc.set_linewidth(2)
        ax_pol.add_collection(lc)

    

    #"""
    avgSPA = helper.getSPA(targetDir, lobes = [1])
    avgSPA_allLobes = helper.getSPA(targetDir, lobes = [1,2,3,4])
    print(f'avgSPA_allLobes: {avgSPA_allLobes}')
    avg2PA = helper.getNPA(cqlrf_nc, genray_in, 2,lobes = [1])
    print(f'average SPA of the forward lobe: {avgSPA}, avg2PA: {avg2PA}')

    N_para_launch = (genray_in['grill']['anmax(1)'] + genray_in['grill']['anmin(1)'])/2

    cmap = matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(0,1),
            cmap = plt.get_cmap('turbo'))
    cmap.set_array([])
    cticks = np.linspace(0,1,5)

    #ax.set_title(f"Plotting Rays until {(maxDelPwrPlot) * 100} %\n ray power deposition")
    ax_pol.set_aspect('equal')

    lim = np.array([xlim, ylim]).T
    limiterPath  = Path(lim)

    helper.drawFluxSurfaces(ax_pol, gfileDict = gfileDict, rhosToPlot = [.2,.4,.6,.8], colors = 'k', limPath = limiterPath)
    helper.drawFluxSurfaces(ax_pol, gfileDict = gfileDict, rhosToPlot = [1], colors = 'k', limPath = limiterPath)
    ax_pol.plot(xlim, ylim, color = 'grey', lw = 3)

    ax_pol.set_ylabel("Z (m)", labelpad = -10)
    ax_pol.set_xlabel("R (m)")

    return avgSPA

def main():
    maxDelPwrPlot = .9 #what portion of ray power must have been damped before we stop plotting that ray

    #in row major order
    fig, axs = plt.subplots(2,6,figsize = (14,10))

    machine = 'NTPT'
    power = 10#MW
    NPara_target = -2.3
    thgrill = 120

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/'

    prefix = 'n'
    if NPara_target > 0:
        prefix = 'p'

    shotNums = ['MANTA']*11
    shotTimes = ['.NT05','.NT04','.NT03','.NT02','.NT01','.NT00','.PT01','.PT02','.PT03', '.PT04','.PT05']

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/'

    fig,ax = plt.subplots()

    for i in range(len(shotNums)):
        shotNum = shotNums[i]
        shotTime = shotTimes[i]

        targetDir = f'{stem}{machine}_{shotNum}{shotTime}/{machine}_{shotNum}{shotTime}_{prefix}{np.abs(NPara_target)}Npara_{thgrill}thgrill_{power}MW'

        SPA = plotRays(axs.flat[i], targetDir, maxDelPwrPlot)
        axs.flat[i].set_title(f'{shotTimes[i][1:]}')

    cmap = matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(0,1),
            cmap = plt.get_cmap('turbo'))
    cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    fig.colorbar(cmap, cax=cbar_ax, shrink=.8).set_label(r"Fractional power remaining in ray")
    plt.show()

main()
