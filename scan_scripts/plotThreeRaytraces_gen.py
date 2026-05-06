"""
Plots three GENRAY raytraces side by side
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib
from matplotlib.path import Path
import netCDF4
import os, sys

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import helperFunctions as helper
import getGfileDict

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

def plotRays(ax_pol, targetDir, maxDelPwrPlot):
    gfileDict = getGfileDict.getGfileDict(targetDir = targetDir)
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    wr  = genray_nc.variables["wr"][:] #major radius of the ray at each point along the trace
    wz  = genray_nc.variables["wz"][:] #height of the ray at each point along the trace
    delpwr= genray_nc.variables["delpwr"][:] #power in the ray at each point
    wr *= .01; wz*=.01 #convert to m from cm

    radialVariable = (np.copy(genray_nc.variables["spsi"]))#radial variable. I think it's rho_pol in my case. Doesn't really matter in this application
    norm = plt.Normalize(0, 1)
    
    #plot the ray using a LineCollection which allows the colormap to be applied to each ray
    for ray in range(len(wr)):
        if ray is not 10:
            continue
        #if not(genray_in['grill'][f'anmax({1})'] >= nparas[ray][0] >= genray_in['grill'][f'anmin({1})']):
        #    continue

        bounceIndex = helper.findBounceIndex(radialVariable[ray],bounceToFind = 3)
        mostPowerDep = helper.findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]/delpwr[ray][0]) #find the index of the last ray point we want to plot
        endingIndex = min(bounceIndex,mostPowerDep)
        #if rays are in the forwrad lobe, use them to calcualte first and second pass absorption
        points = np.array([wr[ray][:endingIndex], wz[ray][:endingIndex]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap='turbo')
        # Set the values used for colormapping
        lc.set_array(delpwr[ray][:endingIndex]/delpwr[ray][0])
        lc.set_linewidth(2)
        ax_pol.add_collection(lc)

    avgSPA = helper.getSPA(targetDir)
    print(f'avgSPA: {avgSPA}')

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
    ax_pol.set_ylim(min(ylim)-.05, max(ylim)+.05)
    ax_pol.set_xlim(min(xlim)-.05, max(xlim)+.05)

    return avgSPA

def main():
    maxDelPwrPlot = .9 #what portion of ray power must have been damped before we stop plotting that ray

    #in row major order
    #fig = plt.figure(figsize = (10,5.25))
    fig, axs = plt.subplots(1, 3, figsize = (10,5.25))


    time = '.04135'
    shot = '203619'
    machine = 'DIIID'

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/gridTests/DIIID_{shot}{time}_expSpectrum'
    print(f'starting making scans')
    
    targetDirs = [
        f'{stem}_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.001prmt4',
        f'{stem}_1Zeff_cqlHighRes_30000nrelt_0.0025prmt6_0.001prmt4',
        f'{stem}_1Zeff_cqlHighRes_30000nrelt_0.001prmt6_0.001prmt4',
    ]

    labels = [
        r'prmt6=0.005',
        r'prmt6=0.0025',
        r'prmt6=0.001',
    ]
    shotNum = f'{shot}{time}'
    targetDirs =[
        f'{stem}_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.001prmt4',
        f'{stem}_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.0005prmt4',
        f'{stem}_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.0001prmt4',
        
        ]
    
    labels = [
              'prmt4 = 0.001', 
              'prmt4 = 0.0005', 
              'prmt4 = 0.0001', 
              ]


    for i, targetDir in enumerate(targetDirs):
        SPA = plotRays(axs[i], targetDir, maxDelPwrPlot)
        axs[i].set_title(labels[i])
        #axs.flat[i].set_title(f'SPA = {SPA:.2f}\n N|| = {NPara_for}')
    #fig.suptitle(rf'$N_{{||,LCFS}} = {NPara_for}$')

    fig.tight_layout(rect=[0, 0, 0.9, 1])

    cmap = matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(0,1),
            cmap = plt.get_cmap('turbo'))
    bbox = axs[2].get_position()

    # Manually add colorbar axis: narrow and same height as ax3's inner box
    # Place it slightly to the right of ax3
    cbar_width = 0.02
    pad = 0.01
    cax = fig.add_axes([
        bbox.x1 + pad,       # left
        bbox.y0,             # bottom
        cbar_width,          # width
        bbox.height          # height
    ])
    fig.colorbar(cmap, cax=cax).set_label(r"Fractional power remaining in ray")
    
    plt.show()

main()
