"""
Plots two CQL3D ray trajectories side by side that share a colorbar
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
    cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    wr  = cqlrf_nc.variables["wr"][:] #major radius of the ray at each point along the trace
    wz  = cqlrf_nc.variables["wz"][:] #height of the ray at each point along the trace
    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    wr *= .01; wz*=.01 #convert to m from cm

    radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))#radial variable. I think it's rho_pol in my case. Doesn't really matter in this application
    norm = plt.Normalize(0, 1)
    
    #plot the ray using a LineCollection which allows the colormap to be applied to each ray
    for ray in range(len(wr)):
        if ray % 10 > 0:
            continue

        bounceIndex = helper.findBounceIndex(radialVariable[ray],bounceToFind = 2)
        mostPowerDep = helper.findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]/delpwr[ray][0]) #find the index of the last ray point we want to plot
        endingIndex = min(bounceIndex,mostPowerDep)
        #if rays are in the forwrad lobe, use them to calcualte first and second pass absorption
        points = np.array([wr[ray][:endingIndex], wz[ray][:endingIndex]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap='turbo')
        # Set the values used for colormapping
        lc.set_array(delpwr[ray][:endingIndex]/delpwr[ray][0])
        lc.set_linewidth(3)
        ax_pol.add_collection(lc)
    #"""
    avgSPA = helper.getSPA(targetDir)
    print(f'avgSPA_allLobes: {avgSPA}')

    #N_para_launch = (genray_in['grill']['anmax(1)'] + genray_in['grill']['anmin(1)'])/2

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
    #ax_pol.plot(xlim, ylim, color = 'grey', lw = 3)

    ax_pol.set_ylabel("Z (m)", labelpad = -10)
    ax_pol.set_xlabel("R (m)")

    d = 0.05
    #ax_pol.set_yticks([-.75,-.5,-.25,0,.25,.5,.75])
    ax_pol.set_ylim([np.min(ylim)-d, np.max(ylim)+d])
    ax_pol.set_xlim([np.min(xlim)-d, np.max(xlim)+d])

    #"""
    d = 0.2
    #ax_pol.set_yticks([-.75,-.5,-.25,0,.25,.5,.75])
    ax_pol.set_ylim([np.min(zbbbs)-d, np.max(zbbbs)+d])
    ax_pol.set_xlim([np.min(rbbbs)-d, np.max(rbbbs)+d])
    #"""

    return avgSPA

def main():
    maxDelPwrPlot = .9 #what portion of ray power must have been damped before we stop plotting that ray
    machine = 'NTPT'

    if machine == 'DIIID':
        time = '04130'
        shot = '203619'
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}.{time}/Npara_height_scan/{machine}_{shot}.{time}'

        targetDirs =[
            f'{stem}_n2.9Npara_0.0grillHeight_1MW',
            f'{stem}_n2.9Npara_0.0grillHeight_1MW_LFS',
            ]
        
    elif machine == 'NTPT':    
        fakeDevice = 'ARC'

        if fakeDevice == 'DIIID':
            fakeShot = '147634'
            power = 1

            npara = 2.8
            grillHeight = -.25
            PT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.{fakeShot}PT/{machine}_{fakeDevice}.{fakeShot}PT_n{np.abs(npara)}Npara_{grillHeight}grillHeight_{power}MW'
            NT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.{fakeShot}NT/{machine}_{fakeDevice}.{fakeShot}NT_n{np.abs(npara)}Npara_{grillHeight}grillHeight_{power}MW'

            #"""
            fakeShot = '193765'
            power = 1

            npara = 2.5
            grillHeight = .4
            PT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.{fakeShot}PT/{machine}_{fakeDevice}.{fakeShot}PT_p{np.abs(npara)}Npara_{grillHeight}grillHeight_{power}MW'
            NT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.{fakeShot}NT/{machine}_{fakeDevice}.{fakeShot}NT_p{np.abs(npara)}Npara_{grillHeight}grillHeight_{power}MW'
            #"""

        elif fakeDevice == 'ARC':
            fakeShot = 'V3A'
            power = 10

            npara = -1.5
            grillHeight = 1.0
            PT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.{fakeShot}PT/{machine}_{fakeDevice}.{fakeShot}PT_n{np.abs(npara)}Npara_{grillHeight}grillHeight_{power}MW'
            NT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.{fakeShot}NT/{machine}_{fakeDevice}.{fakeShot}NT_n{np.abs(npara)}Npara_{grillHeight}grillHeight_{power}MW'

        #fig.suptitle(r'$N_{||}$ = ' + f'{npara}, ' + r'$Z_{launcher}$ = ' + f'{grillHeight} m')
        #fig.subplots_adjust(right=0.85, wspace=0.3)

        # Add a manually placed colorbar
        
        targetDirs = [PT_target, NT_target]

    if machine == 'NTPT':
        if fakeShot == '193765':
            fig, axs = plt.subplots(1,2,figsize = (9,5))    
        if fakeShot == '147634':
            fig, axs = plt.subplots(1,2,figsize = (7.5,5))
        if fakeShot == 'V3A':
            fig, axs = plt.subplots(1,2,figsize = (7.5,5))
    else:
        fig, axs = plt.subplots(1,2,figsize = (7.5,5))

    for i, targetDir in enumerate(targetDirs):
        SPA = plotRays(axs.flat[i], targetDir, maxDelPwrPlot)

    cmap = matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(0,1),
            cmap = plt.get_cmap('turbo'))
    
    fig.colorbar(cmap, ax=axs[1]).set_label(r"Fractional power remaining in ray")
    fig.tight_layout()

    #plt.savefig('toka_ARC_HFS_LFS.jpeg',dpi=300)
    plt.show()

main()
