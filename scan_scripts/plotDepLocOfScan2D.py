import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os, sys
import matplotlib
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getInputFileDictionary
import helperFunctions as helper

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.rc('axes', labelsize = 18)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 16)


import numpy as np
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import netCDF4
print(f'past imports')


machine = 'NTPT'

if machine == 'DIIID':
    NPara_fors = np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
    thgrills = np.array([120,140,160,180,200,220,240])
    grillHeights = np.round(np.linspace(-.5,.5,13),3)
    #NPara_fors = np.array([2.3,2.5,2.7,2.9,3.1])#-1*np.array([1.5,1.7,1.9,2.1,2.3,2.5])
    #thgrills = np.array([240,220,200,180,160,140,120])
    time = '.05500'
    shot = '193765'
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/Npara_height_scan/{machine}_{shot}{time}'

if machine == 'NTPT':
    time = '.V3APT'
    shot = 'ARC'
    if '193765' in time:
        grillHeights = np.round(np.linspace(-.5,.5,11),3)
        NPara_targets = np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
        power = 1

    elif '147634' in time:
        grillHeights = np.round(np.linspace(-.75,.75,13),3)
        NPara_targets = -1*np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
        power = 1

    elif 'V3A' in time:
        grillHeights = np.round(np.linspace(-1.75,1.75,15),3)
        NPara_targets = -1*np.array([1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2])

        power = 10


    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/LFSVersion/{machine}_{shot}{time}'

def plotDepMatrix():
    depMatrix = np.zeros((len(NPara_targets), len(grillHeights)))
    counter = 0
    for i in range(len(NPara_targets)):
        NPara_for = NPara_targets[i]
        prefix = 'n'
        if np.sign(NPara_for) > 0:
            prefix = 'p'
        for j in range(len(grillHeights)):
            grillHeight = grillHeights[j]
            targetDir = f'{stem}_{prefix}{np.abs(NPara_for)}Npara_{grillHeights[j]}grillHeight_{power}MW'

            try:
                print(f'targetDir: {targetDir}')
                SPA = helper.getSPA(targetDir)[0]
                print(f'SPA = {SPA}')
                print(f'{SPA >= 0.9}')

                if SPA >= 0.9:
                    depMatrix[i,j],_ = helper.getAvgCurrentLocAndTotal(targetDir)
                    counter += 1
                else:
                    depMatrix[i,j] = np.nan
            except Exception as e:
                import traceback
                print(traceback.format_exc())
                depMatrix[i,j] = np.nan

    print(counter)
    print(f'min, max: {np.nanmin(depMatrix), np.nanmax(depMatrix)}')
    fig,ax = plt.subplots(figsize=(7.5,5.5))
    p = ax.pcolormesh(grillHeights, NPara_targets, depMatrix,shading = 'nearest', 
                       cmap='inferno_r', vmin=0.5, vmax = 1)
    
    print(f'min loc at: {np.nanmin(depMatrix)}')

    ax.set_ylabel(r'N$_{||,LCFS}$')
    ax.set_xlabel(r'$Z_{launcher}$ (m)')
    ax.set_yticks(NPara_targets)
    ax.set_xticks(grillHeights[::2])
    ax.tick_params(axis='x', rotation=30)
    ax.yaxis.get_label().set_fontsize(16)

    triString = ''
    if 'PT' in time:
        triString = 'Positive'
    else:
        triString = 'Negative'
    ax.set_title(f'{triString} triangularity {time[1:-2]}-like LFS Launch')
    #ax.set_title(r'Positive triangularity V3A-like ($\delta = 0.5$)')

    cbar = fig.colorbar(p, ax = ax, shrink = .9, pad = .01)
    #cbar.set_label(r'$\int \rho_{pol} \cdot J_{LH}\, dA \left/ \int J_{LH}\, dA \right.$')
    cbar.set_label(r'Current centroid ($\rho_{pol}$)')
    

    fig.tight_layout()
    if 'PT' in time:
        plt.savefig('toka_V3A_PT_LFS_depLoc.jpeg',dpi=300)
        pass
    if 'NT' in time:
        #plt.savefig('toka_V3A_NT_depLoc.jpeg',dpi=300)
        pass

    plt.show()

plotDepMatrix()