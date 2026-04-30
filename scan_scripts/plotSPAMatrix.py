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

plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.rc('axes', labelsize = 20)
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
    NPara_targets = -1*np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
    thgrills = np.array([120,140,160,180,200,220,240])
    grillHeights = np.round(np.linspace(-.75,.75,13),3)
    #NPara_fors = np.array([2.3,2.5,2.7,2.9,3.1])#-1*np.array([1.5,1.7,1.9,2.1,2.3,2.5])
    #thgrills = np.array([240,220,200,180,160,140,120])
    time = '.04565'
    shot = '147634'
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/Npara_height_scan/{machine}_{shot}{time}'
    power = 1

if machine == 'NTPT':
    time = '.193765NT'
    shot = 'DIIID'
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
        NPara_targets = -1*np.array([1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2])
        power = 10


    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{machine}_{shot}{time}'

def plotSPAMatrix():
    SPAmatrix = np.zeros((len(NPara_targets), len(grillHeights)))

    for i in range(len(NPara_targets)):
        NPara_for = NPara_targets[i]
        prefix = 'n'
        if np.sign(NPara_for) > 0:
            prefix = 'p'
        for j in range(len(grillHeights)):
            #thgrill = thgrills[j]
            grillHeight = grillHeights[j]
            targetDir = f'{stem}_{prefix}{np.abs(NPara_for)}Npara_{grillHeights[j]}grillHeight_{power}MW'

            try:
                print(f'targetDir: {targetDir}')
                SPA = helper.getSPA(targetDir)[0]

                SPAmatrix[i,j] = SPA
                print(f'SPA = {SPAmatrix[i,j]}')
            except:
                print(f'failed')
                SPAmatrix[i,j] = np.nan


    
    fig,ax = plt.subplots(figsize=(7.25,4.8))
    p = ax.pcolormesh(grillHeights, NPara_targets, SPAmatrix*100,shading = 'nearest', 
                       cmap='viridis', vmin=0, vmax = 100)
    
    ax.set_ylabel(r'N$_{||,LCFS}$')
    #ax.set_xlabel(r'time (ms)')
    ax.set_xlabel(r'$Z_{launcher}$ (m)')
    ax.set_yticks(NPara_targets)
    ax.set_xticks(grillHeights[::2])
    #ax.yaxis.get_label().set_fontsize(16)
    #plt.xticks(rotation=-60)

    cbar = fig.colorbar(p, ax = ax, shrink = 1, pad = .01)
    #cbar.set_label(r'SPA$_{forward}$')
    cbar.set_label(r'Single pass absorption (percent)', fontsize = 16)
    
    #ax.set_title(r'$\delta = -0.5$')
    #ax.set_title(rf'{shot}{time}, $N_{{||}}$ and thgrill scan')
    #ax.set_title(f'{machine} shot {shot} at {time[2:]} ms')
    triString = ''
    if 'PT' in time:
        triString = 'Positive'
    else:
        triString = 'Negative'
    ax.set_title(f'{triString} triangularity {time[1:-2]}-like')

    fig.tight_layout()
    plt.savefig(f'toka_{time[1:-2]}_{time[-2:]}_SPA_shaped.jpg',dpi=300)

    plt.show()

plotSPAMatrix()