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
from cql3d_scripts import helperFunctions as helper

def plotDepOfScan():
    stem = '/home/grantr/scratch/genray_batch/WEST_shots/WEST_56898.6000/Npara_thgrill_scan/WEST_56898.6000'
    NPara_forwards = [-1.5,-1.6,-1.7,-1.8,-1.8,-1.9,-2,-2.1,-2.2,-2.3,-2.4]
    thgrills = [220]

    fig, ax = plt.subplots()
    depLoc = np.zeros(len(NPara_forwards))
    spas = np.zeros(len(NPara_forwards))
    for i in range(len(NPara_forwards)):
        Tscale = 3
        thgrill = 220
        targetDir = f'{stem}_n{np.abs(NPara_forwards[i])}Npara_thgrill{thgrill}_Tscale{Tscale}'
        print(f'targetDir: {targetDir}')
        cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
        genray_in = getInputFileDictionary.getInputFileDictionary('genray',pathprefix=f'{parentdir}/', targetDir = targetDir)
        spas[i] = helper.getSPA_forward(cqlrf_nc, genray_in)
        depLoc[i] = helper.getPeakFirstPassDepostion(cqlrf_nc, genray_in)

    ax.plot(NPara_forwards, depLoc, lw = 2)
    ax2=ax.twinx()
    ax2.plot(NPara_forwards, spas, lw = 2, linestyle = 'dashed', color = 'k')
    ax2.set_ylabel('SPA')
    ax.set_title(f'thgrill = {thgrill}')
    ax.set_ylabel('location of most absorption in first pass')
    plt.show()





plotDepOfScan()