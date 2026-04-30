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

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

def plotSPAMatrix():


    NPara_fors = [-1.5,-1.6,-1.7,-1.8,-1.9,-2,-2.1,-2.2,-2.3,-2.4,-2.5,-2.6,-2.7]
    thgrill = 220
    #NPara_for = -2.55
    #stem = f'/home/grantr/scratch/genray_batch/WEST_shots/WEST_56898.6000/n{np.abs(NPara_for)}Npara_thgrill220_revLobeNparaAndPower_scan/WEST_56898.6000_n{np.abs(NPara_for)}Npara_thgrill220_500kW'
    #NPara_revs = [-3,-3.1,-3.2,-3.3,-3.4,-3.5,-3.6,-3.7,-3.8,-3.9,-4]
    #revPowers = [100000,150000,200000,300000,500000]
    SPA_scan = np.zeros(len(NPara_fors))
    print(f'starting making scans')
    for i in range(len(NPara_fors)):
        stem = f'/home/grantr/scratch/genray_batch/WEST_shots/WEST_56898.6000/Npara_thgrill_scan/WEST_56898.6000_n{np.abs(NPara_fors[i])}Npara'

        Tscale = 3
        targetDir = f'{stem}_thgrill{thgrill}_Tscale{Tscale}'
        print(f'targetDir: {targetDir}')
        cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
        genray_in = getInputFileDictionary.getInputFileDictionary('genray',pathprefix=f'{parentdir}/', targetDir = targetDir)
        SPA_scan[i] = helper.getSPA_forward(cqlrf_nc, genray_in)

    #SPAmatrix[SPAmatrix < 0.95] = np.NAN

    fig,ax = plt.subplots()
    ax.grid()
    ax.scatter(NPara_fors, SPA_scan, marker ='s', s = 40)
    ax.set_xlabel(r'N$_{||,for}$')
    ax.set_ylabel(r'SPA')
    ax.set_xticks([-1.5,-1.7,-1.9,-2.1,-2.3,-2.5,-2.7])
    ax.set_ylim([0,1])
    ax.set_yticks(np.linspace(0,1,11))
    ax.set_axisbelow(True)
    #ax.tick_params(axis='x', labelrotation=25)
    
    #ax.set_ylabel(r'N$_{||,2}$')
    #ax.set_xlabel(r'P$_{||,2}$ (kW)')

    ax.set_title(rf'Launching 500 kW, no reverse lobe')
    fig.tight_layout()

    plt.show()

plotSPAMatrix()