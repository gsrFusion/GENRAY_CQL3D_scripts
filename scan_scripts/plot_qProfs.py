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

import getGfileDict
import helperFunctions as helper

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

def plotSPAMatrix():
    
    shot = '195555'
    times = ['01200','02000','03000','04000']

    
    timeLabels = [int(time) for time in times]

    stem = f'/home/grantr/scratch/genray_batch/DIIID_shots/DIIID_{shot}'

    fig, ax = plt.subplots(figsize = (8*.8,6*.8))
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    for j in range(len(times)):
        targetDir = f'{stem}.{times[j]}/DIIID_{shot}.{times[j]}_p{2.7}Npara_300kW'
        print(f'targetDir: {targetDir}')
        
        cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
        rya = cql_nc.variables["rya"][:]
        q_prof = cql_nc.variables["qsafety"][:]

        ax.plot(rya, np.abs(q_prof), color = colors[j], lw = 2, label = timeLabels[j])#plot LCFS

    ax.legend(loc = 'upper left')
    
    ax.set_ylabel("safety factor")
    ax.set_xlabel(r"$\rho_{pol}$")

    ax.set_title(rf'Shot {shot}')
    fig.tight_layout()

    plt.show()

plotSPAMatrix()