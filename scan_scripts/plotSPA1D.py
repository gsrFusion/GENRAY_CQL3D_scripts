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

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)


def plotSPA1D():
    machine = 'DIIID'
    shot = '195555'
    tribots = ['', 'Tribot0.3', 'Tribot0.4', 'Tribot0.5', 'Tribot0', 'RevB']
    trilabels = [0.2, 0.3, 0.4, 0.5, 0, '0.2, rev B']

    stem = f'/home/grantr/scratch/genray_batch/DIIID_shots/DIIID_{shot}'

    NPara_fors = np.array([2.3, 2.5, 2.7, 2.9, 3.1])

    fig, ax = plt.subplots(figsize = (8*.8,6*.8))
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    #for j in range(len(times)):
    for j in range(len(tribots)):
        SPA_scan = np.zeros(len(NPara_fors))
        for i in range(len(NPara_fors)):
            prefix = 'n'
            if NPara_fors[i] > 0:
                prefix = 'p'
            if tribots[j] == 'RevB':
                prefix = 'n'
            targetDir = f'{stem}{tribots[j]}.03000/DIIID_{shot}{tribots[j]}.03000_{prefix}{NPara_fors[i]}Npara_300kW'

            cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
            genray_in = getInputFileDictionary.getInputFileDictionary('genray', targetDir = targetDir)
            SPA_scan[i] = helper.getSPA(cqlrf_nc, genray_in, lobes = [1])

        ax.scatter(NPara_fors, SPA_scan, marker ='s', s = 40, label = r'$\delta = $' + f'{trilabels[j]}', color = colors[j], zorder= 10)

    ax.grid()
    
    ax.set_xlabel(r'sign(-B)*N$_{||,for}$')
    ax.set_ylabel(r'SPA')
    ax.set_xticks(NPara_fors)
    ax.set_ylim([0,1.05])
    ax.legend(loc = 'best', ncol = 2)
    ax.set_title(rf'Launching 300 kW')
    fig.tight_layout()

    plt.show()

plotSPA1D()