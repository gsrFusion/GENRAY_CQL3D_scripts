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

def plotSPAMatrix():
    
    shot = '195557'
    #times = ['01400', '02500', '04500', '05750']
    #times = ['01500', '03000', '04000', '05000']
    #times = ['01200', '01750', '02500', '03500' , '04500']
    #times = ['01500', '02500', '03500', '04500']
    #times = ['01300', '02500', '03650', '04900']
    #times = ['02000','03000','04500']
    #times = ['01500','03000','04500']
    #times = ['01500','02750','03300']
    times = ['01200','02000','03000','04500']
    timeLabels = [int(time) for time in times]

    stem = f'/home/grantr/scratch/genray_batch/DIIID_shots/DIIID_{shot}'

    fig, ax = plt.subplots(figsize = (8*.8,6*.8))

    for j in range(len(times)):
        targetDir = f'{stem}.{times[j]}/DIIID_{shot}.{times[j]}_p{2.7}Npara_300kW'
        print(f'targetDir: {targetDir}')

        inputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d', targetDir=targetDir)

        enescal = 1

        try:
            enescal = inputFileDict['setup']['enescal']
        except:
            pass
        n_e = inputFileDict['setup']['enein(1,1)']*1e6*enescal
        rhos = inputFileDict['setup']['ryain']

        ax.plot(rhos, np.array(n_e/1e19), label = f't = {timeLabels[j]}', linewidth = 3)

    ax.legend(loc = 'best')

    ax.set_xlabel(r'$\rho_{pol}$', fontsize = 20)
    ax.set_xlim([0,1])
    ax.set_ylim([0,10])
    ax.set_ylabel(r'n$_e$ (10$^{19}$/m$^{-3}$)')
    ax.grid()

    ax.set_title(rf'CQL3D {shot} n$_{{e}}$ profiles'+f'\n enescal = {enescal}')
    fig.tight_layout()
    plt.show()

plotSPAMatrix()