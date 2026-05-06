###
# Plot the power deposition profiles predicted by GENRAY for several simulations
###


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.collections import LineCollection

import matplotlib
import os, sys
from scipy.signal import find_peaks
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getGfileDict
import helperFunctions as helper
import getInputFileDictionary
import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
machine = getTargetInfo.getMachine()
#netCDF4.Dataset(f'{targetDir}/genray.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 14)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 12)

machine = 'DIIID'

if machine == 'DIIID':
    #"""
    shotNum = '203619.04130'

    case = '203619 4130 prmt6 LH scan'
    if case == '203619 4130 prmt6 LH scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            f'{stem2}_1e-2prmt6LH',
            f'{stem2}_2.5e-3prmt6LH',
            f'{stem1}_second',
        ]
        
        labels = [
            r'prmt6LH = 1e-2',
            r'prmt6LH = 2.5e-3',
            r'prmt6LH = 5e-3',
        ]

elif machine == 'NTPT':
    shotNum = 'DIIID.147634PT05'

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}'

    targetDirs =[
                f'{stem}/{machine}_{shotNum}_n2.8Npara_140thgrill_1MW',
                f'{stem}/{machine}_{shotNum}_n2.8Npara_140thgrill_1MW_highResGenTest',
                f'{stem}/{machine}_{shotNum}_n2.8Npara_140thgrill_1MW_highResGenTest_101NR',
                f'{stem}/{machine}_{shotNum}_n2.8Npara_140thgrill_1MW_highResGenTest_51NR',
                ]
    
    labels = ['nrelt = 6000,prmt6 = 0.005,NR=201', 
              'nrelt = 30000,prmt6 = 0.001,NR=201',
              'nrelt = 30000,prmt6 = 0.001,NR=101',
              'nrelt = 30000,prmt6 = 0.001,NR=51',]

#labels = ['Zeff = 1', 'Zeff = 1.25', 'Zeff = 1.3']
fig,ax = plt.subplots(figsize=(10,5))

#ax.set_title('prmt4 = 0.001, rksteps = 35000')

#adds the ray traces to ax
def plotPowerDep(targetDir, label = ''):
    print(targetDir)
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')
    powden_e = genray_nc.variables["powden_e"][:]
    rho_bin_center = genray_nc.variables["rho_bin_center"][:]
    power_total = genray_nc.variables["power_total"][:]
    print(f'absorbed {power_total*1e-7} W')
    ax.plot(rho_bin_center, powden_e, lw = 3, label = label)


for i, targetDir in enumerate(targetDirs):
    plotPowerDep(targetDir,label = labels[i])

ax.set_ylabel("power (erg/s cm^3)")
ax.set_xlabel("rho_p")
ax.legend()
fig.tight_layout()
plt.show()


    

