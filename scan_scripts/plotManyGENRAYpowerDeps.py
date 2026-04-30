###
# Plot the ray trajectories and damping as predicted by GENRAY
# since it's a linear code, the damping is very wrong for LH, but the ray trajectory can be useful
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

    scanType = 'rhopsi'

    if scanType == 'rhopsi':
        stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_203619.04130/DIIID_203619.04130_expSpectrum_2Zeff'
        targetDirs = [
            f'{stem}/DIIID_203619.04130_expSpectrum_2Zeff_second',
            f'{stem}/scans/DIIID_203619.04130_expSpectrum_2Zeff_second_0.99rhopsi0'
        ]

        labels = [
            'rhopsi0 = 1',
            'rhopsi0 = 0.99',

        ]

    if scanType == 'testNewGENRAY':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_expSpectrum_1Zeff_10000nrelt'
        targetDirs = [
            f'{stem}_0.005prmt6_1e-7prmt4_4nthin_39nnkpar_id16newTest',
            #f'{stem}_0.005prmt6_1e-7prmt4_4nthin_39nnkpar',
        ]

        labels = [
            'old Sam genray',
            'new Sam Genray',

        ]
    if scanType == 'expSpectrum nnkpar=33 prmt4':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_expSpectrum_1Zeff_10000nrelt'
        targetDirs = [
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_33nnkpar',
            f'{stem}_0.005prmt6_0.0001prmt4_4nthin_33nnkpar',
        ]

        labels = [
            'nthin = 4, nnkpar = 33, prmt4 = 1e-3',
            'nthin = 4, nnkpar = 33, prmt4 = 1e-4',

        ]

    if scanType == 'expSpectrum nnkpar=30':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}'
        targetDirs = [
            f'{stem}_expSpectrum_1Zeff_10000nrelt_0.004prmt6_1e-8prmt4_4nthin_30nnkpar',
            f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_35nnkpar',
            f'{stem}_expSpectrum_1Zeff_10000nrelt_0.0045prmt6_1e-8prmt4_4nthin_35nnkpar',
            f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_37nnkpar',
            #[f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_2e-5prmt4_4nthin_39nnkpar'],

        ]
        labels = [
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 30, prmt4=1e-8, ',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 35, prmt4=1e-8',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 35, prmt4=1e-8,prmt6=0.0045',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 37, prmt4=1e-8',
        ]

    if scanType == 'expSpectrum nnkpar':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_expSpectrum_1Zeff_10000nrelt'
        targetDirs = [
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_15nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_20nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_25nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_30nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_33nnkpar',
        ]

        labels = [
            'nthin = 4, nnkpar = 15',
            'nthin = 4, nnkpar = 20',
            'nthin = 4, nnkpar = 25',
            'nthin = 4, nnkpar = 30',
            'nthin = 4, nnkpar = 33',

        ]

    if scanType == 'old prmt4 scan':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/gridTests/{machine}_{shotNum}_1Zeff_30000nrelt'
        targetDirs = [
            f'{stem}_0.005prmt6_0.001prmt4',
            f'{stem}_0.005prmt6_0.0001prmt4',
            f'{stem}_0.005prmt6_0.00001prmt4',
            f'{stem}_0.005prmt6_0.000001prmt4',
            f'{stem}_0.005prmt6_0.0000001prmt4',
            
        ]

        labels = [
            'prmt6 = 0.005, prmt4 = 1e-3',
            'prmt6 = 0.005, prmt4 = 1e-4',
            'prmt6 = 0.005, prmt4 = 1e-5',
            'prmt6 = 0.005, prmt4 = 1e-6',
            'prmt6 = 0.005, prmt4 = 1e-7',

        ]

    if scanType == 'old prmt6 scan':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/gridTests/{machine}_{shotNum}_1Zeff_30000nrelt'
        targetDirs = [
            f'{stem}_0.005prmt6_0.001prmt4',
            f'{stem}_0.0025prmt6_0.001prmt4',
            f'{stem}_0.001prmt6_0.001prmt4',
            f'{stem}_0.0008prmt6_0.001prmt4',
            
        ]

        labels = [
            'prmt6 = 0.005, prmt4 = 0.001',
            'prmt6 = 0.0025, prmt4 = 0.001',
            'prmt6 = 0.001, prmt4 = 0.001',
            'prmt6 = 0.0008, prmt4 = 0.001',

        ]

    if scanType == 'nnkpar':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_1Zeff_1Lobe_20000nrelt'
        targetDirs = [
            #f'{stem}_0.005prmt6_0.001prmt4_4nthin_7nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_10nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_15nnkpar',
            #f'{stem}_0.005prmt6_0.001prmt4_8nthin_15nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_20nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_30nnkpar',
            f'{stem}_0.005prmt6_0.0005prmt4_4nthin_30nnkpar',
            f'{stem}_0.005prmt6_0.001prmt4_4nthin_35nnkpar',
            f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_1Zeff_1Lobe_10000nrelt_0.005prmt6_0.001prmt4_4nthin_35nnkpar',
            
        ]

        labels = [
            #'nthin = 4, nnkpar = 7',
            'nthin = 4, nnkpar = 10, prmt4 = 0.005',
            'nthin = 4, nnkpar = 15, prmt4 = 0.005',
            #'nthin = 8, nnkpar = 15',
            'nthin = 4, nnkpar = 20, prmt4 = 0.005',
            'nthin = 4, nnkpar = 30, prmt4 = 0.005',
            'nthin = 4, nnkpar = 30, prmt4 = 0.0005',
            'nthin = 4, nnkpar = 35, prmt4 = 0.005',
            'nthin = 4, nnkpar = 35, prmt4 = 0.005, nrelt = 10000',

        ]

    if scanType == 'nthin':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_1Zeff_1Lobe_20000nrelt_0.005prmt6_0.001prmt4'
        targetDirs = [
            f'{stem}_4nthin_7nnkpar',
            f'{stem}_4nthin_10nnkpar',
            f'{stem}_8nthin_10nnkpar',
            f'{stem}_12nthin_10nnkpar',
            
        ]

        labels = [
            'nthin = 4, nnkpar = 7',
            'nthin = 4, nnkpar = 10',
            'nthin = 8, nnkpar = 10',
            'nthin = 12, nnkpar = 10',

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


    

