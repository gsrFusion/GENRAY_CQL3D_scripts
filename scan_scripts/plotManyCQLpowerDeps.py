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

    shotNum = '203619.04135'

    
    scanType = 'timeslices'

    if scanType == 'timeslices':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots'
        targetDirs = [
            f'{stem}/{machine}_203619.04120/{machine}_203619.04120_expSpectrum_2Zeff_second',
            f'{stem}/{machine}_203619.04135/{machine}_203619.04135_expSpectrum_2Zeff',
            f'{stem}/{machine}_203619.04160/{machine}_203619.04160_expSpectrum_2Zeff'
        ]
        labels = [
            '4120',
            '4135',
            '4160'
        ]

    if scanType == 'nnkpar':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_1Zeff_1Lobe_20000nrelt_0.005prmt6_0.001prmt4'
        targetDirs = [
            f'{stem}_4nthin_15nnkpar',
            f'{stem}_8nthin_15nnkpar',
            f'{stem}_4nthin_20nnkpar',
            f'{stem}_4nthin_30nnkpar',
            f'{stem}_4nthin_35nnkpar',
            f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_1.25Zeff_1Lobe_10000nrelt_0.005prmt6_0.001prmt4_4nthin_35nnkpar',
        ]

        labels = [
            'nthin = 4, nnkpar = 15',
            'nthin = 8, nnkpar = 15',
            'nthin = 4, nnkpar = 20',
            'nthin = 4, nnkpar = 30',
            'nthin = 4, nnkpar = 35',
            'nthin = 4, nnkpar = 35, Zeff = 1.25',
        ]

    if scanType == 'expSpectrum nnkpar=30':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}'
        targetDirs = [
            f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_2e-6prmt4_4nthin_30nnkpar_1e-4prmt4ECE_correctWall_newYuri'
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.004prmt6_1e-8prmt4_4nthin_30nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_35nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.0045prmt6_1e-8prmt4_4nthin_35nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_37nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-7prmt4_4nthin_39nnkpar_noDebug',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_39nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_43nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_45nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_50nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_55nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_1e-8prmt4_4nthin_60nnkpar',
            #[f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_2e-5prmt4_4nthin_39nnkpar'],

        ]
        labels = [
            'Zeff = 1, nnkpar = 30, prmt4=2e-6',
            #'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 30, prmt4=1e-8, ',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 35, prmt4=1e-8',
            #'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 35, prmt4=1e-8,prmt6=0.0045',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 37, prmt4=1e-8',
            #'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 39, prmt4=1e-7',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 39, prmt4=1e-8',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 43, prmt4=1e-8',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 45, prmt4=1e-8',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 50, prmt4=1e-8',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 55, prmt4=1e-8',
            'expSpectrum, Zeff = 1, nthin = 4, nnkpar = 60, prmt4=1e-8',
        ]

    if scanType == 'nthin':
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/numRaysTest/{machine}_{shotNum}_1Zeff_1Lobe_20000nrelt_0.005prmt6_0.001prmt4'
        targetDirs = [
           
        ]

        labels = [

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
fig,ax = plt.subplots()

#ax.set_title('prmt4 = 0.001, rksteps = 35000')

#adds the ray traces to ax
def plotPowerDep(targetDir, label = ''):
    print(targetDir)
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    powrft = cql_nc.variables["powrft"][:][-1,:]
    rya = cql_nc.variables["rya"][:]

    ax.plot(rya, powrft, lw = 3, label = label)


for i, targetDir in enumerate(targetDirs):
    plotPowerDep(targetDir,label = labels[i])

ax.set_ylabel("power (W/cm^3)")
ax.set_xlabel("rho_p")
ax.legend()
fig.tight_layout()
plt.show()


    

