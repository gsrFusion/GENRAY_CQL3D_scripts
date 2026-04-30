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

    ray = 10

    stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}'
    
    targetDirs =[
                f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.001prmt4',
                f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.0025prmt6_0.001prmt4',
                f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.001prmt6_0.001prmt4',
                f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.0008prmt6_0.001prmt4',
                ]
    
    labels = [
              'prmt6 = 0.005, prmt4 = 0.001', 
              'prmt6 = 0.0025, prmt4 = 0.001', 
              'prmt6 = 0.001, prmt4 = 0.001', 
              'prmt6 = 0.0008, prmt4 = 0.001', 
              ]
    

    targetDirs =[
            f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.001prmt4',
            f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.0005prmt4',
            f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.0001prmt4',
            f'{stem}/gridTests/DIIID_{shotNum}_expSpectrum_1Zeff_cqlHighRes_30000nrelt_0.005prmt6_0.00005prmt4',
            
            ]
    
    labels = [
              'prmt6 = 0.005, prmt4 = 0.001', 
              'prmt6 = 0.005, prmt4 = 0.0005', 
              'prmt6 = 0.005, prmt4 = 0.0001', 
              'prmt6 = 0.005, prmt4 = 0.00005', 
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
def plotPowerDep(targetDir, ray, label = ''):
    print(targetDir)
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')
    delpwr= genray_nc.variables["delpwr"][:] #power in the ray at each point
    spsi = (np.copy(genray_nc.variables["spsi"]))
    ws = genray_nc.variables["ws"][:] # poloidal length along ray

    maxPowerToPlot = 0.95
    mostPowerDep = helper.findNearestIndex(1 - maxPowerToPlot, delpwr[ray]/delpwr[ray][0])

    ax.plot(ws[ray][:mostPowerDep], delpwr[ray][:mostPowerDep], lw = 3, label = label)


for i, targetDir in enumerate(targetDirs):
    plotPowerDep(targetDir, ray, label = labels[i])

ax.set_ylabel("delpwr")
ax.set_xlabel("poloidal distance along ray")
ax.legend()
fig.tight_layout()
plt.show()


    

