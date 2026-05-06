###
# Plots delpwr for a particular ray across several simulations
# Liekly only useful if the exact same spectrum is launched
###


import numpy as np
import matplotlib.pyplot as plt
import os, sys
import netCDF4

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import helperFunctions as helper

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

    ray = 10

    stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_203619.04130/DIIID_203619.04130_expSpectrum_2Zeff/DIIID_203619.04130_expSpectrum_2Zeff_'
    stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_203912.02700/DIIID_203912.02700_expSpectrum_'
    targetDirs = [f'{stem1}first',
                  f'{stem2}second',]


    labels = [
              '203619',
              '203912' 
              ]
    
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


    

