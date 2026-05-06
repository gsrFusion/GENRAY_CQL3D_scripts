"""
Plots a matrix where the x and y axes designate a simulation run, and the colorbar is the driven current
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getInputFileDictionary
import helperFunctions as helper

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.rc('axes', labelsize = 18)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 14)


import numpy as np
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import netCDF4
print(f'past imports')

machine = 'DIIID'

minSPA = 0.8

if machine == 'DIIID':
    time = '.04525'
    shot = '147634'
    NPara_fors = -1*np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])


Imatrix = np.zeros((len(NPara_fors), len(NPara_fors)))

for i in range(len(NPara_fors)):
    N1 = NPara_fors[i]
    prefix = 'n'
    if N1 > 0:
        prefix = 'p'

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/twoColorScans/{machine}_{shot}{time}'
    for j in range(i+1):
        N2 = NPara_fors[j]
        targetDir = f'{stem}_{prefix}{np.abs(N2)}Npara_{prefix}{np.abs(N1)}Npara_1MW'

        print(f'targetDir: {targetDir}')
        cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
        genray_in = getInputFileDictionary.getInputFileDictionary('genray_LH', targetDir = targetDir)

        SPA, _, _ = helper.getSPA(targetDir)
        if SPA[0] < minSPA:
            Imatrix[i,j] = np.nan
            Imatrix[j,i] = np.nan
        else:
            cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
            curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
            darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
            rya = cql_nc.variables["rya"][:]
            totalCD = np.sum(curr*darea)

            Imatrix[i,j] = totalCD
            Imatrix[j,i] = totalCD

fig,ax = plt.subplots()
cmap = 'viridis'
pcolor = ax.pcolormesh(NPara_fors, NPara_fors, Imatrix ,shading = 'nearest',cmap=cmap)

print(f'average: {np.nanmean(Imatrix)}')

#"""
cbar = fig.colorbar(pcolor, ax = ax, shrink = .9, pad = .01)
cbar.set_label(r'$I_{LH} (MA)$')
#"""
ax.set_ylabel(r'N$_{||,1}$')
ax.set_ylabel(r'N$_{||,2}$')
ax.set_yticks(NPara_fors)
ax.set_xticks(NPara_fors)

fig.tight_layout()
plt.show()
