###
# Plots the electron and ion densities and temperatures according to the cql3d input file
###

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__);dname = os.path.dirname(abspath);os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__));parentdir = os.path.dirname(currentdir);sys.path.append(parentdir)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
import helperFunctions as helper
import getInputFileDictionary
cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d')

import numpy as np
import matplotlib.pyplot as plt
import netCDF4

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
rya = np.ma.getdata(cql_nc.variables["rya"][:])
darea = np.ma.getdata(cql_nc.variables["darea"][:])

_, T_e = helper.getCQLTe(rho_pol = rya)
_, n_e = helper.getCQLne(rho_pol = rya)

dischargeNumber = os.getcwd().split('/')[-1]
if len(dischargeNumber) != 6:
    dischargeNumber = os.getcwd().split('/')[-2]


#####Setup and do plotting#####
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)
plt.rc('legend',fontsize=16)

fig, ax = plt.subplots(figsize = (8*.8,6*.8))

ax2 = ax.twinx()

Tcolor = 'tab:orange'
Ncolor = 'tab:blue'

ax.plot(rya, T_e/n_e, label = r'$T_e/n_e$', linewidth = 3,color = 'tab:red')
ax2.plot(rya, T_e*darea/n_e, label = r'$T_e \Delta_{area}/n_e$', linewidth = 3,color = 'tab:blue')

ax.yaxis.label.set_color('tab:red')
ax2.yaxis.label.set_color('tab:blue')


ax.set_ylabel(r'$T_e/n_e$')
ax2.set_ylabel(r'$T_e \Delta_{area}/n_e$')

ax.set_xlabel(r'$\rho_{pol}$', fontsize = 20)
ax.set_xlim([0,1])
ax.grid()

fig.suptitle(f'Shot {shotNum}')
fig.tight_layout(rect=[0, 0.0, 1, 1.05])
plt.show()
