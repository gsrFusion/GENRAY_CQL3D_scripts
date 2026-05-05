###
# Simple script for plotting the q profile as a function of rho_pol
###
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)

import matplotlib.pyplot as plt
import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

print(f'targetDir: {targetDir}')
cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
rya = cql_nc.variables["rya"][:]
q_prof = cql_nc.variables["qsafety"][:]

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

fig,ax = plt.subplots()
ax.plot(rya, q_prof, lw = 2)
ax.axhline(1, color = 'k')
ax.set_xlabel(r'$\rho_{pol}$')
ax.set_ylabel(f'safety factor')
ax.grid()
fig.tight_layout()
plt.show()