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


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib
import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
from omfit_classes import omfit_eqdsk

print(f'targetDir: {targetDir}')
#cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
#rya = cql_nc.variables["rya"][:]
#q_prof = cql_nc.variables["qsafety"][:]

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

eq_147634_PT = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/g203912.02700_deepPort')
print(f'{eq_147634_PT["BCENTR"], eq_147634_PT["CURRENT"], }')
eq_147634_PT.plot()
plt.show()

fig,ax = plt.subplots()
ax.plot(np.sqrt(eq_147634_PT['fluxSurfaces']['levels']), np.abs(eq_147634_PT['fluxSurfaces']['avg']['q']), lw = 3)
#ax.plot(rya, np.abs(q_prof), lw = 2)
ax.axhline(1, color = 'k')
ax.set_xlabel(r'$\rho_{pol}$')
ax.set_ylabel(f'safety factor')
ax.grid()
fig.tight_layout()
plt.show()