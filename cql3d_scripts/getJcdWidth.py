###
# Simple script to look at current profile and print out FWHM in meters
###

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os, sys
from scipy.interpolate import interp1d
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getGfileDict
import getTargetInfo
import helperFunctions as helper

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

gfileDict = getGfileDict.getGfileDict()

targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

#LH current density
curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
peakCurr = np.max(curr)
#rya are the radial points (rho_pol) where the distribution function is calculated
rya = np.ma.getdata(cql_nc.variables["rya"][:])
#convert these rho_pol points to major radius points on the LFS midplane
R_lfs = helper.convertRhopolToRmidplane(rya, targetDir, side = 'LFS')

fwhm = helper.getJcdWidth(fracOfPeak = 0.5, mainPeak = True)
print(f'fwhm: {fwhm} m')

fig,ax = plt.subplots()
ax.plot(R_lfs, curr, lw = 2, label = r'$J_{LH}$')
ax.set_xlabel(r'$R_{LFS}$')
ax.set_ylabel(r'current density J (MA/$m^2$)')
fig.tight_layout()
plt.show()

