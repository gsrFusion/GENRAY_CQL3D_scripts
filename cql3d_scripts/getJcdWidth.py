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
gfileDict = getGfileDict.getGfileDict()

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

import getInputFileDictionary
cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d')

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

import helperFunctions as helper
from scipy.interpolate import interp1d

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
peakCurr = np.max(curr)
rya = np.ma.getdata(cql_nc.variables["rya"][:])
R_lfs = helper.convertRhopolToRmidplane(rya, targetDir, side = 'LFS')
print(R_lfs)
R_lfs_interp = np.linspace(R_lfs[0], R_lfs[-1],500)
curr_interp = interp1d(R_lfs, curr)(R_lfs_interp)

mask = np.where(curr_interp >= peakCurr/2)
mainPeakPoints = R_lfs_interp[mask]
print(f'mainPeakPoints: {mainPeakPoints}')
depWidth = np.abs(mainPeakPoints[0] - mainPeakPoints[-1])
print(f'fwhm = {depWidth}')

avgJcd = np.sum(curr_interp[mask])/(len(curr_interp[mask]))
print(f'avg value: {avgJcd}')

fig,ax = plt.subplots()
ax.plot(rya, curr, lw = 2, label = r'$J_{LH}$')
#ax.axvline(R_lfs_interp[mask][0], lw = 2, color = 'k', linestyle = 'dashed')
#ax.axvline(R_lfs_interp[mask][-1], lw = 2, color = 'k', linestyle = 'dashed')
ax.set_xlabel(r'$R_{LFS}$')
ax.set_ylabel(r'current density J (MA/$m^2$)')
fig.tight_layout()
plt.show()
