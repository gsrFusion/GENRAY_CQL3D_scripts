import numpy as np
import matplotlib.pylab as plt
import netCDF4
import os, sys
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
shotNum = getTargetInfo.getShotNum()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

plt.rc('xtick', labelsize = 15)
plt.rc('ytick', labelsize = 15)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)
plt.rc('legend',fontsize=16)

#distribution function
f = np.ma.getdata(cql_nc.variables["f"])*1e6#convert to 1/m^3
#pitch angles mesh at which f is defined
pitchAngleMesh = np.ma.getdata(cql_nc.variables["y"][:])
rya = np.ma.getdata(cql_nc.variables["rya"][:])

fastEnergy = 100 #keV

normalizedVel = np.ma.getdata(cql_nc.variables["x"][:])
cql3dEnergies = np.ma.getdata(cql_nc.variables["enerkev"][:])

ne_fast = np.zeros(len(rya))
indices =  []
if len(indices) == 0:
    indices = np.where(cql3dEnergies < fastEnergy)[0]
minCQL3DEnergyIndex = indices[-1]
fRelevant = f[:, minCQL3DEnergyIndex:, :]

for i in range(0, len(rya)):
    #this integrates the distribution function at each radial position
    #this is effectively a spherical jacobian 
    integFOverVel = np.ma.getdata(np.trapz(fRelevant[i,:,:]*normalizedVel[minCQL3DEnergyIndex:, None]**2, 
        normalizedVel[minCQL3DEnergyIndex:, None], axis = 0))
    ne_fast[i] = 2*np.pi*np.trapz(integFOverVel*np.sin(pitchAngleMesh[i]), pitchAngleMesh[i], axis = 0)

rho_pol, n_e = helper.getCQLne()

fig,ax = plt.subplots()
ax.plot(rya, ne_fast, label = r'n$_e$(E > ' + f'{fastEnergy} keV)', lw = 2.5)
ax.plot(rho_pol, n_e, label = r'n$_e$(E > 0 keV)', lw = 2.5)
ax.set_yscale('log')
ax.set_ylim([1e13,1e20])
ax.set_ylabel('electon density (1/m^3)')
ax.set_xlabel(r'$\rho_{pol}$')
ax.legend()
fig.tight_layout()
plt.show()

