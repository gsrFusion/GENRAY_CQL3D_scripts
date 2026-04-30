###
# Plot the ray trajectories and damping as predicted by GENRAY
# since it's a linear code, the damping is very wrong for LH, but the ray trajectory can be useful
###


import numpy as np
import matplotlib.pyplot as plt

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

import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
machine = getTargetInfo.getMachine()
genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')#netCDF4.Dataset(f'{targetDir}/genray.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 14)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)


#adds the ray traces to ax
def main():
    
    fig,ax = plt.subplots()

    ws  = genray_ece_nc.variables['ws'][:]/100 #major radius of the ray at each point along the trace, in m
    radialVariable = (np.copy(genray_ece_nc.variables["spsi"])) #rho_pol of the ray at each point along the ray trace

    for ray in range(len(ws)):
        #delpwr[ray,:] = delpwr[ray,:]/delpwr[ray,0] #normalize the ray power to that ray's starting power
        ax.plot(ws[ray],radialVariable[ray],lw = 3)

    ax.set_ylabel(r'$\rho_{pol}$')
    ax.set_xlabel(r'Poloidal distance along ray (m)')
    fig.tight_layout()
    plt.show()

main()
