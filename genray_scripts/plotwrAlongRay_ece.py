"""
Plots the ray traces and the RF power deposition density
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os, sys
from scipy.integrate import cumtrapz

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getTargetInfo

targetDir = getTargetInfo.getTargetDir()
print(targetDir)
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('figure', titlesize = 22)
plt.rc('legend', fontsize = 14)

#plots either the toroidal and/or poloidal ray trajectories
def main():
    wal_emis_nc = genray_ece_nc.variables['wal_emis_nc'][:]
    wsn_nc = genray_ece_nc.variables['wsn_nc'][:]
    wr_em_nc = genray_ece_nc.variables['wr_em_nc'][:]
    wfreq_nc = genray_ece_nc.variables['wfreq_nc'][:]
    wr = genray_ece_nc.variables['wr'][:]*1e2
    ws = genray_ece_nc.variables['ws'][:]

    fig,ax = plt.subplots()

    for i in range(len(ws)):
        ax.plot(ws[i],wr[i])

    ax.set_xlabel('distance along ray')
    ax.set_ylabel(r'major radius (m)')
    #ax.set_xlim([0,170])
    #ax.set_ylim([1e-10,10])
    #ax.set_ylim([1e-11,5e-7])
    #ax.set_yscale('log')
    #ax.legend(ncol = 1)
    fig.tight_layout()
    plt.show()


main()
