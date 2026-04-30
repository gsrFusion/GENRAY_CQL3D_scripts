"""
Plots the ray traces and the RF power deposition density
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
    wj_emis_nc = genray_ece_nc.variables['wj_emis_nc'][:]
    wal_emis_nc = genray_ece_nc.variables['wal_emis_nc'][:]
    wsn_nc = genray_ece_nc.variables['wsn_nc'][:]
    wr_em_nc = genray_ece_nc.variables['wr_em_nc'][:]
    rho_pol = genray_ece_nc.variables['spsi'][:]
    wfreq_nc = genray_ece_nc.variables['wfreq_nc'][:]
    fig,ax = plt.subplots(figsize=(10,6))


    colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:cyan','tab:olive',
              'tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:cyan','tab:olive']
    linestyles = ['solid','solid','solid','solid','solid','solid','solid','solid','solid','solid',
                  'dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed',]

    for i in range(len(wj_emis_nc)):
        if i == 12 or i ==13:
            ax.plot(wr_em_nc[i], wj_emis_nc[i], label = f'f = {wfreq_nc[i]:.1f} GHz', color = colors[i], linestyle = linestyles[i])

    ax.set_xlabel('major radius ray')
    ax.set_ylabel('emission coefficient (erg*sec/cm**3)')
    ax.set_xlim([1.4,2.35])
    #ax.set_ylim([-.5e-9,4.45e-8])
    ax.set_ylim([1e-11,5e-7])
    ax.set_yscale('log')
    ax.legend(ncol = 2)
    fig.tight_layout()
    plt.show()


main()
