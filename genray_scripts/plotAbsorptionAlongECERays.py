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
    wal_emis_nc = genray_ece_nc.variables['wal_emis_nc'][:]*100 #convert to 1/m
    wsn_nc = genray_ece_nc.variables['wsn_nc'][:]/100 #convert to m
    wr_em_nc = genray_ece_nc.variables['wr_em_nc'][:]
    wfreq_nc = genray_ece_nc.variables['wfreq_nc'][:]
    fig,ax = plt.subplots()


    colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:cyan','tab:olive',
              'tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:cyan','tab:olive']
    linestyles = ['solid']*20#['solid','solid','solid','solid','solid','solid','solid','solid','solid','solid',
                 # 'dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed','dashed',]

    for i in range(len(wal_emis_nc)):
        if i == 0:#i == 12 or i ==13 or i == 14:
            ax.plot(wsn_nc[i][wr_em_nc[i] > 1], wal_emis_nc[i][wr_em_nc[i] > 1], lw = 3, label = f'f = {wfreq_nc[i]:.1f} GHz', color = colors[i], linestyle = linestyles[i])

    ax.set_xlabel('Distance along ray (m)')
    ax.set_ylabel('Absorption coefficient (1/m)')
    #ax.set_xlim([1.4,2.4])
    ax.set_ylim([1e-6,1e4])
    #ax.set_ylim([1e-11,5e-7])
    ax.set_yscale('log')
    ax.legend(ncol = 1)
    fig.tight_layout()
    plt.savefig('203619_absorption.jpeg',dpi=300)
    plt.show()


main()
