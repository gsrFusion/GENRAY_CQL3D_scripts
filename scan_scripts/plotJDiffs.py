##
# Plots the current profiles from a series of simulations
##


import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from scipy.interpolate import interp1d
import os, sys
import matplotlib
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getInputFileDictionary
import helperFunctions as helper

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 14)
plt.rc('legend', fontsize = 12)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()
import shotToEqdsk
#from omfit_classes import omfit_eqdsk

def plotSPAMatrix():
    machine = 'DIIID'


    if machine == 'DIIID':
        case = '203619 4130 E'
        if case == '203619 4130 E':
            shotNum = '203619.04130'

            stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff'

            targetDirs = [
                [f'{stem}',f'{stem}_noLH'],
                [f'{stem}_n0.025Fld',f'{stem}_n0.025Fld_noLH'],
                [f'{stem}_n0.05Fld',f'{stem}_n0.05Fld_noLH'],
                [f'{stem}_n0.075Fld',f'{stem}_n0.075Fld_noLH'],
            ]
            labels = [
                r'E = 0.0 V/m',
                r'E = -0.025 V/m',
                r'E = -0.05 V/m',
                r'E = -0.075 V/m',
            ]
            colors = ['tab:blue', 'tab:red', 'tab:green', "tab:purple", 'tab:cyan', 'tab:orange', 'tab:brown', 'y']
        else:
            pass
        fig,ax = plt.subplots()
        
        #fig,ax = plt.subplots(figsize=(5.7,5))
        
        print(f'starting making scans')

        for i,targetDirTuple in enumerate(targetDirs):
            cql_nc_yesLH = netCDF4.Dataset(f'{targetDirTuple[0]}/cql3d.nc','r')
            curr_yesLH = cql_nc_yesLH.variables["curr"][-1,:]*1e4/1e6#MA/m^2
            rya = np.ma.getdata(cql_nc_yesLH.variables["rya"][:])

            cql_nc_noLH = netCDF4.Dataset(f'{targetDirTuple[1]}/cql3d.nc','r')
            curr_noLH = cql_nc_noLH.variables["curr"][-1,:]*1e4/1e6#MA/m^2

            diff = curr_yesLH - curr_noLH

            ax.plot(rya, diff, lw = 3, label = labels[i], color = colors[i])#label = labels[i], color = colors[i])##r'N$_{||}$' + f' = {NPara_for}')

        ax.set_xlabel(r'$\rho_{pol}$')
        ax.set_ylabel(r'J(P$_{LH}$ = 54 kW) - J(P$_{LH}$ = 0 kW) ')
        ax.set_xlim([0,1])
        ax.legend(loc = 'best')
        fig.tight_layout()
        plt.show()

plotSPAMatrix()