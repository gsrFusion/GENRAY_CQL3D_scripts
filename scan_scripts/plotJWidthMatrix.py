"""
Plots a matrix where the x and y axes designate a GENRAY/CQL3D run and the colorbar is the FWHM of the current density profile
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

import helperFunctions as helper

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.rc('axes', labelsize = 18)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 14)

fakeShot = '.193765NT'

if '193765' in fakeShot:
    fakeMachine = 'DIIID'
    NPara_fors = np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
    grillHeights = np.round(np.linspace(-.5,.5,11),3)
    power = 1

if '147634' in fakeShot:
    fakeMachine = 'DIIID'
    NPara_fors = -1*np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
    grillHeights = np.round(np.linspace(-.75,.75,13),3)
    power = 1

if 'V3A' in fakeShot:
    fakeMachine = 'ARC'
    grillHeights = np.round(np.linspace(-1.75,1.75,15),3)
    NPara_fors = -1*np.array([1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2])
    power = 10

machine = 'NTPT'

def plotJFWHMMatrix():
    JFWHMMatrix = np.zeros((len(NPara_fors), len(grillHeights)))

    minSPA = 0.9

    avgWidth = 0
    counter = 0

    R_lfs = None
        
    for i in range(len(NPara_fors)):
        NPara_for = NPara_fors[i]
        prefix = 'n'
        if np.sign(NPara_for) > 0:
            prefix = 'p'
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeMachine}{fakeShot}/{machine}_{fakeMachine}{fakeShot}'
        for j in range(len(grillHeights)):
            targetDir = targetDir = f'{stem}_{prefix}{np.abs(NPara_for):.1f}Npara_{grillHeights[j]}grillHeight_{power}MW'
            
            if R_lfs is None:
                cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
                rya = np.ma.getdata(cql_nc.variables["rya"][:])
                R_lfs = helper.convertRhopolToRmidplane(rya, targetDir, side = 'LFS')

            print(f'targetDir: {targetDir}')
            SPA = helper.getSPA(targetDir)[0]
            print(f'SPA: {SPA}')
            if SPA > minSPA:
                JFWHMMatrix[i,j] = helper.getJcdWidth(targetDir, R_lfs = R_lfs)
                counter += 1
                avgWidth += JFWHMMatrix[i,j]
            else:
                JFWHMMatrix[i,j] = np.nan
            print(f'width = {JFWHMMatrix[i,j]}')

    print(f'min FWHM: {np.nanmin(JFWHMMatrix)}')
    print(f'mean FWHM: {np.nanmean(JFWHMMatrix)}')
    print(f'max FWHM: {np.nanmax(JFWHMMatrix)}')
    fig,ax = plt.subplots()
    p=ax.pcolormesh(grillHeights, NPara_fors, JFWHMMatrix ,shading = 'nearest',cmap='viridis')
    cbar = fig.colorbar(p, ax = ax, shrink = .9, pad = .01)
    cbar.set_label(r'FWHM of J$_{LH}$ ($\rho_{pol}$)')
    
    ax.set_ylabel(r'$N_{||,LCFS}$')
    ax.set_xlabel(r'Grill vertical location (m)')

    ax.set_title(rf'{machine} {fakeMachine}{fakeShot}, P$_{{LH, for}}$ = 1 MW')

    fig.tight_layout()
    plt.show()

plotJFWHMMatrix()