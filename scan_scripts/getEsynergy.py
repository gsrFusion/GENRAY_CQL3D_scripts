import numpy as np
import matplotlib.pyplot as plt
import netCDF4
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
plt.rc('legend', fontsize = 14)

def getEsynergy():
    time = '04000'
    shot = '200388'
    machine = 'DIIID'

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}.{time}/3modules/{machine}_{shot}.{time}_'

    endings = ['Efield', 'p2.7Npara_300kW_3modules_Efield', 'p2.7Npara_300kW_3modules']

    J = {}

    print(f'starting making scans')
    for i in range(len(endings)):
        ending = endings[i]
        targetDir = f'{stem}{ending}'
        print(f'targetDir: {targetDir}')

        cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
        curr = cql_nc.variables["curr"][-1,:]
        rya = cql_nc.variables["rya"][:]
        J[ending] = curr


    fig,ax = plt.subplots()
    #ax.grid()
    ax.plot(rya, J['Efield'] + J['p2.7Npara_300kW_3modules'], lw=2, label = 'RF + Efield after CQL3D')
    ax.plot(rya, J['p2.7Npara_300kW_3modules'], lw=2, label = 'RF')
    ax.plot(rya, J['p2.7Npara_300kW_3modules_Efield'], lw=2, label = 'RF + Efield in CQL3D')
    ax.plot(rya, J['p2.7Npara_300kW_3modules_Efield'] - J['Efield'], lw=2, label = 'diff')
    ax.plot(rya, J['Efield'], lw=2, label = 'E field')
    ax.set_xlabel(r'$\rho_{pol}$')
    ax.set_ylabel(r'J (A/cm^2)')
    ax.legend()
    ax.set_title(f'Shot {shot}.{time}')
    fig.tight_layout()
    plt.show()

getEsynergy()