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

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

def plotSPAMatrix():
    NPara_fors = -1*np.array([2.5,2.7,2.9])
    time = '04565'
    shot = '147634'
    machine = 'DIIID'
    I = np.zeros(len(NPara_fors))
    SPAs = np.zeros(len(NPara_fors))
    print(f'starting making scans')

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}.{time}/Npara_thgrill_scan/{machine}_{shot}.{time}_'
    print(f'starting making scans')
    for i in range(len(NPara_fors)):
        prefix = 'n'
        if NPara_fors[i] > 0:
            prefix = 'p'
        targetDir = f'{stem}{prefix}{np.abs(NPara_fors[i])}Npara_180thgrill_1MW'
        print(f'targetDir: {targetDir}')

        cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
        curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
        darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
        totalCD = np.sum(curr*darea)
        
        I[i] = totalCD

        SPA_allLobes,_,_ = helper.getSPA(targetDir=targetDir)
        SPAs[i] = SPA_allLobes[0]

    #SPAmatrix[SPAmatrix < 0.95] = np.NAN

    fig,ax = plt.subplots()
    #ax.grid()
    ax.plot(NPara_fors, SPAs, lw=2, color = 'tab:blue')
    ax.yaxis.label.set_color('tab:blue')
    ax.set_xlabel(r'N$_{||,forward}$')
    ax.set_ylabel(r'SPA')

    ax2 = ax.twinx()
    ax2.plot(NPara_fors, I*1000,lw=2, color = 'tab:orange')
    ax2.yaxis.label.set_color('tab:orange')
    ax2.set_ylabel(r'I$_{LH}$ (kA)')       

    ax.set_xticks(NPara_fors)
    ax.set_xticklabels(ax.get_xticks(), rotation = 30)
    ax.set_ylim([0,1.01])

    ax.set_title(f'Shot {shot}')
    fig.tight_layout()
    plt.show()

plotSPAMatrix()