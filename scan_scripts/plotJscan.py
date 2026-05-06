##
# Plots the current profiles from a series of simulations
##


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
import getTargetInfo

plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 14)
plt.rc('legend', fontsize = 12)

targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()
#from omfit_classes import omfit_eqdsk

def plotSPAMatrix():
    machine = 'DIIID'

    if machine == 'NTPT':
        shotNums = ['DIIID.193765PT','DIIID.193765NT',
                    ]
        
        stem = f'/home/grantr/symlinks/genray_batch/NTPT_shots/NTPT'

        targetDirs = []

        endingDirs = ['p2.9Npara_0.2grillHeight_1MW', 
                      'p2.9Npara_0.2grillHeight_1MW', 
                      ]
        colors = ['darkturquoise', 'orange']
        labels = [
            r'PT 193765-like', 
            r'NT 193765-like',

            ]

        fig,ax = plt.subplots()

        for i in range(len(shotNums)):
            shotNum = shotNums[i]
            targetDir = f'{stem}_{shotNum}/NTPT_{shotNum}_{endingDirs[i]}'
        
            cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
            curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#MA/m^2
            rya = np.ma.getdata(cql_nc.variables["rya"][:])
            darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
            totalCD = np.sum(curr*darea)
            print(f'width: {helper.getJcdWidth(targetDir)}, max: {np.max(curr)}')
            print(f'totalCD: {totalCD*1e3} kA')
            ax.plot(rya, curr, lw = 3, color = colors[i], label = labels[i])

        ax.set_xlim([0,1])
        ax.legend( loc = 'best')

        ax.set_ylabel(r'$J_{LH}$ $(MA/m^2)$')
        ax.set_xlabel(r'$\rho_{pol}$')
        ax.set_title(r'$N_{||,LCFS} = 2.9$, $Z_{launcher} = 0.2$')
        fig.tight_layout()
        plt.savefig('DIIID_twoColor.jpeg', dpi=300)
        plt.show()

    if machine == 'DIIID':
        case = '147634 2 color'

        if case == '147634 2 color':
            shotNum = '147634.04525'

            stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/twoColorScans/{machine}_{shotNum}'
            shotNums = [shotNum]*3

            targetDirs = [
                f'{stem}_n2.7Npara_n2.7Npara_1MW',
                f'{stem}_n2.7Npara_n3.0Npara_1MW',
                f'{stem}_n3.0Npara_n3.0Npara_1MW',
            ]
            labels = [
                r'$N_{||} = -2.7$',
                r'$N_{||} = -2.7, -3.0$',
                r'$N_{||} = -3.0$',
            ]
            colors = ['tab:blue', 'tab:green','tab:red', 'tab:green', "tab:red", 'tab:cyan', 'tab:grey', 'tab:brown', 'y', 'tab:pink']


        if case == '4130 Zeff scan':
            shotNum = '203619.04130'
            stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
            stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_1.6Zeff/DIIID_{shotNum}_expSpectrum_1.6Zeff'
            stem3 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2.4Zeff/DIIID_{shotNum}_expSpectrum_2.4Zeff'
            stem4 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff_id2/DIIID_{shotNum}_expSpectrum_2Zeff_id2'

            targetDirs = [
                f'{stem1}_second',
                f'{stem2}_second',
                f'{stem3}_second',
                f'{stem4}_second',
            ]
            labels = [
                r'Z$_{eff}$ = 1.6',
                r'Z$_{eff}$ = 2',
                r'Z$_{eff}$ = 2.4',
                r'Z$_{eff}$ = 2, id = 2',

            ]
            colors = [ 'tab:red', 'tab:blue','seagreen', 'tab:purple', 'tab:grey']

        else:
            pass
        fig,ax = plt.subplots()
        
        #fig,ax = plt.subplots(figsize=(5.7,5))
        
        print(f'starting making scans')

        for i,targetDir in enumerate(targetDirs):
            print(targetDir)
            try:
                SPA,_,_ = helper.getSPA(targetDir)
                print(f'SPA: {SPA}')

                if SPA[0] > 0:
                    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
                    curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#MA/m^2
                    rya = np.ma.getdata(cql_nc.variables["rya"][:])
                    darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
                    totalCD = np.sum(curr*darea)
                    print(f'{totalCD*1000} kA')

                    label = labels[i] + f'\nI$_{{LH}}$ = {totalCD*1e3:.2f} kA'

                    ax.plot(rya, curr, lw = 3, label = label, color = colors[i])
            except Exception as e:
                print(e)
                pass

        ax.set_xlabel(r'Minor radius ($\rho_{pol}$)')#ax.set_xlabel(r'$\rho_{pol}$')
        ax.set_ylabel(r'LH-driven current density (MA/$m^2$)')
        ax.set_title(f'Shot {shotNum.split(".")[0]}, {shotNum.split(".")[1][1:]} ms', loc = 'right')
        ax.set_xlim([0,1])
        ax.legend(loc = 'best')
        fig.tight_layout()
        #plt.savefig('DIIID_twoColor.jpeg', dpi=300)

        plt.show()

plotSPAMatrix()