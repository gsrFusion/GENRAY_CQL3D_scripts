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

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.rc('axes', labelsize = 18)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 14)

machine = 'NTPT'

minSPA = 0.9

def makePlot():
    global NPara_fors
    global grillHeights
    global power
    fig,ax = plt.subplots(figsize=(7.5,4.8))

    fakeDevice = 'ARC'

    if fakeDevice == 'DIIID':
        shotNums = ['DIIID.147634PT', 'DIIID.147634NT', 'DIIID.193765PT', 'DIIID.193765NT']
        labels = ['PT 147634-like', 'NT 147634-like', 'PT 193765-like', 'NT 193765-like']
        colors = ['mediumblue', 'crimson' ,'darkturquoise', 'orange']
        ax.set_ylim([0,220])


    elif fakeDevice == 'ARC':
        shotNums = ['ARC.V3APT', 'ARC.V3ANT']
        labels = ['PT V3A-like', 'NT V3A-like']
        colors = ['teal', 'goldenrod']
        ax.set_ylim([0,325])
        #ax.set_ylim([0,1500])

    fig_curr, ax_curr = plt.subplots(figsize=(7.5,4.8))
    ax_curr.set_xlabel(r'LH current density centroid ($\rho_{pol}$)')
    ax_curr.set_ylabel(r'$\langle 1/ N_{||}^2 \rangle$')
    #ax_curr = None

    for k, shotNum in enumerate(shotNums):
        if k is not 0:
            pass
        if '147634' in shotNum:
            NPara_fors = -1*np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
            grillHeights = np.round(np.linspace(-.75,.75,13),3)
            power = 1
            
        elif '193765' in shotNum:
            NPara_fors = np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
            grillHeights = np.round(np.linspace(-.5,.5,11),3)
            power = 1

        elif 'V3A' in shotNum:
            grillHeights = np.round(np.linspace(-1.75,1.75,15),3)
            NPara_fors = -1*np.array([1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2])[:-1]
            power = 10
        addCurrentCurve(shotNum, ax, colors[k], labels[k], ax_curr)
    if ax_curr is not None:
        ax_curr.set_xlim([0,1])
        ax_curr.legend()
        fig_curr.tight_layout()
    #fig_curr.show()
    ax.set_ylabel('LH driven current (kA)')
    ax.set_xlabel(r'LH current density centroid ($\rho_{pol}$)')#(r'$\int \rho_{pol} \cdot J_{LH}\, dA \left/ \int J_{LH}\, dA \right.$')
    ax.set_xlim([0,1])
    ax.legend(loc = 'upper left')
    fig.tight_layout()
    #plt.savefig('toka_ARC_effic.jpeg',dpi=300)
    plt.show()

def addCurrentCurve(shotNum, ax, color, label, ax_curr):
    centroids = []
    totalCDs = []
    avgEffics = []

    #from omfit_classes import omfit_eqdsk
    from scipy.interpolate import interp1d

    Rs = None
    rho_pol_gfile = None

    for i in range(len(NPara_fors)):
        NPara_for = NPara_fors[i]
        prefix = 'n'
        if NPara_for > 0:
            prefix = 'p'

        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/{machine}_{shotNum}'
        stem2 = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}'
        
        if 'LFS' in label:
            stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shotNum}/LFSVersion/{machine}_{shotNum}'
        
        for j in range(len(grillHeights)):
            """
            if i == 0 and j == 0:
                eqdskName = f'{stem2}/gNTPT_{shotNum[:-2]}_'
                if 'NT' in shotNum:
                    eqdskName += 'NT_TM_10032026'
                else:
                    eqdskName += 'PT_TM_10032026'
                print(eqdskName)
                geqdsk = omfit_eqdsk.OMFITgeqdsk(eqdskName)
                Rs = geqdsk['fluxSurfaces']['avg']['R']
                rho_pol_gfile = np.sqrt(geqdsk['fluxSurfaces']['levels'])
            """
            grillHeight = grillHeights[j]
            targetDir = f'{stem}_{prefix}{np.abs(NPara_for)}Npara_{grillHeight}grillHeight_{power}MW'

            print(f'targetDir: {targetDir}')

            cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

            SPA = helper.getSPA(targetDir)[0]

            if SPA > minSPA:
               
                weightedAvg_Npara = helper.getAvgEfficMetricAtDamping(targetDir = targetDir, lobe = 1)

                curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
                rya = np.ma.getdata(cql_nc.variables["rya"][:])#convert to m^2

                darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
                totalCD = np.sum(curr*darea)
                print(f'totalCD: {totalCD}')
                loc, CD = helper.getAvgCurrentLocAndTotal(targetDir)
                #newRs = interp1d(rho_pol_gfile, Rs)(rya)
                #totalCD *= interp1d(rho_pol_gfile, Rs)(loc)

                centroids.append(loc)
                totalCDs.append(totalCD)
                avgEffics.append(weightedAvg_Npara)


    print(f'min, max: {np.nanmin(centroids), np.nanmax(centroids)}')

    centroids = np.array(centroids)
    totalCDs = np.array(totalCDs)
    avgEffics = np.array(avgEffics)
    totalCDs = totalCDs[np.argsort(centroids)]*1e3
    avgEffics = avgEffics[np.argsort(centroids)]
    centroids = np.sort(centroids)
    print(np.max(totalCDs))
    ax.scatter(centroids, totalCDs, marker = 'D', facecolors = 'none', edgecolors = color, s=75, label = label)

    if ax_curr is not None:
        ax_curr.scatter(centroids, avgEffics, facecolors = 'none', edgecolors = color, s=75, label = label)



makePlot()
