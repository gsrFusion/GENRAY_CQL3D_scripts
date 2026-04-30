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
plt.rc('legend', fontsize = 16)

import numpy as np
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import netCDF4
print(f'past imports')

machine = 'NTPT'

if machine == 'NTPT':
    time = '.147634NT'
    shot = 'DIIID'
    if '193765' in time:
        grillHeights = np.round(np.linspace(-.5,.5,11),3)
        NPara_targets = np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
        power = 1

    elif '147634' in time:
        grillHeights = np.round(np.linspace(-.75,.75,13),3)
        NPara_targets = -1*np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])
        power = 1

    elif 'V3A' in time:
        grillHeights = np.round(np.linspace(-1.75,1.75,15),3)
        NPara_targets = -1*np.array([1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2])

        power = 10


    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{machine}_{shot}{time}'

m_e = 9.109e-31
m_D = 3.343e-27
eps_0 = 8.854e-12
e = 1.6e-19
c = 3e8
w = 2*np.pi*4.6e9

def getSDP(n, B):
    w_ce = -e*B/m_e
    w_cD = e*B/m_D
    w_pe2 = e**2*n/(m_e*eps_0)
    w_pD2 = e**2*n/(m_D*eps_0)

    S = 1 - (w_pe2/(w**2 - w_ce**2)) - (w_pD2/(w**2 - w_cD**2))
    D = (w_ce/w)*(w_pe2/(w**2 - w_ce**2)) + (w_cD/w)*(w_pD2/(w**2 - w_cD**2))
    P = 1 - w_pe2/w**2 - w_pD2/w**2
    return S, D, P

def doesModeConv(targetDir):
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')
    cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

    nparas = np.copy(genray_nc.variables["wnpar"]) #n_|| of the ray at each point along the ray trace
    nperps = np.copy(genray_nc.variables["wnper"]) #n_perp of the ray at each point along the ray trace
    wr  = genray_nc.variables["wr"][:]/100 #major radius of the ray at each point along the trace
    wr = np.ma.getdata(wr)
    wz  = genray_nc.variables["wz"][:]/100 #height of the ray at each point along the trace
    wz = np.ma.getdata(wz)
    ws  = genray_nc.variables["ws"][:]/100 #poloidal distance along ray
    delpwr= np.copy(cqlrf_nc.variables["delpwr"]) #power in the ray at each point

    #dielectric tensor components along the ray
    cweps11 = genray_nc.variables["cweps11"][:]
    cweps11 = cweps11[0] + 1j*cweps11[1]
    cweps12 = genray_nc.variables["cweps12"][:]
    cweps12 = cweps12[0] + 1j*cweps12[1]
    cweps21 = genray_nc.variables["cweps21"][:]
    cweps21 = cweps21[0] + 1j*cweps21[1]
    cweps33 = genray_nc.variables["cweps33"][:]
    cweps33 = cweps33[0] + 1j*cweps33[1]

    minRatioToPlot = .05#.1 #when the ratio of ray power to ray starting power is below this number, the trace ends
    summed = 0
    for i in range(len(nparas)):
        npara = nparas[i]
        nperp = nperps[i]
        N = np.sqrt(npara**2 + nperp**2)

        delpwrRatios = delpwr[i]/np.max(delpwr[i])
        mostPowerDep = helper.findNearestIndex(minRatioToPlot, delpwrRatios)

        S, D, P = np.real(cweps11[i]), np.imag(cweps21[i]), np.real(cweps33[i])
        
        #"""
        A = S
        B = -((S - npara**2)*(P + S) - D**2)
        C = P*((S - npara**2)**2 - D**2)
        Delta = B**2-4*A*C

        #fast wave solution for N_perp^2
        fastNperp2 = (-B - np.sqrt(Delta))/(2*A)
        #slow wave solution for N_perp^2
        slowNperp2 = (-B + np.sqrt(Delta))/(2*A)
        
        fastNperp = np.sqrt(fastNperp2)
        fastNperp[fastNperp2 < 0] = np.inf
        slowNperp = np.sqrt(slowNperp2)
        slowNperp[slowNperp2 < 0] = np.inf
        
        distToFast = np.abs(np.abs(nperp)- fastNperp)
        distToSlow = np.abs(np.abs(nperp) - slowNperp)
        
        slowOrFast = np.zeros(len(npara))#zero where the wave is slow
        slowOrFast[distToFast < distToSlow] = 1

        summed += np.sum(slowOrFast)

    return summed


def plotIfModeConv():
    convMatrix = np.zeros((len(NPara_targets), len(grillHeights)))
    counter_conv = 0
    counter_total = 0
    for i in range(len(NPara_targets)):
        NPara_for = NPara_targets[i]
        prefix = 'n'
        if np.sign(NPara_for) > 0:
            prefix = 'p'
        for j in range(len(grillHeights)):
            targetDir = f'{stem}_{prefix}{np.abs(NPara_for)}Npara_{grillHeights[j]}grillHeight_{power}MW'

            try:
                print(f'targetDir: {targetDir}')
                SPA = helper.getSPA(targetDir)[0]
                print(f'SPA = {SPA}')
                print(f'{SPA >= 0.9}')

                if SPA >= 0.9:
                    summed = doesModeConv(targetDir)
                    counter_total+=1
                    if summed > 0:
                        convMatrix[i,j] = 1
                        counter_conv += 1
                else:
                    convMatrix[i,j] = np.nan
            except Exception as e:
                import traceback
                print(traceback.format_exc())
                convMatrix[i,j] = np.nan

    print(f'number of cases mode converting: {counter_conv}, total cases: {counter_total}')
    fig,ax = plt.subplots(figsize=(7.5,5.5))
    p = ax.pcolormesh(grillHeights, NPara_targets, convMatrix,shading = 'nearest', 
                       cmap='inferno_r', vmin=0, vmax = 1)
    
    print(f'min loc at: {np.nanmin(convMatrix)}')

    ax.set_ylabel(r'N$_{||,LCFS}$')
    ax.set_xlabel(r'$Z_{launcher}$ (m)')
    ax.set_yticks(NPara_targets)
    ax.set_xticks(grillHeights[::2])
    ax.tick_params(axis='x', rotation=30)
    ax.yaxis.get_label().set_fontsize(16)

    triString = ''
    if 'PT' in time:
        triString = 'Positive'
    else:
        triString = 'Negative'
    ax.set_title(f'{triString} triangularity {time[1:-2]}-like')
    #ax.set_title(r'Positive triangularity V3A-like ($\delta = 0.5$)')

    cbar = fig.colorbar(p, ax = ax, shrink = .9, pad = .01)
    #cbar.set_label(r'$\int \rho_{pol} \cdot J_{LH}\, dA \left/ \int J_{LH}\, dA \right.$')
    cbar.set_label(r'Current centroid ($\rho_{pol}$)')
    

    fig.tight_layout()
    if 'PT' in time:
        #plt.savefig('toka_V3A_PT_depLoc.jpeg',dpi=300)
        pass
    if 'NT' in time:
        #plt.savefig('toka_V3A_NT_depLoc.jpeg',dpi=300)
        pass

    plt.show()

plotIfModeConv()