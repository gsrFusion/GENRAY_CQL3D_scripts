import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from scipy.interpolate import interp1d
import os, sys
import matplotlib
import pickle
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

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()
from omfit_classes import omfit_eqdsk
import shotToEqdsk
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def multi_gauss(x, *params):
    y = np.zeros(len(x))
    num_gaussians = len(params)//3
    for i in range(num_gaussians):
        amp = params[i*3]
        cen = params[i*3+1]
        sigma = params[i*3+2]
        y += amp*np.exp(-(x-cen)**2/(2*sigma**2))

    return y

def writeInformation():
    #NPara_fors = -1*np.array([2.67,2.8,2.9,3.0])
    
    time = '04400'
    shot = '180403'
    machine = 'DIIID'
    targetSuffix = '_1MW'#_0.2Width_0.1deltaT'

    maxAllowableGauss = 1

    JPeakCenters = []
    JPeakMags = []
    JPeakFWHMs = []
    JPeakNparas = []

    if shot == '184281':
        if targetSuffix == '_2MW':
            NPara_fors = -1*np.array([2.95,2.97,3.0,3.05,3.1,3.15,3.2])
        elif targetSuffix == '_2MW_0.1Width_0.1deltaT':
            NPara_fors = -1*np.array([2.95,2.97,3.0,3.05,3.1,3.2])
        elif targetSuffix == '_2MW_0.075Width_0.1deltaT':
            NPara_fors = -1*np.array([2.95,2.97,3.0,3.05,3.1,3.15,3.2])
        elif targetSuffix == '_2MW_0.05Width_0.1deltaT':
            NPara_fors = -1*np.array([2.95,2.97,3.0,3.05,3.1,3.15,3.2])
        else:
            NPara_fors = -1*np.array([2.95,2.97,3.0,3.05,3.1,3.15,3.2])


    if shot == '199605':
        if targetSuffix == '_2MW_0.1Width_0.1deltaT':
            NPara_fors = -1*np.array([2.95,3.0,3.1,3.15,3.2,3.25,3.3,3.4])
        if '_1MW' in targetSuffix:
            NPara_fors = -1*np.array([3.1,3.15,3.2,3.25,3.3])
        if '_2MW' == targetSuffix:
            NPara_fors = -1*np.array([3.1,3.15,3.2,3.25,3.3])
        else:
            NPara_fors = -1*np.array([2.9,2.95,3.0,3.1,3.15,3.2,3.25,3.3,3.4])



    if shot == '180403':
        if targetSuffix == '_1MW':
            #NPara_fors = -1*np.array([2.8,2.89,2.9,2.91,2.92,2.93,2.95,3])
            NPara_fors = -1*np.array([2.8,2.89,2.91,2.93,2.95,])
        else:
            NPara_fors = -1*np.array([2.8,2.85,2.9,2.91,2.92,2.93,2.95,3])


    if shot == '199375':
        if targetSuffix == '_1MW':
            NPara_fors = np.array([2.3,2.5,2.7,2.9,3.1])
        else:
            NPara_fors = -1*np.array([2.8,2.85,2.9,2.91,2.92,2.93,2.95,3])


    topshotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}.{time}'
    particularShotStem = f'{topshotDir}/{machine}_{shot}.{time}_'

    eqdskName = shotToEqdsk.getEqdskName(f'{shot}.{time}', machine = machine)

    gfile = omfit_eqdsk.OMFITgeqdsk(f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}.{time}/{eqdskName}')

    omfit_rho_pol = np.sqrt(gfile['fluxSurfaces']['levels'])
    rho_interp = np.linspace(0,1,1000)
    rho_pol = rho_interp
    #rho_pol = np.linspace(.01,.99,300)
    #rho_pol = np.load(f'/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/NTM_rho_pol.npy')

    R_LFSmidplane = helper.convertRhopolToRmidplane(rho_interp, targetDir = topshotDir, side = 'LFS')
    #R_LFSmidplane = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['a'])(rho_interp)
    #ws = np.linspace(0.01, .3, 100)

    ws = np.load(f'/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/NTM_ws_m.npy')

    q2FluxWidths = np.zeros(len(ws))

    q_prof = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['q'])(rho_pol)

    q2Index = helper.findNearestIndex(2, np.abs(q_prof))
    print(f'q=2 occurs at rho_pol = {rho_pol[q2Index]}')
    LFSmidplaneR = helper.convertRhopolToRmidplane(rho_pol)
    R_q2 = LFSmidplaneR[q2Index]

    for l in range(len(ws)):
        leftHalfFluxWidthIndex = helper.findNearestIndex(-ws[l]/2, LFSmidplaneR-R_q2)
        rightHalfFluxWidthIndex = helper.findNearestIndex(ws[l]/2, LFSmidplaneR-R_q2)
        q2FluxWidths[l] = rho_pol[rightHalfFluxWidthIndex]**2 - rho_pol[leftHalfFluxWidthIndex]**2
    """
    fig,ax = plt.subplots()
    ax.plot(ws, q2FluxWidths)
    ax.set_ylabel('q2 flux width')
    ax.set_xlabel('w (m)')
    plt.show()
    """
    fig,ax = plt.subplots(figsize = (7.5,4.8))

    for i in range(len(NPara_fors)):
        NPara_for = NPara_fors[i]
        prefix = 'n'
        if NPara_for > 0:
            prefix = 'p'
        targetDir = f'{particularShotStem}{prefix}{np.abs(NPara_for)}Npara{targetSuffix}'
        print(f'targetDir: {targetDir}')

        cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
        curr = cql_nc.variables["curr"][-1,:]*1e4#convert to A/m^2
        rya = np.ma.getdata(cql_nc.variables["rya"][:])
        
        SPA = helper.getSPA(targetDir)[0]
        print(f'SPA: {SPA}')
        if SPA > 0.9:
            percentError = np.inf
            numGauss = 1

            peaks,_ = find_peaks(curr, distance = 3, height = np.max(curr)*0.2, width = None)
            peakHeights = curr[peaks]
            sortedIndices = np.argsort(peakHeights)
            topPeakHeights = np.flip(peakHeights[sortedIndices])
            topPeakIndices = np.flip(peaks[sortedIndices])

            p0 = []
            bounds = []
            while numGauss <= maxAllowableGauss and np.abs(percentError) > .26:# and numGauss <= len(topPeakHeights):
                
                #"""
                if numGauss > len(topPeakHeights):
                    p0.extend([topPeakHeights[0], rya[topPeakIndices[0]], .02])
                    bounds.extend([[topPeakHeights[0]*.2,topPeakHeights[0]*1.1],
                                [0,1],
                                [.016, .5]])
                else:
                    p0.extend([topPeakHeights[numGauss-1], rya[topPeakIndices[numGauss-1]], .02])
                    bounds.extend([[topPeakHeights[numGauss-1]*.7,topPeakHeights[numGauss-1]*1.1],
                                [0,1],
                                [.016, .5]])
                popt,pcov = curve_fit(multi_gauss, rya, curr, p0=p0, bounds = np.array(bounds).T) 

                errorMask = np.where(curr > 1e3)
                percentError = (multi_gauss(rya[errorMask], *popt)-curr[errorMask]) / curr[errorMask]
                percentError = np.sum(percentError)/len(percentError)
                print(f'percent error: {percentError}')
                #"""
                numGauss += 1
            numGauss-=1#to account for the addition at the end of the final loop
            print(f'numGauss: {numGauss}')
            print(f'popt: {popt}')
            if numGauss == 1:
                j = 0
            else:
                mags = popt[0::3]
                print(f'mags: {mags}')
                maxIndex = np.argmax(mags)
                j = maxIndex

            JPeakNparas.append(NPara_for)

            JPeakMags.append(popt[j*3])
            JPeakCenters.append(popt[j*3+1])
            singleFit = multi_gauss(rho_interp, *popt[j*3:j*3+3])
            mask = np.where(singleFit > 0.5 * popt[j*3])
            #ax.plot([rho_interp[mask][0], rho_interp[mask][-1]], [0.5 * popt[j*3],0.5 * popt[j*3]])
            
            width = R_LFSmidplane[mask][-1]-R_LFSmidplane[mask][0] 
            JPeakFWHMs.append(width)

            gaussFit = multi_gauss(rho_interp, *popt)

            #ax.plot(rya, 4e5*i+curr, lw = 2, label = r'N$_{||}$' + f' = {NPara_for}')
            ax.plot(rya, curr/1e6, lw = 3, label = r'N$_{||}$' + f' = {NPara_for}', zorder = 100)
            #ax.scatter(rho_interp, gaussFit, lw = 2, linestyle = 'dashed', marker = 'D')

            ax.set_xlim([.55,.9])
            ax.set_ylim([-1e4/1e6,4.5e5/1e6])

            

    ax.set_xlabel(r'$\rho_{pol}$')
    ax.set_ylabel(r'$   J_{LH} (MA/m^2$)')
    ax.axvline(rho_pol[q2Index], lw = 2, color = 'k', linestyle = 'dashed', label = 'q=2')

    ax.legend()
    #ax.set_ylim([0,500])
    ax.set_title(f'Shot {shot}, 1 MW')
    fig.tight_layout()
    #plt.show()

    widthFunc = interp1d(JPeakCenters, JPeakFWHMs, kind = 'linear', bounds_error = False, fill_value = np.inf)
    magFunc = interp1d(JPeakCenters, JPeakMags, kind = 'linear', bounds_error = False, fill_value = 0)
    widths = widthFunc(rho_pol)
    print(f'widthFunc(.74) : {widthFunc(.74)}')
    peaks = magFunc(rho_pol)
    print(f'JPeakCenters: {JPeakCenters}')
    print(f'JPeakFWHMs: {JPeakFWHMs}, JPeakMags: {JPeakMags}')
    print(f'width at q2, mag at q2: {widthFunc(rho_pol[q2Index]), magFunc(rho_pol[q2Index])}')

    plt.savefig('NTM_locs.jpeg',dpi=300)


    fig,ax = plt.subplots()
    ax.plot(rho_pol, peaks)
    ax.scatter(rho_pol, np.zeros(len(rho_pol)))
    ax.scatter(JPeakCenters, JPeakMags)
    ax.set_ylabel('J_LH (A/m^2)')
    ax.set_xlabel('rho_pol')
    ax.axvline(rho_pol[q2Index], lw = 2, color = 'k', linestyle = 'dashed', label = 'q=2')
    ax.scatter([rho_pol[q2Index]], [magFunc(rho_pol[q2Index])])
    #plt.show()
    
    fig,ax = plt.subplots()
    ax.plot(rho_pol, widths)
    ax.scatter(rho_pol, np.zeros(len(rho_pol)))
    ax.scatter(JPeakCenters, JPeakFWHMs)
    ax.set_ylabel('J_LH width (m)')
    ax.set_xlabel('rho_pol')
    ax.axvline(rho_pol[q2Index], lw = 2, color = 'k', linestyle = 'dashed', label = 'q=2')

    ax.scatter([rho_pol[q2Index]], [widthFunc(rho_pol[q2Index])])
    plt.show()

    resultDict = {'q2FluxWidths_psiN': q2FluxWidths,
                  'JPeakCenters_rhop' : JPeakCenters,
                  'JPeakFWHMfits_m': JPeakFWHMs,
                  'JPeakMagsfits_Aperm2': JPeakMags,
                  'JPeakNparas' : JPeakNparas
                  #'JFWHMs_m' : widths,
                  #'JMags_Aperm2' : peaks,
                  }

    storageStem = f'/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/{machine}_{shot}.{time}/'

    """
    with open(f'{storageStem}{machine}_{shot}.{time}{targetSuffix}_resultDict.pkl', 'wb') as handle:
        pickle.dump(resultDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    np.save(f'/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/NTM_rho_pol.npy',rho_pol)
    #np.save(f'/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/NTM_ws_m.npy',ws)
    """
writeInformation()