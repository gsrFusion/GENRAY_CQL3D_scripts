import numpy as np
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import setupInputFiles
import generateNparaSpectrum
import netCDF4
import helperFunctions
print(f'past imports')

time = '.02700'
shot = '203912'
machine = 'DIIID'
whatCode = 'test'

calcThgrill = False

if machine == 'DIIID':
    #these module list for ordered module 1 to module 8
    if shot == 'test':
        modulePhases_DIIID = np.array([260.0, 173.0, 87.0, 0.0, 273.0, 187.0, 100.0, 14.0])
        modulePowers = np.array([1, 1, 1, 1, 1, 1, 1, 1])
        #isNormalBt = True

    if shot == 'phaseTest':
        modulePhases_DIIID = np.array([0,0,87.0, 0.0, 273.0, 187.0, 100.0, 14.0])
        modulePowers = np.array([0,0, 1, 1, 0,0, 1, 0])

    if shot == '203619' and '.04130' in time:
        modulePhases_DIIID = np.array([0, 0, 0, 0, 4.5, 177.13, 0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 59.18,48.54, 0.0, 0.0])
        calcThgrill = True
        isNormalBt = False

    if shot == '203912' and time == '.02700':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.6, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = True

    if shot == '203917' and time == '.03700':
        modulePhases_DIIID = np.array([0, 0, 43.9, 0.0, 327.6, 0, 77.0, 0])
        modulePowers = np.array([0.0, 0.0, 8.49, 27.0, 38.6, 0.0, 41.0, 0.0])
        isNormalBt = True
        
    if shot == '203917' and time == '.02800':
        modulePhases_DIIID = np.array([0, 0, 47.6, 0.0, 328.2, 0, 74.833, 0])
        modulePowers = np.array([0.0, 0.0, 7.0, 25.6, 36.0, 0.0, 38.0, 0.0])
        isNormalBt = True

    if shot == '203917' and time == '.03000':
        modulePhases_DIIID = np.array([0, 0, 47.6, 0.0, 328.2, 0, 74.833, 0])
        modulePowers = np.array([0.0, 0.0, 7.0, 25.6, 36.0, 0.0, 38.0, 0.0])
        isNormalBt = True

    if shot == '206629' and time == '.01980':
        modulePhases_DIIID = np.array([0, 0, 87.5, 0.0, 0, 0, 96.5, 0])
        modulePowers = np.array([0.0, 0.0, 43.11, 40.26, 0.0, 0.0, 36.6, 0.0])
        isNormalBt = True

    if shot == '206636' and time == '.01960':
        modulePhases_DIIID = np.array([0, 0, 87.5, 0.0, 0, 0, 96.5, 0])
        modulePowers = np.array([0.0, 0.0, 34.9, 44.94, 0.0, 0.0, 0, 0.0])
        isNormalBt = True

    modulePowers_normed = modulePowers/np.max(modulePowers)
    print(f'modulePowers: {modulePowers.tolist()}')
    modulePhases_DIIID%=360
    modulePhases_Andrew = modulePhases_DIIID
    modulePhases_Andrew %= 360

    modulePhases_Andrew_rad = np.pi*modulePhases_Andrew/180

    #from the phasing, get the spectrum
    peakNparas, peakEdges, directivities = generateNparaSpectrum.generateSpectrum(target_npara=None, #if not supplying the phase shift, what is your target N||
                    modulePhaseShift = modulePhases_Andrew_rad, #RADIANS, if not supplying the target N||, what is the phase shift
                    analytic = False, #if you want to calculate things analytically. Only valid for equal powers in all WGs
                    modulePowerRatio = modulePowers_normed, #power ratios between in the modules
                    #wgPowerRatio = None, # power ratios in the WGs within each module
                    numLobes = 3,
                    doPlot = 'spectrum',) #whether or not to plot
    
    #may be useful if you want to ignore the reverse lobe
    """
    peakNparas = np.array([peakNparas[0], peakNparas[2], peakNparas[3]])
    peakEdges = np.array([peakEdges[0], peakEdges[2], peakEdges[3]])
    directivities = np.array([directivities[0], directivities[2], directivities[3]])
    """
    print(peakNparas)
    print(peakEdges)
    print(directivities)

    #generateSpectrum assumes rev Bt. If normal Bt, then we need to flip the sign of N||
    if isNormalBt:
        peakNparas = -peakNparas
        peakEdges = -1*np.flip(peakEdges,axis = 1)

    totalPower_kW = np.sum(modulePowers)
    pwrFactor = 1
    #This factor of 0.72 is an estimate of the resistive losses between the klystrons and the plasmas
    #Comes from table 6 of https://doi.org/10.1016/j.fusengdes.2020.111762
    pwrscale = 0.72*totalPower_kW/1000*pwrFactor
    print(f'starting input file helper')
    intermediateDir = ''#'60nnkpar_1e-8prmt4_1e-7prmt4ECE/'
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{intermediateDir}{machine}_{shot}{time}'

    targetDir = f'{stem}_expSpectrum_RaymondTest'

    """
    k = 1
    while k < len(directivities):
        print(k)
        if directivities[k]/directivities[0] < .10:
            print(directivities[k]/directivities[0] )
            break
        k += 1
    print(f'k: {k}, len direct: {len(directivities)}')
    directivities = directivities[:k]
    N_para_peaks = peakNparas[:k]
    N_para_edges = peakEdges[:k]
    """    

    N_para_peaks = peakNparas
    N_para_edges = peakEdges
    powerInLobes = directivities*1e6#for 1 MW of forward power
        
    print(f'N_para_peaks: {N_para_peaks}')

    thgrill = 189
    """
    if the z axis is not located at (or very near) Z=0, thgrill needs to be recalculated
    """
    if calcThgrill:
        print(f'yes')
        thgrill = np.inf
    #"""
    helper = setupInputFiles.InputFileHelper(targetDir,  
        waveType = 'LH',
        makeDir = True, overwrite = True, doPlot = True,
        nScale = 1, TScale = 1, ZeffScale = 1,  
        numCQLToFokkerPlanck = 50, ndens = 101, njene= 101, 
        isScoping = False, eqsym = 'average',
        thgrill=thgrill, powerInLobes = powerInLobes,  N_para_edges = N_para_edges, 
        pwrScale = pwrscale, N_para_peaks = N_para_peaks,
    )
    helper.copySetupAndClean() 
    if whatCode == 'CQL3D':
        os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runCQL.sh {targetDir}')
    elif whatCode == 'both':
        os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGENThenCQL.sh {targetDir}')
    elif whatCode == 'GENRAY':
        os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGEN.sh {targetDir}')
    #"""
