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

time = '.01980'
shot = '206629'
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

    if shot == '203920' and time == '.03100':
        modulePhases_DIIID = np.array([0, 0, 44.6, 0.0, 327.1, 0, 72.5, 0])
        modulePowers = np.array([0.0, 0.0, 7.17, 23.75, 33.8, 0.0, 38.9, 0.0])
        isNormalBt = True

    if shot == '205015' and time == '.02500':
        modulePhases_DIIID = np.array([0, 0, 0, 0.0, 350.5, 0, 173.0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 32.0, 41.0, 0.0, 24.0, 0.0])
        isNormalBt = False

    if shot == '201877-05Ip' and time == '.01500':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.60000000000002, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = False

    if shot == '201877-07Ip' and time == '.01500':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.60000000000002, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = False

    if shot == '201877' and time == '.01500':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.60000000000002, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = False

    if shot == '200858' and time == '.02000':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.60000000000002, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = False

    if shot == '200858' and time == '.01200':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.60000000000002, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = False

    if shot == '203617' and time == '.04150':
        modulePhases_DIIID = np.array([0, 0, 0, 0, 7.0, 181.0, 0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 79.0, 62.0, 0.0, 0.0])
        calcThgrill = True
        isNormalBt = False

    if shot == '203617' and time == '.04120':
        modulePhases_DIIID = np.array([0, 0, 0, 0, 6.42, 179.32, 0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 45.0, 41.0, 0.0, 0.0])
        calcThgrill = True
        isNormalBt = False

    if shot == '203617' and time == '.04140':
        calcThgrill = True
        isNormalBt = False
        modulePhases_DIIID = np.array([0, 0, 0, 0, 4.82, 177.94, 0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 69.37, 55.68, 0.0, 0.0])

    if shot == '203619' and time == '.04160':
        modulePhases_DIIID = np.array([0, 0, 0, 0, 6.7,180.1, 0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 88.95,70.2, 0.0, 0.0])
        isNormalBt = False

    if shot == '203619' and '.04130' in time:
        modulePhases_DIIID = np.array([0, 0, 0, 0, 4.5, 177.13, 0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 59.18,48.54, 0.0, 0.0])
        calcThgrill = True
        isNormalBt = False

    if shot == '203574' and time == '.01600':
        modulePhases_DIIID = np.array([0, 0, 0, 0, 312.0, 181.0, 146.0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 11.0, 65.0, 19.0, 0.0])
        isNormalBt = True

    if shot == '203575' and time == '.02400':
        modulePhases_DIIID = np.array([0, 0, 0, 0, 311.0, 264.0, 146.0, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 0.0, 11.0, 43.0, 19.5, 0.0])
        isNormalBt = True

    if shot == '204655' and time == '.02400':
        modulePhases_DIIID = np.array([0, 0, 0, 0.0, 348.0, 0, 179.5, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 20.0, 50.0, 0.0, 13.6, 0.0])
        isNormalBt = True

    if shot == '204654' and time == '.01780':
        modulePhases_DIIID = np.array([0, 0, 0, 0.0, 347.9, 0, 178.9, 0])
        modulePowers = np.array([0.0, 0.0, 0.0, 14.65, 44.0, 0.0, 14.0, 0.0])
        isNormalBt = True

    if shot == '203915' and time == '.02200':
        modulePhases_DIIID = np.array([0, 0, 47.0, 0.0, 329.0, 0, 76.0, 0])
        modulePowers = np.array([0.0, 0.0, 6.2, 25.2, 36.7, 0.0, 38.3, 0.0])
        isNormalBt = True

    if shot == '203913' and time == '.02200':
        modulePhases_DIIID = np.array([0, 0, 49.6, 0.0, 331.0, 0, 78.0, 0])
        modulePowers = np.array([0.0, 0.0, 4.93, 23.0, 36.7, 0.0, 39.9, 0.0])
        isNormalBt = True

    if shot == '204535' and time == '.01800':
        modulePhases_DIIID = np.array([0, 0, 15.0, 0.0, 209.0, 0, 13.0, 0])
        modulePowers = np.array([0.0, 0.0, 28.0, 33.5, 24.0, 0.0, 41.0, 0.0])
        isNormalBt = True

    if shot == '203912' and time == '.02700':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.6, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = True

    if shot == '203912' and time == '.02780':
        modulePhases_DIIID = np.array([0, 0, 50, 0, 331.2, 0, 75.74,0])
        modulePowers = np.array([0, 0, 5.4, 23.4, 37,0,38.9,0])
        calcThgrill = True
        isNormalBt = True

    if shot == '203912' and time == '.02840':
        modulePhases_DIIID = np.array([0, 0, 47.5, 0.0, 329.4, 0, 76.2, 0])
        modulePowers = np.array([0.0, 0.0, 5.4, 23.3, 37.6, 0.0, 39.5, 0.0])
        isNormalBt = True

    if shot == '203692' and time == '.02300':
        modulePhases_DIIID = np.array([0, 0, 44.0, 0.0, 328.5, 0, 72.60000000000002, 0])
        modulePowers = np.array([0.0, 0.0, 5.7, 23.2, 37.0, 0.0, 39.0, 0.0])
        isNormalBt = False

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

    if shot == '203924' and time == '.03250':
        modulePhases_DIIID = np.array([0, 0, 42.0, 0.0, 328.0, 0, 71.2, 0])
        modulePowers = np.array([0.0, 0.0, 4.3, 27.0, 33.2, 0.0, 40.0, 0.0])
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
    #modulePhases_Andrew = modulePhases_DIIID-modulePhases_DIIID[0]
    modulePhases_Andrew = np.zeros(8)
    """
    for i in range(len(modulePhases_Andrew)):
        if i == 0:
            modulePhases_Andrew[i] = 0
        else:
            modulePhases_Andrew[i] = modulePhases_DIIID[i] - modulePhases_DIIID[i-1] + modulePhases_Andrew[i-1]
    print(f'new: {modulePhases_Andrew}, orig; {modulePhases_DIIID-modulePhases_DIIID[0]}')
    """
    modulePhases_Andrew = modulePhases_DIIID
    #modulePhases_Andrew = modulePhases_Andrew - modulePhases_Andrew[2]
    modulePhases_Andrew %= 360

    modulePhases_Andrew_rad = np.pi*modulePhases_Andrew/180

    peakNparas, peakEdges, directivities, P_total = generateNparaSpectrum.generateSpectrum(target_npara=None, #if not supplying the phase shift, what is your target N||
                    modulePhaseShift = modulePhases_Andrew_rad, #RADIANS, if not supplying the target N||, what is the phase shift
                    analytic = False, #if you want to calculate things analytically. Only valid for equal powers in all WGs
                    modulePowerRatio = modulePowers_normed, #power ratios between in the modules
                    #wgPowerRatio = None, # power ratios in the WGs within each module
                    numLobes = 4,
                    doPlot = 'spectrum',) #whether or not to plot
    
    peakNparas = np.array([peakNparas[0], peakNparas[2], peakNparas[3]])
    peakEdges = np.array([peakEdges[0], peakEdges[2], peakEdges[3]])
    directivities = np.array([directivities[0], directivities[2], directivities[3]])

    print(peakNparas)
    print(peakEdges)
    print(directivities)

    car = los
    if isNormalBt:
        peakNparas = -peakNparas
        peakEdges = -1*np.flip(peakEdges,axis = 1)
    totalPower_kW = np.sum(modulePowers)
    pwrFactor = 1
    pwrscale = 0.72*totalPower_kW/1000*pwrFactor
    print(f'starting input file helper')
    intermediateDir = ''#'60nnkpar_1e-8prmt4_1e-7prmt4ECE/'
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{intermediateDir}{machine}_{shot}{time}'

    targetDir = f'{stem}_expSpectrum'

    #jesus christ this shouldn't be legal
    """
    k = 1
    while k < len(directivities):
        print(k)
        if directivities[k]/directivities[0] < .10:
            print(directivities[k]/directivities[0] )
            break
        k += 1
    """    
    #print(f'k: {k}, len direct: {len(directivities)}')
    #directivities = directivities[:k]
    #N_para_peaks = peakNparas[:k]
    #N_para_edges = peakEdges[:k]
    N_para_peaks = peakNparas
    N_para_edges = peakEdges
    powerInLobes = directivities*1e6#for 1 MW of forward power
        
    print(f'N_para_peaks: {N_para_peaks}')

    thgrill = 189
    if calcThgrill:
        print(f'yes')
        thgrill = np.inf
    #"""
    helper = setupInputFiles.InputFileHelper(targetDir,  
    waveType = 'LH',
    makeDir = True, overwrite = True, doPlot = True,
    nScale = 1, TScale = 1, ZeffScale = 1,  
    numCQLToFokkerPlanck = 50, ndens = 101, njene= 101, 
    includeE = False, isScoping = False, eqsym = 'average',
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
