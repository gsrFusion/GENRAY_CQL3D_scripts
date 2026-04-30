import numpy as np
import os, sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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
import getGfileDict
import shotToEqdsk
import helperFunctions as helper

print(f'past imports')

machine = 'NTPT'

#times = ['04525', '02800','04300','04500','04000','03500','03300','02300','03500','03200']
#shots = ['147634','172550','186514','186651','180758','195081','176878','167502','175273','176034']

if machine == 'NTPT':
    time = '.V3APT'
    shot = 'ARC'
    whatCode = 'both'

    print(f'starting making scans')
    intermediateDir = 'LFSVersion/'
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{intermediateDir}{machine}_{shot}{time}'
    if '193765' in time:
        grillHeights = np.round(np.linspace(-.5,.5,11),3)
        NPara_targets = np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])

    elif '147634' in time:
        grillHeights = np.round(np.linspace(-.75,.75,13),3)
        NPara_targets = -1*np.array([2.5,2.6,2.7,2.8,2.9,3.0,3.1])

    elif 'V3A' in time:
        grillHeights = np.round(np.linspace(-1.75,1.75,15),3)[:1]
        NPara_targets = -1*np.array([1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2])[:1]#-1*np.array([1.4,1.6,1.8,2.0,2.2])


    if shot == 'DIIID':
        power = 1
    elif shot == 'ARC':
        power = 10

    pwrscale = power
    thgrills = np.zeros(len(grillHeights))

    eqdskName = shotToEqdsk.getEqdskName(f'{shot}{time}', machine)
    gfileDict = getGfileDict.getGfileDict(f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/')

    LCFS_R = gfileDict['rbbbs']
    LCFS_Z = gfileDict['zbbbs']

    argmaxZ = np.argmax(LCFS_Z)
    HFS_mask = np.where(LCFS_R <= LCFS_R[argmaxZ])[0]
    LFS_mask = np.where(LCFS_R >= LCFS_R[argmaxZ])[0]
    HFS_R = LCFS_R[HFS_mask]
    HFS_Z = LCFS_Z[HFS_mask]

    LFS_R = LCFS_R[LFS_mask]
    LFS_Z = LCFS_Z[LFS_mask]
    Rmaxis = gfileDict['rmaxis']

    for j, grillHeight in enumerate(grillHeights):
        nearestIndex = helper.findNearestIndex(grillHeight, LFS_Z)
        nearZ = LFS_Z[nearestIndex]
        nearR = LFS_R[nearestIndex]

        deltaR = nearR - Rmaxis

        rad = np.arctan2(nearZ, deltaR)
        deg = rad*180/np.pi % 360
        thgrills[j] = np.round(deg,3)

        print(f'thgrill: {thgrills[j]}')

        fig,ax = plt.subplots()
        ax.scatter(LFS_R, LFS_Z)
        ax.scatter([Rmaxis], [0])
        ax.axhline(-1.75)
        radius = np.sqrt(nearZ**2 + nearR**2)
        ax.scatter([np.cos(rad)+Rmaxis],[np.sin(rad)])
        ax.set_aspect('equal')
        plt.show()


        car = los

        
    for i in range(len(NPara_targets)):
        NPara_target = NPara_targets[i]
        prefix = 'n'
        factor = 1
        if NPara_target > 0:
            prefix = 'p'
            factor = -1

        N_para_peaks = np.array([NPara_target])#-2.7
        N_para_edges = np.array([[N_para_peaks[0]-.2, N_para_peaks[0] + .2]])#None
        powerInLobes = np.array([1e6])

        for j in range(len(thgrills)):
            thgrill = thgrills[j]

            targetDir = f'{stem}_{prefix}{np.abs(NPara_target)}Npara_{grillHeights[j]}grillHeight_{power}MW'
            print(f'targetDir: {targetDir}')
            doPlot = False
            if i + j > 0:
                doPlot = False

            helper = setupInputFiles.InputFileHelper(targetDir,  
            waveType = 'LH',
            makeDir = True, overwrite = True, doPlot = doPlot,
            numCQLToFokkerPlanck = 50, ndens = 101, njene= 101, 
            includeE = False, isScoping = True, eqsym = 'average',
            thgrill=thgrill, powerInLobes = powerInLobes,  N_para_edges = N_para_edges, 
            pwrScale = pwrscale, N_para_peaks = N_para_peaks,
            )
            helper.copySetupAndClean() 
            if whatCode == 'CQL3D':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runCQL.sh {targetDir}')
            elif whatCode == 'both':
                os.system(f'cd {targetDir} && bash -s < /home/grantr/codes/GENRAY_CQL3D_scripts/runGENThenCQL.sh {targetDir}')
            elif whatCode == 'GENRAY':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGEN.sh {targetDir}')


if machine == 'FENIX':

    print(f'starting making scans')
    for i in range(len(NPara_fors)):
        NPara_for = NPara_fors[i]
        prefix = 'n'
        if NPara_for > 0:
            prefix = 'p'

        """  
        targetNpara = NPara_for
        inputTarget = targetNpara
        
        if targetNpara > 0:
            inputTarget *= - 1
        N_para_peaks, N_para_edges, directivities = generateNparaSpectrum.generateSpectrum(inputTarget, analytic = False, powerRatio = [1,1,1,0,0,0,0,0], doPlot = False)
        if targetNpara > 0:
            N_para_peaks = -N_para_peaks
            N_para_edges = -1*np.flip(N_para_edges,axis = 1)

        for k in range(len(directivities)):
            if directivities[k]/directivities[0] < .05:
                directivities = directivities[:k]
                N_para_peaks = N_para_peaks[:k]
                N_para_edges = N_para_edges[:k]
                break

        powerInLobes = directivities*1e6#for 1 MW of forward power
        """
        N_para_peaks = np.array([NPara_for])
        N_para_edges = np.array([[N_para_peaks[0]-.2, N_para_peaks[0] + .2]])

        powerInLobes = np.array([1e6])#for 1 MW of forward power
        #for j in range(len(thgrills)):
            #thgrill = thgrills[j]
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/npara_thgrill_scan/{machine}_{shot}{time}'
        for j in range(len(thgrills)):
            thgrill = thgrills[j]
            doPlot = False
            if i == 0 and j == 0:
                doPlot = False
            targetDir = f'{stem}_{prefix}{np.abs(NPara_for):.2f}Npara_{thgrill}thgrill_1MW'#_3modules_{nScales[j]}nScale_{TScales[j]}Tscale'
            print(f'targetDir: {targetDir}')
            pwrScale = 1
            innerGap = 0

            helper = setupInputFiles.InputFileHelper(targetDir,  
                    waveType = 'LH',
                    makeDir = True, overwrite = True, doPlot = doPlot,
                    nScale = 1, TScale = 1, ZeffScale = 1,  
                    numCQLToFokkerPlanck = 50, ndens = 101, njene= 101,
                    includeE = False, isScoping = True, eqsym = 'average',
                    thgrill=thgrill, powerInLobes = powerInLobes,  N_para_edges = N_para_edges, pwrScale = 1, N_para_peaks = N_para_peaks,
                    )
            helper.copySetupAndClean() 


            if whatCode == 'CQL3D':
                os.system(f'ssh grantr@eofe7.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runCQL.sh {targetDir}')
            elif whatCode == 'both':
                os.system(f'ssh grantr@eofe7.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGENThenCQL.sh {targetDir}')
            elif whatCode == 'GENRAY':
                os.system(f'ssh grantr@eofe7.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGEN.sh {targetDir}')
                #"""
if machine == 'DIIID':
    time = '.05500'
    shot = '193765'
    whatCode = 'both'

    print(f'starting making scans')
    intermediateDir = 'Npara_height_scan/'
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{intermediateDir}{machine}_{shot}{time}'
    #power = 100#kW
    #100#1000#100#200 #kw
    #pwrscale = power/1000#/1000
    pwrscale = 1
    NPara_targets = np.array([2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1])
    
    grillHeights = np.round(np.linspace(-.5,.5,11),3)
    thgrills = np.zeros(len(grillHeights))

    eqdskName = shotToEqdsk.getEqdskName(f'{shot}{time}', machine)
    gfileDict = getGfileDict.getGfileDict(f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/')

    LCFS_R = gfileDict['rbbbs']
    LCFS_Z = gfileDict['zbbbs']

    argmaxZ = np.argmax(LCFS_Z)
    argminZ = np.argmin(LCFS_Z)
    HFS_mask = np.where(((LCFS_Z >= 0) & (LCFS_R <=LCFS_R[argmaxZ]))|((LCFS_Z <= 0) & (LCFS_R <=LCFS_R[argminZ])))[0]
    HFS_R = LCFS_R[HFS_mask]
    HFS_Z = LCFS_Z[HFS_mask]

    HFS_Z_interp = np.linspace(np.min(HFS_Z), np.max(HFS_Z),500)
    HFS_R_interp = interp1d(HFS_Z, HFS_R)(HFS_Z_interp)

    Rmaxis = gfileDict['rmaxis']
    Zmaxis = gfileDict['zmaxis']

    for j, grillHeight in enumerate(grillHeights):
        nearestIndex = helper.findNearestIndex(grillHeight, HFS_Z_interp)
        nearZ = HFS_Z_interp[nearestIndex]
        nearR = HFS_R_interp[nearestIndex]

        deltaR = nearR - Rmaxis
        deltaZ = nearZ - Zmaxis

        rad = np.arctan2(deltaZ, deltaR)
        deg = rad*180/np.pi % 360
        thgrills[j] = np.round(deg,3)

        """
        print(f'grill height: {grillHeight}, angle: {deg}')
        fig,ax = plt.subplots()
        ax.scatter(HFS_R_interp, HFS_Z_interp)
        ax.scatter(LCFS_R, LCFS_Z)
        ax.scatter([nearR],[nearZ])
        ax.scatter([Rmaxis],[Zmaxis],color = 'k')

        ax.scatter([np.cos(rad)+Rmaxis],[np.sin(rad)])

        ax.set_aspect('equal')
        plt.show()
        """

    for i in np.arange(len(NPara_targets)):
        NPara_target = NPara_targets[i]

        prefix= 'n'
        factor = 1
        if NPara_target > 0:
            prefix = 'p'
            factor = -1
        #"""
        N_para_peaks = np.array([NPara_target])
        N_para_edges = np.array([[NPara_target-.2, NPara_target+.2]])
        directivities = np.array([1])
        #"""
        """
        N_para_peaks, N_para_edges, directivities,_ = generateNparaSpectrum.generateSpectrum(target_npara = factor*NPara_target,
                                                                                                #modulePhaseShift = (factor*NPara_target+1.63229)/-.22314,
                                                                                                analytic = True, 
                                                                                                doPlot = False, 
                                                                                                num_module=8, 
                                                                                                #w_spacer=0.005, 
                                                                                                #delta = 1e-3/2,
                                                                                                )
        if NPara_target > 0:
            N_para_peaks = -N_para_peaks
            N_para_edges = -1*np.flip(N_para_edges,axis = 1)

        for k in range(len(directivities)):
            if directivities[k]/directivities[0] < .10:
                directivities = directivities[:k]
                N_para_peaks = N_para_peaks[:k]
                N_para_edges = N_para_edges[:k]
                break
        """

        powerInLobes = directivities*1e6#for 1 MW of forward power
        #targetDir = f'{stem}_{prefix}{np.abs(NPara_target)}Npara_1MW'
        #"""
        for j in np.arange(len(grillHeights)):
            grillHeight = grillHeights[j]
            thgrill = thgrills[j]
        
            targetDir = f'{stem}_{prefix}{np.abs(NPara_target)}Npara_{grillHeight}grillHeight_1MW'
          
            doPlot = True
            if i > 0 or j > 0:
                doPlot = False
            """
            for j in range(len(islandWidths)):
                width = islandWidths[j]
                for k in range(len(deltaTs)):
                    deltaT = deltaTs[k]

                    targetDir = f'{stem}_{prefix}{np.abs(NPara_target)}Npara_{power}MW_{width}Width_{deltaT}deltaT'

                    if i + j + k > 0:
                        doPlot = False
                    print(f'targetDir: {targetDir}')
                    
                    islandParamDict={'width':width, 'deltaT':deltaT, 'islandq' : 2}
            """
            helper = setupInputFiles.InputFileHelper(targetDir,  
                waveType = 'LH',
                makeDir = True, overwrite = True, doPlot = doPlot,
                nScale = 1, TScale = 1, ZeffScale = 1,  
                numCQLToFokkerPlanck = 50, ndens = 101, njene= 101, 
                includeE = False, isScoping = True, eqsym = 'average',
                thgrill=thgrill, powerInLobes = powerInLobes,  N_para_edges = N_para_edges, 
                pwrScale = pwrscale, N_para_peaks = N_para_peaks,
                #makeIsland = True,islandParamDict=islandParamDict,
                #rya = rya,
                )
            helper.copySetupAndClean() 
            if whatCode == 'CQL3D':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runCQL.sh {targetDir}')
            elif whatCode == 'both':
                os.system(f'ssh grantr@vlogin001 "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGENThenCQL.sh {targetDir}')
            elif whatCode == 'GENRAY':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGEN.sh {targetDir}')

if machine == 'KSTAR':
    time = '.009250'
    shot = '39608'
    whatCode = 'both'

    print(f'starting making scans')
    intermediateDir = 'Npara_thgrill_scan_5GHz/'
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{intermediateDir}{machine}_{shot}{time}'
    pwrscale = 1
    NPara_targets = -1*np.array([2.1,2.3,2.5,2.7,2.9])
    thgrills = np.array([140,160,180,200,220])

    for i in np.arange(len(NPara_targets)):
        NPara_target = NPara_targets[i]

        prefix= 'n'
        factor = 1
        if NPara_target > 0:
            prefix = 'p'
            factor = -1

        N_para_peaks = np.array([NPara_target])
        N_para_edges = np.array([[NPara_target-.2, NPara_target+.2]])
        directivities = np.array([1])

        powerInLobes = directivities*1e6#for 1 MW of forward power
        for j in np.arange(len(thgrills)):
            thgrill = thgrills[j]
            
            targetDir = f'{stem}_{prefix}{np.abs(NPara_target)}Npara_{thgrill}thgrill_1MW_5GHz'
            doPlot = True
            if i > 0 or j > 0:
                doPlot = False
           
            helper = setupInputFiles.InputFileHelper(targetDir,  
                waveType = 'LH',
                makeDir = True, overwrite = True, doPlot = doPlot,
                nScale = 1, TScale = 1, ZeffScale = 1,  
                numCQLToFokkerPlanck = 50, ndens = 101, njene= 101, 
                includeE = False, isScoping = True, eqsym = 'average',
                thgrill=thgrill, powerInLobes = powerInLobes,  N_para_edges = N_para_edges, 
                pwrScale = pwrscale, N_para_peaks = N_para_peaks,
                frqncy = 5e9,
                )
            helper.copySetupAndClean() 
            if whatCode == 'CQL3D':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runCQL.sh {targetDir}')
            elif whatCode == 'both':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGENThenCQL.sh {targetDir}')
            elif whatCode == 'GENRAY':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGEN.sh {targetDir}')


if machine == 'MANTA':
    print(f'starting making scans')
    NPara_fors = np.arange(-2.25,-1.25 + .1,.1)
    thgrills = np.arange(110, 250 + 10, 10)
    for i in range(len(NPara_fors)):
        NPara_for = NPara_fors[i]
        prefix = 'n'
        if NPara_for > 0:
            prefix = 'p'

        N_para_peaks = np.array([NPara_for])
        N_para_edges = np.array([[N_para_peaks[0]-.2, N_para_peaks[0] + .2]])

        powerInLobes = np.array([1e6])#for 1 MW of forward power
        #for j in range(len(thgrills)):
            #thgrill = thgrills[j]
        stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/npara_thgrill_scan/{machine}_{shot}{time}'
        for j in range(len(thgrills)):
            thgrill = thgrills[j]
            doPlot = False
            if i +j > 0:
                doPlot = False
            targetDir = f'{stem}_{prefix}{np.abs(NPara_for):.2f}Npara_{thgrill}thgrill_1MW'#_3modules_{nScales[j]}nScale_{TScales[j]}Tscale'
            print(f'targetDir: {targetDir}')
            pwrScale = 1
            innerGap = 0

            helper = setupInputFiles.InputFileHelper(targetDir,  
                    waveType = 'LH',
                    makeDir = True, overwrite = True, doPlot = False,
                    nScale = 1, TScale = 1, ZeffScale = 1,  
                    numCQLToFokkerPlanck = 50, ndens = 101, njene= 101,
                    includeE = False, isScoping = True, eqsym = 'average',
                    thgrill=thgrill, powerInLobes = powerInLobes,  N_para_edges = N_para_edges, pwrScale = 1, N_para_peaks = N_para_peaks,
                    )
            helper.copySetupAndClean() 


            if whatCode == 'CQL3D':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runCQL.sh {targetDir}')
            elif whatCode == 'both':
                os.system(f'ssh vlogin001 "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGENThenCQL.sh {targetDir}')
            elif whatCode == 'GENRAY':
                os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGEN.sh {targetDir}')
                    #"""