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
print(f'past imports')

whatCode = 'both'


machine = 'NTPT'

if machine == 'NTPT':
    print(f'starting making scans')
    intermediateDir = ''
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/'
    power = 10#MW
    #100#1000#100#200 #kw
    pwrscale = power#/1000
    NPara_target = -2.1
    thgrill = 180
    
    prefix = 'n'
    factor = 1
    if NPara_target > 0:
        prefix = 'p'
        factor = -1

    N_para_peaks = np.array([NPara_target])#-2.7
    N_para_edges = np.array([[N_para_peaks[0]-.2, N_para_peaks[0] + .2]])#None
    powerInLobes = np.array([1e6])

    shotNums = ['MANTA']*11
    shotTimes = ['.NT05', '.NT04','.NT03','.NT02','.NT01','.NT00','.PT01','.PT02','.PT03', '.PT04', '.PT05']

    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/'

    for i in range(len(shotNums)):
        shotNum = shotNums[i]
        shotTime = shotTimes[i]
        

        targetDir = f'{stem}{machine}_{shotNum}{shotTime}/{machine}_{shotNum}{shotTime}_{prefix}{np.abs(NPara_target)}Npara_{thgrill}thgrill_{power}MW'
    
        doPlot = True
        if i > 0:
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
            os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGENThenCQL.sh {targetDir}')
        elif whatCode == 'GENRAY':
            os.system(f'ssh grantr@orcd-vlogin001.mit.edu "cd {targetDir} && bash -s" < /home/grantr/codes/GENRAY_CQL3D_scripts/runGEN.sh {targetDir}')

