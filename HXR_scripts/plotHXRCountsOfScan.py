import numpy as np
import CountMatrix
import os, sys
import netCDF4
import matplotlib.pyplot as plt

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import helperFunctions as helper
import getInputFileDictionary

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)

def getAvgSigCounts(targetDir):
    chords = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])
    _, countMatrix = CountMatrix.getCountMatrix(chords, targetDir = targetDir, attenuate = True, includeResponseFunc = False, E_pMin=50, E_pMax=250)
    countsPerChord = np.sum(countMatrix, axis = 2)
    avgNonthermalCounts = np.sum(countsPerChord[1])/len(countsPerChord[1])

    return avgNonthermalCounts

def getCurrent(targetDir):
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
    darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
    totalCD = np.sum(curr*darea)
    
    return totalCD

def createScan():
    rootDir = f'/home/grantr/scratch/genray_batch/DIIID_shots/'
    folderList = fast_scandir(rootDir)
    
    fig, ax = plt.subplots()

    mostCountsDir = ''
    mostCounts = -1

    for folder in folderList:
        try:
            if '300kW' in folder and not ('190316' in folder) and '2.7' in folder:
                cqlrf_nc = netCDF4.Dataset(f'{folder}/cql3d_krf001.nc','r')
                cql_nc = netCDF4.Dataset(f'{folder}/cql3d.nc','r')
                genray_in = getInputFileDictionary.getInputFileDictionary('genray', targetDir = folder)
                SPA = helper.getSPA(cqlrf_nc, genray_in)
                if SPA > 0.95:

                    current  = helper.getCurrent(cql_nc)
                    avgSigCounts = getAvgSigCounts(folder)

                    print(f'{folder}\n counts: {avgSigCounts}, SPA: {SPA}, current: {current}')

                    if avgSigCounts > mostCounts:
                        mostCounts = avgSigCounts
                        mostCountsDir = folder

                    ax.scatter([current], [avgSigCounts], color = 'tab:blue', s = 10)
        except:
            pass

    print(f'most counts occured in sim: {mostCountsDir}')

    ax.set_ylabel('Avg nonthermal HXR count rate')
    ax.set_xlabel(r'$I_{LH}$ (MA)')
    ax.set_title(f'All 300 kW sims with SPA > 0.95')
    fig.tight_layout()
    plt.show()

def fast_scandir(dirname):
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders  


createScan()