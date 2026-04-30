"""
Plots the ray traces and the RF power deposition density
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from omfit_classes import omfit_eqdsk
from scipy.interpolate import interp1d

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import getGfileDict
gfileDict = getGfileDict.getGfileDict()
import helperFunctions as helper
import getInputFileDictionary
import getTargetInfo
import shotToEqdsk


def getAverageCDFigureMerit(targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()

    splitted = targetDir.split('/')[-1].split('_')
    shotNum = splitted[1]
    machine = splitted[0]

    eqdskName = shotToEqdsk.getEqdskName(shotNum, machine = machine)
    eqdsk = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')

    genray_in = getInputFileDictionary.getInputFileDictionary('genray', targetDir=targetDir)
    cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d', targetDir=targetDir)

    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

    darea = cql_nc.variables["darea"][:]*1e-4#convert to m^2
    dvol = cql_nc.variables["dvol"][:]*1e-6#convert to m^3
    curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
    rya = np.ma.getdata(cql_nc.variables["rya"][:])

    currRelContribution = curr*darea
    totalCurrent = np.sum(currRelContribution)

    powrft = cql_nc.variables["powrft"][-1,:] #W/cm^3 = MW/m^3
    powrft[powrft<1e-5] = 0

    powrftRelContribution = powrft*dvol
    totalPower = np.sum(currRelContribution)

    avgR = eqdsk['fluxSurfaces']['avg']['R']
    eqdskRho_p = np.sqrt(eqdsk['fluxSurfaces']['levels'])
    avgR = interp1d(eqdskRho_p, avgR)(rya)

    avgRRelContribution = np.sum(curr*darea*avgR)/totalCurrent

    _, ne = helper.getCQLne(targetDir = targetDir, rho_pol = rya)
    neRelContribution = np.sum(curr*darea*ne)/totalCurrent
    

    figureOfMerit = (neRelContribution/1e20) * avgRRelContribution * totalCurrent / totalPower 
    print(figureOfMerit)
    """#figureOfMerit[powrft == 0] = 0
    print(avgRRelContribution)
    print(neRelContribution)
    fig,ax = plt.subplots()
    ax.plot(rya, powrft)
    ax.plot(rya, curr)
    ax.twinx().plot(rya, figureOfMerit)
    plt.show()"""
    return figureOfMerit

getAverageCDFigureMerit()