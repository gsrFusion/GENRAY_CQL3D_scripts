import numpy as np
from scipy.interpolate import interp1d
import netCDF4

import AttenMat
import DetectorResponse
import DetectorInformation

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getTargetInfo
import getInputFileDictionary


def getCountMatrix(chords, targetDir = None, attenuate = True, includeResponseFunc = True, E_pMin= 50, E_pMax = 250):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d', targetDir=targetDir)

    en_ = np.ma.getdata(cql_nc.variables["en_"][:])#energy bins used by CQL3D for the XR detector
    nv = cqlinput['setup']['nv']#number of sight lines

    dE = en_[1]-en_[0]
    etendues =DetectorInformation.getEtendues()
      
    minIndex = findNearestIndex(E_pMin, en_)
    maxIndex = findNearestIndex(E_pMax, en_)
    #(thermal/nonthermal, chord, energy)
    eflux = np.ma.getdata(cql_nc.variables["eflux"][:])
    eflux = eflux[:,:,minIndex:maxIndex+1]
    en_OfInterest = en_[minIndex:maxIndex+1]

    energyBinsMesh, void= np.meshgrid(en_OfInterest, np.zeros(eflux.shape[1]))
    
    etenduesMesh = np.stack([etendues for i in range(len(en_OfInterest))], axis = -1)
    efluxNormed = (eflux* 624150974000 * etenduesMesh/energyBinsMesh)*dE
    efluxNormed = efluxNormed[:,chords-1, :]
    if attenuate:
        attenMat = AttenMat.getAttenFunc(chords, en_)
        #print(f'attenmat: {attenMat.shape}, efluxNo: {efluxNormed.shape}')
        efluxNormed *= attenMat[:,minIndex:maxIndex+1]

    #efluxNormed = np.zeros(efluxNormed.shape)
    #efluxNormed[:,:,:] = 1

    if includeResponseFunc is False:
        return en_OfInterest, efluxNormed
    else:
        respedMatrix = DetectorResponse.applyDetectorResponse(chords, efluxNormed, en_)[:,:,minIndex:maxIndex+1]
        print(f'diff between w/wo resp, in count matrix: {np.sum(np.abs(respedMatrix - efluxNormed[:,:,minIndex:maxIndex+1]))}')
        return en_OfInterest, respedMatrix


def findNearestIndex(value, array):
    idx = (np.abs(array - value)).argmin()
    return idx