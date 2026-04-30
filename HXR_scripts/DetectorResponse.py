import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from numpy import exp, sqrt
from scipy.special import erfc

import DetectorEfficiency

###
# First apply detector efficiencies
# Then apply response function matrices
###

def applyDetectorEfficiency(en_, efluxNormed):
    dE = en_[1]-en_[0]
    energies = en_-dE/2
    efficiencies = DetectorEfficiency.getEfficiencies(energies)
    assert np.max(efficiencies < 1)
    efluxInteract = efluxNormed * efficiencies
    
    return efluxInteract

#chordToChannelDict - dictionary going from chord to which channel is installed
# the chords relevant to interactMatrix
# interactMatrix - result of applying attenuation plates and detector efficiencies to eflux
#       these are the incident photon energies that actually interact with the detectors
# en_ - the energies at which the cql3d count rates were calculated at
def applyPHS(chords, interactMatrix, en_):
    import getChordToChannelDict
    chordToChannelDict = getChordToChannelDict.getChordToChannelDict()
    fig,ax = plt.subplots()
    newCountMatrix = np.zeros(interactMatrix.shape)
    for chordIndex in range(len(chords)):
        chord = chords[chordIndex]
        channel = chordToChannelDict[chord]
        responseMatrix = np.load(f'/home/grantr/codes/HXRCameraResponseFunctions/responseFunctionMatrices/PHSMatrix_ch{channel}.npy')
        energyForRespMatrix = np.load(f'/home/grantr/codes/HXRCameraResponseFunctions/responseFunctionMatrices/energyBasis_ch{channel}.npy')

        assert np.sum(en_ - energyForRespMatrix) == 0

        #normalize so the sum over each response function is 1
        responseMatrix = responseMatrix/responseMatrix.sum(axis=1)[:,None]

        thermalCountsOfChord = interactMatrix[0,chordIndex,:]
        nonThermalCountsOfChord = interactMatrix[1,chordIndex,:]

        newTCountsOfChord = np.zeros(len(thermalCountsOfChord))
        newNTCountsOfChord = np.zeros(len(nonThermalCountsOfChord))

        for energyIndex in range(len(en_)):
            newTCountsOfChord[energyIndex] = np.sum(thermalCountsOfChord[:]*responseMatrix[:,energyIndex])
            newNTCountsOfChord[energyIndex] = np.sum(nonThermalCountsOfChord[:]*responseMatrix[:,energyIndex])
            
            ax.plot(en_, (thermalCountsOfChord[energyIndex]*responseMatrix[energyIndex,:]))
        newCountMatrix[0, chordIndex, :] = newTCountsOfChord
        newCountMatrix[1, chordIndex, :] = newNTCountsOfChord
        print(f'np.sum(newCountMatrix): {np.sum(newCountMatrix)}, np.sum(interactMat: {np.sum(interactMatrix)})')

    assert np.round(np.sum(interactMatrix),10) == np.round(np.sum(newCountMatrix),10)

    print(f'diff between w/wo resp, in det resp: {np.sum(np.abs(interactMatrix - newCountMatrix))}')

    return newCountMatrix

def applyDetectorResponse(chords, efluxNormed, en_):
    interactMatrix = applyDetectorEfficiency(en_, efluxNormed)
    measuredMatrix = applyPHS(chords, interactMatrix, en_)
    return measuredMatrix