import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

cqlinput = proj.model1.genray_cql3d.cql3d.cqlinput
#number of detectors
nv = cqlinput.get_contents('setup', 'nv')[0]
#min photon energy we are looking for
emin = cqlinput.get_contents('setup', 'enmin')[0]
#max photon energy we are looking for
emax = cqlinput.get_contents('setup', 'enmax')[0]

nc = proj.model1.genray_cql3d.cql3d.cql3d_nc
#energy bins used by CQL3D for the XR detector
en_ = np.ma.getdata(nc.get_contents("variables","en_"))
#count rates of the chords. eflux[0] is thermal bremsstrahlung and eflux[1] is the nonthermal 
eflux = nc.get_contents("variables","eflux")

print(f"len(en_):{len(en_)}")

def getEtendues():
    if nv == 66:
        return np.array([0.00012569, 0.00012847, 0.00013057, 0.00013188, 0.00013233, 0.00013188, 0.00013057, 0.00012847, 0.00012569, 0.00012726, 0.00013057, 0.00013323, 0.0001351, 0.00013606, 0.00013606, 0.0001351, 0.00013323, 0.00013057, 0.00012726, 0.00013188, 0.0001351, 0.00013755, 0.00013908, 0.0001396, 0.00013908, 0.00013755, 0.0001351, 0.00013188, 0.00013233, 0.00013606, 0.00013908, 0.00014121, 0.00014231, 0.00014231, 0.00014121, 0.00013908, 0.00013606, 0.00013233, 0.00013606, 0.0001396, 0.00014231, 0.00014401, 0.00014459, 0.00014401, 0.00014231, 0.0001396, 0.00013606, 0.0001351, 0.00013908, 0.00014231, 0.00014459, 0.00014578, 0.00014578, 0.00014459, 0.00014231, 0.00013908, 0.0001351, 0.00013755, 0.00014121, 0.00014401, 0.00014578, 0.00014638, 0.00014578, 0.00014401, 0.00014121, 0.00013755])
        #return np.array([0.0273, 0.0278, 0.0282, 0.0284, 0.0284, 0.0284, 0.0282, 0.0278, 0.0273, 0.0276, 0.0282, 0.0286, 0.0289, 0.0291, 0.0291, 0.0289, 0.0286, 0.0282, 0.0276, 0.0284, 0.0289, 0.0293, 0.0296, 0.0296, 0.0296, 0.0293, 0.0289, 0.0284, 0.0284, 0.0291, 0.0296, 0.0299, 0.0301, 0.0301, 0.0299, 0.0296, 0.0291, 0.0284, 0.0291, 0.0296, 0.0301, 0.0303, 0.0304, 0.0303, 0.0301, 0.0296, 0.0291, 0.0289, 0.0296, 0.0301, 0.0304, 0.0306, 0.0306, 0.0304, 0.0301, 0.0296, 0.0289, 0.0293, 0.0299, 0.0303, 0.0306, 0.0307, 0.0306, 0.0303, 0.0299, 0.0293])
    else:
        return np.array([.0273]*nv)
_etendues = getEtendues()
_imageBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])
_LFS = np.array([2,3,4,5,11,12,13,14,15])
_core = np.array([21,22,23,24,25, 31,32,33,34,35, 41,42,43,44,45])
_HFS = np.array([52,53,54,55, 62,63,64,65])
_all = np.arange(1,nv+1,1)

chords = _imageBand

_etendues = _etendues[chords - 1]
eflux = np.ma.getdata(eflux[:, chords - 1, :])

plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)
plt.rc('axes', labelsize = 14)
plt.rc('figure', titlesize = 16)

genray_in = proj.model1.genray_cql3d.genray.genray_in
cqlinput = proj.model1.genray_cql3d.cql3d.cqlinput

#plot counts/s vs chord and counts/s vs energy bins
def determineEnergyBinning():
    fig, ax = plt.subplots(dpi=100)

    efluxOfInt = proj.model1.genray_cql3d.HXRScripts.getCountMatrix(chords, True, False)

    countsPerEnergy = np.sum(efluxOfInt, axis = 1)
    countsPerEnergy = countsPerEnergy[0] + countsPerEnergy[1]
    countsPerEnergy*=.2
    dE = en_[1]-en_[0]
    energies = en_-dE/2

    energyCutoffs = [250]
    index = len(en_)-1
    numInGroup = 0
    
    allowedError = .1

    while index >=0:
        if numInGroup ==0:
            numInGroup += countsPerEnergy[index]
            
        if numInGroup !=0:
            numInGroup += countsPerEnergy[index]
            if 1/np.sqrt(numInGroup) < allowedError:
                numInGroup = 0
                energyCutoffs.append(energies[index])

        index -=1

    energyCutoffs = np.array(energyCutoffs)

    ax.plot(energies, countsPerEnergy)
    for cutOff in energyCutoffs:
        ax.axvline(cutOff)

    binWidths = np.abs(energyCutoffs[1:]-energyCutoffs[:-1])
    print(f"min bin size: {min(binWidths)}, max: {np.max(binWidths)}")
    ax.set_yscale('log')
    fig.show()


def findNearestIndex(value, array):
    idx = (np.abs(array - value)).argmin()
    return idx

ans(determineEnergyBinning())