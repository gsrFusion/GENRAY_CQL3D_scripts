import numpy as np
import matplotlib.pyplot as plt
import CountMatrix
from scipy.interpolate import interp1d
import netCDF4
import os,sys

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getInputFileDictionary
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d')
genray_in = getInputFileDictionary.getInputFileDictionary('genray')

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

#number of detectors
nv = cqlinput['setup']['nv']
#min photon energy we are looking for
emin = cqlinput['setup']['enmin']
#max photon energy we are looking for
emax = cqlinput['setup']['enmax']

#energy bins used by CQL3D for the XR detector
en_ = cql_nc.variables["en_"][:]

_imageBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])
_28 = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65])
_LFS = np.array([2,3,4,5,11,12,13,14,15])
_core = np.array([21,22,23,24,25, 31,32,33,34,35, 41,42,43,44,45])
_HFS = np.array([52,53,54,55, 62,63,64,65])
_all = np.arange(1,nv+1,1).astype(int)

campaign2024 = np.array([2,12,22,32,42,52,62, 3,23,33,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])
_28 = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65])
chords = _28

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

#plots the number of counts seen by each chord, summed over energy   
def addCountsPerChord(ax, efluxOfInt):
    
    #sum over energy bins to get counts/chord
    countsPerChord = np.sum(efluxOfInt, axis = 2)
    assert countsPerChord.shape == (2, len(chords)) or countsPerChord.shape == (3, len(chords))
    
    totalCounts = countsPerChord[0] + countsPerChord[1]     
    ax.scatter(chords, countsPerChord[1], label = "Non Thermal\nBremsstrahlung", color = 'b', s = 11)
    ax.scatter(chords, countsPerChord[0], label = "Thermal\nBremsstrahlung", color = 'g', s = 11)
    #ax.scatter(chords, totalCounts, label = "total counts", color = 'r', s = 11)

    avgNonthermalCounts = np.sum(countsPerChord[1])/len(countsPerChord[1])

    print(f"min xrays per s: {np.min(totalCounts): 4e}")
    print(f"max xrays per s: {np.max(totalCounts): 4e}")
    print(f"avg xrays per s: {np.sum(totalCounts)/len(totalCounts): 4e}")
    print(f"avg nonthermal xrays per s: {avgNonthermalCounts: 4e}")


#plots the number of counts per energy bin, summed over all chords
def addCountsPerEnergy(ax, en_OfInterest, efluxOfInt, linestyle = 'solid'):
    #sum over chords to get energy per bin of en_OfInterest
    countsPerEnergy = np.sum(efluxOfInt, axis = 1)

    ax.plot(en_OfInterest, countsPerEnergy[1], label = "Non Thermal\nBremsstrahlung", color = 'b', lw=2, linestyle = linestyle)
    ax.plot(en_OfInterest, countsPerEnergy[0], label = "Thermal\nBremsstrahlung", color = 'g', lw=2, linestyle = linestyle)

#plots the time resolution of each chord
def addTimeResolutions(ax, efluxOfInt):
    #sum over chords to get energy per bin of en_OfInterest
    countsPerChord = np.sum(efluxOfInt, axis = 2)

    factor = 1.687

    nonThermalCounts = countsPerChord[1]*factor

    timeResolutions = 1/(.1**2 * nonThermalCounts)
    print(f'avgTimeRes: {np.average(timeResolutions)}')
    print(f'median time resolution: {np.median(timeResolutions)}')
    ax.scatter(chords, timeResolutions, label = "Non Thermal\nBremsstrahlung", color = 'b', s = 11)

#plot counts/s vs chord and counts/s vs energy bins
def plotNormalizedCounts():
    fig, axes = plt.subplots(nrows = 2)
#    assert len(dischargeNumber) == 6


    #######read in input files to get simulation parameters#######
    power1 = genray_in["grill"]["powers(1)"]
    power2 = genray_in["grill"]["powers(2)"]
    scaledPower1 = power1 * cqlinput["rfsetup"]["pwrscale(1)"]
    scaledPower2 = power2 * cqlinput["rfsetup"]["pwrscale(1)"]

    absorbedPowerMW = (scaledPower1 + scaledPower2)/1e6
    try:
        absorbedPowerMW += genray_in["grill"]["powers(3)"]*cqlinput["rfsetup"]["pwrscale(1)"]/1e6
    except:
        pass

    n_para_f = (genray_in["grill"]["anmax(1)"] + genray_in["grill"]["anmin(1)"])/2
    n_para_r = (genray_in["grill"]["anmax(2)"] + genray_in["grill"]["anmin(2)"])/2
    ##############################################################


    #######get the data to plot######
    #get the energy points of interest and the corresponding count matrix
    E_pMax = 250; E_pMin = 30
    en_OfInterest, countMatrix = CountMatrix.getCountMatrix(chords, attenuate = True, includeResponseFunc = False, E_pMin=E_pMin, E_pMax=E_pMax)
    _, countMatrix_noResp = CountMatrix.getCountMatrix(chords, attenuate = True, includeResponseFunc = False, E_pMin=E_pMin, E_pMax=E_pMax)

    print(f'{en_OfInterest.shape, countMatrix.shape}')

    addCountsPerChord(axes[0], countMatrix)
    #addTimeResolutions(axes[1], countMatrix)
    addCountsPerEnergy(axes[1], en_OfInterest, countMatrix)
    addCountsPerEnergy(axes[1], en_OfInterest, countMatrix_noResp, linestyle = 'dashed')

    print(f'average nonthermal counts per chord between {E_pMin} and {E_pMax} keV: {np.sum(countMatrix[1])/len(countMatrix[1])}')
    print(f'average thermal counts per chord between {E_pMin} and {E_pMax} keV: {np.sum(countMatrix[0])/len(countMatrix[0])}')
    ################################


    #######do plotting stuff#######
    #fig.tight_layout(rect=[0, 0.03, 1, 0.93])
    fig.suptitle(f"Discharge: {shotNum},  Coupled Power: {absorbedPowerMW : .2f} MW" + 
            f"\n$n_{{||,f}}$ = {n_para_f: 0.2f},  $n_{{||,r}}$ = {n_para_r: 0.2f},"+\
            r"  $E_{{p}} \in $" + f"[{E_pMin}, {E_pMax}] keV")
    #fig.set_size_inches(7,8)
    #"""
    axes[0].set_xlabel("Chord number")
    axes[0].set_ylabel("Counts/s")
    axes[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axes[0].set_xlim([1,nv])
    axes[0].legend(loc = 'best', fontsize = 10)

    axes[1].set_xlabel("Energy (keV)")
    axes[1].set_ylabel("Counts/s")
    #axes[1].set_ylim(bottom=1e0)
    axes[1].set_xlim([en_OfInterest[0], en_OfInterest[-1]])
    #axes[1].set_yscale('log')
    legend1 = axes[1].legend(loc = 'best', fontsize = 10, ncol = 2, columnspacing = .5)
    legend1.get_frame().set_alpha(None)
    legend1.get_frame().set_facecolor((1, 1, 1, 0.1))
    fig.tight_layout()
    plt.show()
    ##############################
    
plotNormalizedCounts()
