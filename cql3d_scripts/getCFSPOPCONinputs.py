###
# Plots the electron and ion densities and temperatures according to the cql3d input file
###

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__);dname = os.path.dirname(abspath);os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__));parentdir = os.path.dirname(currentdir);sys.path.append(parentdir)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

import getInputFileDictionary
inputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d')
import helperFunctions as helper
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

import netCDF4
import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
print(targetDir)

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

tiscal = 1
tescal = 1
enescal = 1

#if the temperature and density scale factors are in the input file, go get them
try:
    tiscal = inputFileDict['setup']['tiscal']
except:
    pass
try:
    tescal = inputFileDict['setup']['tescal']
except:
    pass
try:
    enescal = inputFileDict['setup']['enescal']
except:
    pass
T_e = inputFileDict['setup']['tein']*tescal
T_i = inputFileDict['setup']['tiin']*tiscal
n_e = inputFileDict['setup']['enein(1,1)']*1e6*enescal


ryain = inputFileDict['setup']['ryain']

_,Zeffin = helper.getCQLZeff(rho_pol = ryain)


neFunc = interp1d(ryain, n_e)
TiFunc = interp1d(ryain, T_i)
TeFunc = interp1d(ryain, T_e)
ZeffFunc = interp1d(ryain, Zeffin)
rya = np.ma.getdata(cql_nc.variables["rya"][:])
dvol = np.ma.getdata(cql_nc.variables["dvol"][:])

totalVol = np.sum(dvol)
volAvgne = np.sum(neFunc(rya)*dvol)/totalVol
volAvgTe = np.sum(TeFunc(rya)*dvol)/totalVol
volAvgTi = np.sum(TiFunc(rya)*dvol)/totalVol
volAvgZeff = np.sum(ZeffFunc(rya)*dvol)/totalVol

temp_peaking = (np.max(T_e)/volAvgTe + np.max(T_i)/volAvgTi)/2

aminor = 0.59
crossSectionArea = np.ma.getdata(cql_nc.variables["area"][:])[-1]/1e4

print(f'volAVgNe: {volAvgne}, volAvgTe: {volAvgTe}, volAvgZeff: {volAvgZeff}')
print(f'volAvgTi/volAvgTe ratio: {volAvgTi/volAvgTe}')
print(f'average temperature peaking: {temp_peaking}')
print(f'areal elongation = {crossSectionArea/(np.pi*aminor**2)}')

## get line averaged density:
rho_pol, ne = helper.getCQLne()
R_LFS = helper.convertRhopolToRmidplane(rhos = rho_pol, side = 'LFS')
R_HFS = helper.convertRhopolToRmidplane(rhos = rho_pol, side = 'HFS')

nes = np.concatenate([np.flip(ne), ne])
nes_mid = (nes[1:] + nes[:-1])/2
Rs = np.concatenate([np.flip(R_HFS), R_LFS])
dRs = Rs[1:] - Rs[:-1]

lineAvgne = np.sum(nes_mid*dRs)/np.sum(dRs)
print(f'line averaged density: {lineAvgne}')

lineAvgne1 = np.trapz(nes, x = Rs)/np.trapz(np.ones(len(nes)), x = Rs)
print(f'line averaged density: {lineAvgne1}')

fig,ax = plt.subplots()
ax.plot(Rs, nes)
ax.axhline(lineAvgne)
ax.axhline(lineAvgne1)
plt.show()
