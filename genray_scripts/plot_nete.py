###
# Plots the electron and ion densities and temperatures according to the genray input file
# The plot_profs.py cql3d script is superior, but this can be used as a sanity check
###

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__);dname = os.path.dirname(abspath);os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__));parentdir = os.path.dirname(currentdir);sys.path.append(parentdir)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

import getInputFileDictionary
inputFileDict = getInputFileDictionary.getInputFileDictionary('genray_LH')

import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)

temtab = inputFileDict['temtab']['prof']
dentab = inputFileDict['dentab']['prof']

T_e = np.zeros(int(len(temtab)/3))
T_i = np.copy(T_e)
n_e = np.copy(T_e)
n_i = np.copy(T_e)

#read in the profiles
for i in range(0,len(T_e)):
    T_e[i] = temtab[3*i]
    T_i[i] = temtab[3*i+1]
    n_e[i] = dentab[3*i]
    n_i[i] = dentab[3*i+1]

tescal = inputFileDict['plasma']['temp_scale(1)']
tiscal = inputFileDict['plasma']['temp_scale(2)']
nescal = inputFileDict['plasma']['den_scale(1)']
niscal = inputFileDict['plasma']['den_scale(2)']

#apply any scalings
T_e = T_e*tescal
T_i = T_i*tiscal
n_e = n_e*nescal
n_i = n_i*niscal
rhosHelper = np.arange(1,len(T_e)+1,1)
rhos = (rhosHelper-1)/(len(T_e)-1)

fig, ax = plt.subplots(figsize = (8*.8,6*.8))

ax.plot(rhos, T_e, label = r'$T_e$', linewidth = 3, color = 'g')
ax.plot(rhos, T_i, label = r'$T_i$', linewidth = 3, color = 'b')
ax.plot(rhos, np.array(n_e/1e19), label = r'$n_e$', linestyle = 'dashed', color = 'k', linewidth = 3)
ax.plot(rhos, np.array(n_i/1e19), label = r'$n_i$', linestyle = 'dotted', color = 'k', linewidth = 3)

ax.legend(loc = 'best', fontsize = 20, ncol=2)

ax.set_xlabel(r'$\rho_{pol}$', fontsize = 20)
ax.set_xlim([0,1])
ax.set_ylim([0,8])
ax.set_ylabel(r'n$_j$ (10$^{19}$/m$^{-3}$), T$_j$ (keV)')
ax.grid()

fig.suptitle(f'GENRAY {shotNum} Profiles\n enescal = {nescal}, Tescal = {tescal}, Tiscal = {tiscal}')
fig.tight_layout(rect=[0, 0.0, 1, 1.05])
plt.show()

