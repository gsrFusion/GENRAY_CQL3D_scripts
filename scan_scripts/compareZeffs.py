"""
Plots the ray traces and the RF power deposition density
"""
import numpy as np
import matplotlib.pyplot as plt

import os, sys
from scipy.signal import find_peaks
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import helperFunctions as helper
import netCDF4

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 14)
plt.rc('legend', fontsize = 12)

    
def main():

    shotNum = '203619.04135'
    stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum'
    
    t1 = f'{stem}'
    t2 = f'{stem}_2Zeff_lower'

    targetDirs = [t1,t2]#[f'{stem}_n01Efield', f'{stem}_n005Efield',f'{stem}', f'{stem}_0.005Efield', f'{stem}_0.01Efield']##
    colors = ['tab:blue','tab:orange']
    fig,ax = plt.subplots()
    ax.set_ylabel(r'$Z_{eff}$ ')
    ax.set_xlabel(r'Frequency (GHz)')

    for i,targetDir in enumerate(targetDirs):
        rya, Zeff = helper.getCQLZeff(targetDir = targetDir)
        ax.plot(rya, Zeff, lw = 2, color = colors[i])

    fig.tight_layout()
    plt.show()

main()
