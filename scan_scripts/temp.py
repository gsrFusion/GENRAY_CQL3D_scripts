import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os, sys
import matplotlib
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getInputFileDictionary
import helperFunctions as helper
import getGfileDict
from omfit_classes import omfit_eqdsk
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)


targetDirs = ['/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_147634.04565/Npara_thgrill_scan/DIIID_147634.04565_n2.5Npara_120thgrill_1MW',
            '/home/grantr/symlinks/genray_batch/NTPT_shots/NTPT_DIIID.147634PT05/NTPT_DIIID.147634PT05_n2.7Npara_180thgrill_1MW',
            ]

eqdsks = [omfit_eqdsk.OMFITgeqdsk(f'{targetDirs[0]}/g147634.04565'),
          omfit_eqdsk.OMFITgeqdsk(f'{targetDirs[1]}/gNTPT_DIIID.147634_0.5delta_tokamaker'),]

fig, ax = plt.subplots(figsize = (5.25,7.1))
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
#for j in range(len(times)):
#for j in range(len(shots)):
for i,targetDir in enumerate(targetDirs):
    #shot = shots[j]
    #time = times[j]
    #stem = f'/home/grantr/scratch/genray_batch/DIIID_shots/DIIID_{shot}.{time}/{shot}.{time}profiles'

   eqdsks[i].plot()



plt.show()