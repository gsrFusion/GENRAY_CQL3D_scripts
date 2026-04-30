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

import getGfileDict
import helperFunctions as helper

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

def plotLFCSs():
    
    targetDirs = [
                '/home/grantr/symlinks/genray_batch/NTPT_shots/NTPT_DIIID.147634PT05Test/NTPT_DIIID.147634PT05Test_n2.8Npara_-0.25grillHeight_1MW',
                '/home/grantr/symlinks/genray_batch/NTPT_shots/NTPT_DIIID.147634PT05/NTPT_DIIID.147634PT05_n2.8Npara_-0.25grillHeight_1MW',
                  ]

    fig, ax = plt.subplots(figsize = (5.25,7.1))
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    #for j in range(len(times)):
    #for j in range(len(shots)):
    for i,targetDir in enumerate(targetDirs):
        #shot = shots[j]
        #time = times[j]
        #stem = f'/home/grantr/scratch/genray_batch/DIIID_shots/DIIID_{shot}.{time}/{shot}.{time}profiles'

        #targetDir = f'{stem}/DIIID_{shot}.{time}_{shot}.{time}profiles_n{2.7}Npara_300kW'
        print(f'targetDir: {targetDir}')
        gfileDict = getGfileDict.getGfileDict(targetDir = targetDir)
        xlim = gfileDict["xlim"] #R points of the wall
        ylim = gfileDict["ylim"] #Z points of the wall
        rbbbs = gfileDict["rbbbs"] #R points of the LCFS
        zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
        ax.plot(xlim, ylim, color = 'r', lw = 2, )#plot wall
        
        helper.drawFluxSurfaces(ax, gfileDict = gfileDict, rhosToPlot = [.2,.4,.6,.8,1], 
                     colors = colors[i], zBounds = None, limPath = None)

    #ax.legend(loc = 'upper right')

    ax.set_aspect('equal')
    ax.set_ylim(min(ylim)*1.05, max(ylim)*1.05)
    ax.set_xlim(min(xlim)*.95, max(xlim)*1.05)
    ax.set_ylabel("Z (m)")
    ax.set_xlabel("R (m)")

    #ax.set_title(rf'Shot {shot}')
    fig.tight_layout()

    plt.show()

plotLFCSs()