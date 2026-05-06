###
# Plots the LCFS of several shots. May be useful to debugging/understanding differences
###

import matplotlib.pyplot as plt
import os, sys
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
                '/home/grantr/symlinks/genray_batch/NTPT_shots/NTPT_DIIID.147634NT/NTPT_DIIID.147634NT_n2.5Npara_-0.5grillHeight_1MW',
                '/home/grantr/symlinks/genray_batch/NTPT_shots/NTPT_DIIID.147634PT/NTPT_DIIID.147634PT_n2.5Npara_-0.5grillHeight_1MW',
                  ]

    fig, ax = plt.subplots(figsize = (5.25,7.1))
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    for i,targetDir in enumerate(targetDirs):
       
        print(f'targetDir: {targetDir}')
        gfileDict = getGfileDict.getGfileDict(targetDir = targetDir)
        xlim = gfileDict["xlim"] #R points of the wall
        ylim = gfileDict["ylim"] #Z points of the wall
        ax.plot(xlim, ylim, color = 'r', lw = 2, )#plot wall
        
        helper.drawFluxSurfaces(ax, targetDir = targetDir, rhosToPlot = [.2,.4,.6,.8,1], 
                     colors = colors[i], zBounds = None, limPath = None)

    ax.set_aspect('equal')
    ax.set_ylim(min(ylim)*1.05, max(ylim)*1.05)
    ax.set_xlim(min(xlim)*.95, max(xlim)*1.05)
    ax.set_ylabel("Z (m)")
    ax.set_xlabel("R (m)")

    fig.tight_layout()

    plt.show()

plotLFCSs()