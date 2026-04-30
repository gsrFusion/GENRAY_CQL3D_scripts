###
# 2D plots of the total and poloidal magnetic fields
###

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d, RectBivariateSpline
import matplotlib.cm as cm
import matplotlib
import netCDF4

import os, sys
from matplotlib.path import Path
from scipy.signal import find_peaks
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getGfileDict

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)
plt.rc('figure', titlesize = 16)

def plotB_pol(ax, targetDir,label):

    #get the variables we'll need for plotting
    rgrid = gfileDict["rgrid"]
    zgrid = gfileDict["zgrid"]
    B_zGrid = gfileDict["bzrz"]
    B_TGrid = gfileDict["btrz"]
    B_rGrid = gfileDict["brrz"]

    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS

    B_polstrength = np.sqrt(np.square(B_zGrid) + np.square(B_rGrid))

    B_polstrength = B_polstrength/np.max(B_polstrength)

    print(np.nanmax(B_polstrength))
    p2 = ax.pcolormesh(rgrid, zgrid, B_polstrength, shading = 'nearest',cmap='viridis', vmin=np.min(0), vmax = np.max(1))

    ax.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS

    ax.set_aspect('equal')
    #ax.set_title(label)
    

    ax.set_ylabel('Z (m)', labelpad = -10)
    ax.set_xlabel('R (m)')
    ax.set_ylim([np.min(zgrid),np.max(zgrid)])
    ax.set_xlim([np.min(rgrid),np.max(rgrid)])
    return p2

machine = 'NTPT'
fakeDevice = 'DIIID'
power = 2

PT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.PT05/{machine}_{fakeDevice}.PT05_n2.3Npara_180thgrill_{power}MW'
#NT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.NT05/{machine}_{fakeDevice}.NT05_n2.1Npara_180thgrill_{power}MW'
NT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/NTPT_AlcatorC.NT00/NTPT_AlcatorC.NT00'

targetDirs = [PT_target, NT_target]

labels= [r'$\delta = 0.5$', r'$\delta = -0.5$']

fig,axes = plt.subplots(1,2,figsize = (7,5.5))

for i, targetDir in enumerate(targetDirs):
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)
    pcm = plotB_pol(axes[i], targetDir, labels[i])  # Assume this returns a QuadMesh (e.g., from pcolormesh)

fig.subplots_adjust(right=0.85, wspace=0.25)

# Add a manually placed colorbar
cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
fig.colorbar(pcm, cax=cbar_ax).set_label(r'|B$_{\theta}$| (normalized)')

#fig.tight_layout()
plt.show()
