"""
Plots the n_para at which significant electron landau damping occurs
This code can take up to ~30s to run. There's probably a better way to do this plotting than I came up with
This method does not do a great job of plotting near the xpoints since the flux surfaces are quite far apart there
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d
import matplotlib.cm as cm
import matplotlib

import os,sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import getInputFileDictionary
import getGfileDict
import helperFunctions as helper
import getTargetInfo


cqlInputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d')
gfileDict = getGfileDict.getGfileDict()
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)

#below are the variables required to work with the magnetic field
rgrid = gfileDict["rgrid"]
zgrid = gfileDict["zgrid"]
magAxisR = gfileDict['rmaxis'] 
magAxisZ = gfileDict['zmaxis'] 

#relevant variables to find the normalized poloidal flux
psirz = gfileDict["psirz"]
psi_mag_axis = gfileDict["ssimag"]
psi_boundary = gfileDict["ssibdry"]
    
psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)
#interpolated function for poloidal flux
psirzNormFunc = interp2d(rgrid, zgrid, psirzNorm.T)

rInterp = np.linspace(np.min(rgrid), np.max(rgrid), 200)
zInterp = np.linspace(np.min(zgrid), np.max(zgrid), 200)#we need to restrict the Z so that it plots flux surfaces inside the LCFS and not in the divertor
psi_NrzInterp = interp2d(rgrid,zgrid, psirzNorm, kind = 'linear')(rInterp, zInterp)

ryain, Tein = helper.getCQLTe()

TeInterpFunc = interp1d(ryain, Tein, kind = 'linear', bounds_error = False, fill_value = np.nan)
Tes = TeInterpFunc(np.sqrt(psi_NrzInterp))

dampingNs = 6.4/np.sqrt(Tes)

fig,ax = plt.subplots(figsize = (5.25,6.5))

toPlot = dampingNs
toPlot[psi_NrzInterp > 1] = np.nan
toPlot[psi_NrzInterp < 0] = np.nan

p2 = ax.pcolormesh(rInterp, zInterp, toPlot,shading = 'nearest',cmap='viridis', vmin=np.nanmin(toPlot), vmax = 5)#np.nanmax(toPlot))
cbar = fig.colorbar(p2, ax = ax, shrink = .5, pad = .01)
cbar.set_label(r"N$_{||, ELD}$")
ax.set_aspect('equal')
ax.set_ylabel('Z (m)', labelpad = -10)
ax.set_xlabel('R (m)')
ax.set_title(f'Shot {shotNum}')
plt.show()
