"""
Plots the minimum propagating n_para according to the accessibility condition inside the LCFS
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)

import os,sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import getGfileDict
gfileDict = getGfileDict.getGfileDict()
import helperFunctions as helper
import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
import netCDF4

genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')#netCDF4.Dataset(f'{targetDir}/genray.nc','r')
wr  = genray_ece_nc.variables['wr_em_nc'][:] #major radius of the ray at each point along the trace, in m
wz  = genray_ece_nc.variables['wz_em_nc'][:]
wfreq_nc = genray_ece_nc.variables['wfreq_nc'][:]

#below are the variables required to work with the magnetic field
rgrid = gfileDict["rgrid"]
zgrid = gfileDict["zgrid"]
magAxisR = gfileDict['rmaxis'] 
magAxisZ = gfileDict['zmaxis'] 
B_zGrid = gfileDict["bzrz"]
B_TGrid = gfileDict["btrz"]
B_rGrid = gfileDict["brrz"]

xlim = gfileDict["xlim"] #R points of the wall
ylim = gfileDict["ylim"] #Z points of the wall

#relevant variables to find the normalized poloidal flux
psirz = gfileDict["psirz"]
psi_mag_axis = gfileDict["ssimag"]
psi_boundary = gfileDict["ssibdry"]
    
psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)
#interpolated function for poloidal flux
psirzNormFunc = interp2d(rgrid, zgrid, psirzNorm.T)

rInterp = np.linspace(np.min(rgrid), np.max(rgrid), 200)
zInterp = np.linspace(np.min(zgrid), np.max(zgrid), 200)
psi_NrzInterp = interp2d(rgrid,zgrid, psirzNorm, kind = 'linear')(rInterp, zInterp)

Bstrength = np.sqrt(np.square(B_zGrid) + np.square(B_TGrid) + np.square(B_rGrid))
interped_B = interp2d(rgrid,zgrid,Bstrength, kind = 'linear')(rInterp, zInterp)

ryain, nein = helper.getCQLne()
ryain, nDin = helper.getCQLnD()

neInterpFunc = interp1d(ryain, nein, kind = 'linear', bounds_error = False, fill_value = np.nan)
nDInterpFunc = interp1d(ryain, nDin, kind = 'linear', bounds_error = False, fill_value = np.nan)

nes = neInterpFunc(np.sqrt(psi_NrzInterp))
nDs = nDInterpFunc(np.sqrt(psi_NrzInterp))

m_e = 9.109e-31 #electron mass
m_D = 3.343e-27 #deuteron mass
q = 1.6e-19 #elementary charge
eps_0 = 8.85e-12 #permitivity of free space

w_pes = np.sqrt(nes*q**2/(m_e*eps_0))
w_pDs = np.sqrt(nDs*q**2/(m_D*eps_0))
w = 2*np.pi*4.6e9
W_ces = q*interped_B/m_e

w_R = W_ces/2 + np.sqrt(W_ces**2/4 + w_pes**2)


rayNum = 14
freq = wfreq_nc[rayNum]
Rs = wr[rayNum]
Zs = wz[rayNum]

#w_R[w_R >= freq*2*np.pi*1e9] = np.nan

fig,ax = plt.subplots(figsize = (4.25,7.1))


toPlot = w_R/(freq*2*np.pi*1e9)
toPlot[toPlot>=1] = np.nan
toPlot[psi_NrzInterp > 1] = np.nan
toPlot[psi_NrzInterp < 0] = np.nan

ax.plot(Rs[Rs > 1], Zs[Rs > 1], lw = 3, color = 'r')
helper.drawFluxSurfaces(ax, colors = 'k')

p2 = ax.pcolormesh(rInterp, zInterp, toPlot,shading = 'nearest',cmap='viridis', vmin=0, vmax = 1)
cbar = fig.colorbar(p2, ax = ax, shrink = .5, pad = .01)
cbar.set_label(r"$\omega_R / \omega_{ECE}$")
ax.plot(xlim, ylim, color = 'grey', lw = 3)#plot wall

ax.set_aspect('equal')
ax.set_ylabel('Z (m)', labelpad = -10)
ax.set_xlabel('R (m)')
ax.set_title(f'Shot {shotNum}')
fig.tight_layout()
plt.show()