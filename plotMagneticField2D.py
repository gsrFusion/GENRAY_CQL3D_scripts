###
# 2D plots of the total and poloidal magnetic fields
###

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d, RectBivariateSpline
import matplotlib.cm as cm
import matplotlib
import netCDF4
import getGfileDict
import helperFunctions as helper
from matplotlib.path import Path
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)
plt.rc('figure', titlesize = 16)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()
print(targetDir)

gfileDict = getGfileDict.getGfileDict()

#get the variables we'll need for plotting
rgrid = gfileDict["rgrid"]
zgrid = gfileDict["zgrid"]
B_zGrid = gfileDict["bzrz"]
B_TGrid = gfileDict["btrz"]
B_rGrid = gfileDict["brrz"]

xlim = gfileDict["xlim"] #R points of the wall
ylim = gfileDict["ylim"] #Z points of the wall
rbbbs = gfileDict["rbbbs"] #R points of the LCFS
zbbbs = gfileDict["zbbbs"] # Z points of the LCFS

B_polstrength = np.sqrt(np.square(B_zGrid) + np.square(B_rGrid))
B_totstrength = np.sqrt(np.square(B_zGrid) + np.square(B_rGrid) + np.square(B_TGrid))
if machine == 'FENIX':
    B_polstrength = B_polstrength[30:-30, 20:-20]
    B_totstrength = B_totstrength[30:-30, 20:-20]
    zgrid = zgrid[30:-30]
    rgrid = rgrid[20:-20]#rgrid[20:-20]
    B_zGrid = B_zGrid[30:-30, 20:-20]
    B_rGrid = B_rGrid[30:-30, 20:-20]

from scipy.interpolate import interp2d, interp1d, RectBivariateSpline
BT_func = RectBivariateSpline(zgrid,rgrid,B_TGrid)
B_func = RectBivariateSpline(zgrid,rgrid,B_totstrength)
print(f'on-axis B: {BT_func(gfileDict["zmaxis"],gfileDict["rmaxis"])}')
print(f'bcentr: {gfileDict["bcentr"]}')
print(f'rmaxis: {gfileDict["rmaxis"]}, rzero: {gfileDict["rzero"]}')

print(f'B_func(.0269, 1.9047): {B_func(.0269, 1.9047)}')

"""
fig,axes = plt.subplots(1,2,figsize = (7,5.5))

#plot total magnetic field
#ratio1 = B_totstrength
p1 = axes[0].pcolormesh(rgrid, zgrid, B_totstrength,shading = 'nearest',cmap='viridis', vmin=np.min(0), vmax = np.max(3))
axes[0].set_aspect('equal')
axes[0].set_title(r'B$_{tot}$')
axes[0].plot(xlim, ylim, 'r', lw = 2)#plot wall
axes[0].plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS
cbar = fig.colorbar(p1, ax = axes[0], shrink = .5, pad = .01)
#cbar.set_label(r"B$_{total}$ (T)")
axes[0].set_ylabel('Z (m)', labelpad = -10)
axes[0].set_xticks([1,1.5,2,2.5])
axes[0].set_xlabel('R (m)')
axes[0].set_ylim([np.min(zgrid),np.max(zgrid)])
axes[0].set_xlim([np.min(rgrid),np.max(rgrid)])

midplaneIndex = int(len(zgrid)/2)

#axes[1].plot(rgrid, B_polstrength[midplaneIndex,:])
p2 = axes[1].pcolormesh(rgrid, zgrid, B_polstrength,shading = 'nearest',cmap='viridis', vmin=np.min(0), vmax = np.max(.4))
axes[1].set_aspect('equal')
axes[1].set_title(r'B$_{pol}$')
axes[1].plot(xlim, ylim, 'r', lw = 2)#plot wall
axes[1].plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS
cbar = fig.colorbar(p2, ax = axes[1], shrink = .5, pad = .01)
cbar.set_label(r"B$_{pol}$")
axes[1].set_ylabel('Z (m)', labelpad = -10)
axes[1].set_xlabel('R (m)')
#axes[1].set_ylim([np.min(zgrid),np.max(zgrid)])
#axes[1].set_xlim([np.min(rgrid),np.max(rgrid)])
"""
f2_ce = 27.9906*B_totstrength*2
minf2 = 83.5
maxf2 = 84.5
B_totstrength[f2_ce < minf2] = np.nan
B_totstrength[f2_ce > maxf2] = np.nan



fig,ax = plt.subplots(figsize = (5,6))
p1 = ax.pcolormesh(rgrid, zgrid, B_totstrength,shading = 'nearest',cmap='viridis', vmin=0, vmax = 4)
#p1 = ax.pcolormesh(rgrid, zgrid, B_polstrength,shading = 'nearest',cmap='viridis', vmin=0, vmax = 0.35)
ax.set_aspect('equal')
#ax.set_title('thermal resonance layer for \n' + r'$2 \Omega_e = $ 83-84 GHZ')
ax.plot(xlim, ylim, 'grey', lw = 3)#plot wall
#ax.scatter(xlim, ylim)#plot wall
lim = np.array([xlim, ylim]).T
limiterPath  = Path(lim)
#helper.drawFluxSurfaces(ax, colors = 'k', rhosToPlot = [.33,.66, 1], limPath = limiterPath)
helper.drawFluxSurfaces(ax, colors = 'k', rhosToPlot = [.43], limPath = limiterPath)
#ax.plot(rbbbs, zbbbs, 'k', lw = 2)#plot LCFS
cbar = fig.colorbar(p1, ax = ax, pad = .01,label='BT')
cbar.set_label(r"B$_{\theta}$ (T)")
ax.set_ylabel('Z (m)', labelpad = -10)
ax.set_xticks([1,1.5,2,2.5])
ax.set_xlabel('R (m)')
ax.set_ylim([np.min(zgrid),np.max(zgrid)])
ax.set_xlim([np.min(rgrid),np.max(rgrid)])
#"""

#fig.suptitle(f'Shot {shotNum}')
fig.tight_layout()
#plt.savefig('147634_Bpol.jpeg', dpi=300)
plt.show()
