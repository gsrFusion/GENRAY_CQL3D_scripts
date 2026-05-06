###
# Makes a RZ plot showing the difference in Bpol between two plasmas
###

import numpy as np
import matplotlib.pyplot as plt
import os, sys
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

machine = 'NTPT'
fakeDevice = 'DIIID'
power = 1

PT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.147634PT/{machine}_{fakeDevice}.147634PT_n2.5Npara_0.0grillHeight_{power}MW'
NT_target = f'/home/grantr/symlinks/genray_batch/NTPT_shots/{machine}_{fakeDevice}.147634NT/{machine}_{fakeDevice}.147634NT_n2.5Npara_0.0grillHeight_{power}MW'

gfileDict_PT = getGfileDict.getGfileDict(targetDir = PT_target)
gfileDict_NT = getGfileDict.getGfileDict(targetDir = NT_target)

#get the variables we'll need for plotting
rgrid_PT = gfileDict_PT["rgrid"]
zgrid_PT = gfileDict_PT["zgrid"]
B_zGrid_PT = gfileDict_PT["bzrz"]
B_TGrid_PT = gfileDict_PT["btrz"]
B_rGrid_PT = gfileDict_PT["brrz"]

rbbbs_PT = gfileDict_PT["rbbbs"] #R points of the LCFS
zbbbs_PT = gfileDict_PT["zbbbs"] # Z points of the LCFS

rgrid_NT = gfileDict_NT["rgrid"]
zgrid_NT = gfileDict_NT["zgrid"]
B_zGrid_NT = gfileDict_NT["bzrz"]
B_TGrid_NT = gfileDict_NT["btrz"]
B_rGrid_NT = gfileDict_NT["brrz"]

rbbbs_NT = -gfileDict_NT["rbbbs"] #R points of the LCFS
zbbbs_NT = -gfileDict_NT["zbbbs"] # Z points of the LCFS

B_polstrength_PT = np.sqrt(np.square(B_zGrid_PT) + np.square(B_rGrid_PT))
B_polstrength_NT = np.sqrt(np.square(B_zGrid_NT) + np.square(B_rGrid_NT))

fig,ax = plt.subplots(figsize = (5.25,6.5))

p2 = ax.pcolormesh(rgrid_PT, zgrid_PT, B_polstrength_PT - B_polstrength_NT,shading = 'nearest',cmap='bwr', vmin=np.min(-.2), vmax = np.max(.2))

ax.plot(rbbbs_PT, zbbbs_PT, 'k', lw = 1.5)#plot LCFS
ax.plot(rbbbs_NT, zbbbs_NT, 'r', lw = 1.5)#plot LCFS


ax.set_aspect('equal')
ax.set_title(r'B$_{pol, PT}$ - B$_{pol, NT}$')
cbar = fig.colorbar(p2, ax = ax, shrink = .5, pad = .01)
ax.set_ylabel('Z (m)', labelpad = -10)
ax.set_xlabel('R (m)')
ax.set_ylim([np.min(zgrid_PT),np.max(zgrid_PT)])
ax.set_xlim([np.min(rgrid_PT),np.max(rgrid_PT)])

fig.tight_layout()
plt.show()
