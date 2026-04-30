import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 22)
plt.rc('legend', fontsize = 14)

import os, sys
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getInputFileDictionary
import getTargetInfo
import getGfileDict
targetDir = getTargetInfo.getTargetDir()
machine = getTargetInfo.getMachine()
genray_in = getInputFileDictionary.getInputFileDictionary('genray')
gfileDict = getGfileDict.getGfileDict()

m_e = 9.109e-31
e = 1.602e-19
c= 2.99e8
gamma = (100*1.602e-16)/(m_e*c**2) + 1
print(f'gamma: {gamma}')
rgrid = gfileDict["rgrid"]

LCFS_mask = (rgrid < np.max(gfileDict['rbbbs'])*(rgrid > gfileDict['rmaxis']))

R_insideLCFS = rgrid[np.where(LCFS_mask)]
zgrid = gfileDict["zgrid"]
magAxisZ = gfileDict['zmaxis'] 
btmid = gfileDict['btmid'] #toroidal field at the magnetic midplane
bzmid = gfileDict['bzmid'] #vertical field at the magnetic midplane

B_mid = np.sqrt(btmid**2 + bzmid**2)
B_mid_insideLCFS = B_mid[np.where(LCFS_mask)]

f_ce = e*B_mid_insideLCFS/(m_e *2*np.pi*1e9)

fig,ax = plt.subplots()

#ax.plot(R_insideLCFS, f_ce, lw = 2, color = 'tab:blue', label = '1st harmonic')
#ax.plot(R_insideLCFS, f_ce/gamma, lw = 2, color = 'tab:blue', linestyle = 'dashed', label = 'relativistic 1st harmonic')

ax.plot(R_insideLCFS, 2*f_ce, lw = 2, color = 'k', label = r'2nd harmonic ($E_{kin} = 0$ keV)')
ax.plot(R_insideLCFS, 2*(f_ce/gamma), lw = 2, color = 'k', linestyle = 'dashed', label = r'2nd harmonic ($E_{kin} = 100$ keV)')

#ax.plot(R_insideLCFS, 3*f_ce, lw = 2, color = 'tab:green', label = '3rd harmonic')
#ax.plot(R_insideLCFS, 3*(f_ce/gamma), lw = 2, color = 'tab:green', linestyle = 'dashed', label = 'relativistic 3rd harmonic')

ax.set_ylabel('Frequency (GHz)')
ax.set_xlabel('Major radius (m)')
ax.set_title(f'Shot 203912', loc = 'right')
ax.set_ylim(bottom = 65, top = 120)

ax.legend()
fig.tight_layout()
plt.show()
