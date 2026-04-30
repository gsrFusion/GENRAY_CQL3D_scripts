###
# Plots the operational window between N_||,acc and N_||,ELD in which the slow LH wave can propagate
# this calculation is done at the midplane and spans from the HFS to the LFS
###

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d, RectBivariateSpline

import os,sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import getInputFileDictionary
cqlInputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d')
import getGfileDict
gfileDict = getGfileDict.getGfileDict()

m_e = 9.109e-31 #electron mass
m_D = 3.343e-27 #deuteron mass
q = 1.6e-19 #elementary charge
eps_0 = 8.85e-12 #permitivity of free space

def getW2_pe(n):#electron plasma frequency squared
    return n*q**2/(m_e*eps_0)
def getW2_pD(n):#deuteron plasma frequency squared
    return n*q**2/(m_D*eps_0)
def getW_ce(B):#electron cyclotron frequency
    return -q*B/m_e

#returns the accessibility criterion
def getAccess(n, B):
    w2_pD = getW2_pD(n); w2_pe = getW2_pe(n); w_ce = getW_ce(B)
    w = 4.6e9*2*np.pi #wave's angular frequency

    return np.sqrt(1 - w2_pD/w**2 + w2_pe/w_ce**2) + np.sqrt(w2_pe)/np.abs(w_ce)

#below are the variables required to work with the magnetic field
rgrid = gfileDict["rgrid"]

LCFS_mask = (rgrid < np.max(gfileDict['rbbbs'])*(rgrid > np.min(gfileDict['rbbbs'])))

R_insideLCFS = rgrid[np.where(LCFS_mask)]
zgrid = gfileDict["zgrid"]
magAxisZ = gfileDict['zmaxis'] 
btmid = gfileDict['btmid'] #toroidal field at the magnetic midplane
bzmid = gfileDict['bzmid'] #vertical field at the magnetic midplane

B_mid = np.sqrt(btmid**2 + bzmid**2)
B_mid_insideLCFS = B_mid[np.where(LCFS_mask)]

psirz = gfileDict["psirz"]
psi_mag_axis = gfileDict["ssimag"]
psi_boundary = gfileDict["ssibdry"]
    
psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)
#interpolated function for poloidal flux
psirzNormFunc = RectBivariateSpline(zgrid,rgrid,psirzNorm)

ryain = cqlInputFileDict["setup"]["ryain"]
Tes = cqlInputFileDict["setup"]["tein"]
ns = cqlInputFileDict["setup"]["enein(1,1)"]*1e6

tescal = 1
try:
    tescal = cqlInputFileDict["setup"]["tescal"]
except:
    pass
TeInterpFunc = interp1d(ryain**2, Tes*tescal, kind = 'cubic')
neInterpFunc = interp1d(ryain**2, ns, kind = 'cubic')

psiMidplane = psirzNormFunc(magAxisZ,R_insideLCFS).flatten()

plt.rc('ytick', labelsize = 16)
plt.rc('xtick', labelsize = 16)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)
plt.rc('legend',fontsize=20)
plt.rcParams['axes.ymargin'] = 0
plt.rcParams['axes.xmargin'] = 0
#plt.rcParams['figure.dpi'] = 200

dampingCond = 5.8/np.sqrt(TeInterpFunc(psiMidplane))
accessCond = getAccess(neInterpFunc(psiMidplane), B_mid_insideLCFS)
print(f'min of damping condition: {np.min(dampingCond)}')
fig, ax = plt.subplots(figsize = (5.25*1.4,2.5*1.4))

psiMidplane_test = psiMidplane
psiMidplane_test[R_insideLCFS < gfileDict['rmaxis']] *= -1

ax.plot(R_insideLCFS, accessCond, color = 'firebrick', linewidth = 2)
ax.plot(R_insideLCFS, dampingCond, color = 'darkblue', linewidth = 2)

ax.fill_between(R_insideLCFS,accessCond,dampingCond,color = 'forestgreen')

ax.fill_between(R_insideLCFS,accessCond,np.min(accessCond),color = 'w')
ax.fill_between(R_insideLCFS,dampingCond,np.max(dampingCond),color = 'w')
#"""
#ax.text(1.07, 0.8, 'HFS', fontsize = 16)
#ax.text(2.2, 0.8, 'LFS',  fontsize = 16)

upperPlottingLim = 3#6
"""
midR = (np.max(R_insideLCFS) + np.min(R_insideLCFS))/2 + .05
inaccessZTop = .3/(2.5*1.4) * upperPlottingLim + np.min(accessCond)

ax.text(midR,inaccessZTop, 'Inaccessible', fontsize = 18, color = 'firebrick', 
         bbox = dict(facecolor = 'w', alpha = 0, edgecolor = 'w'),ha='center', va='center')
         
ax.text(midR,.9*upperPlottingLim, 'Strong Damping', fontsize = 18, color = 'darkblue', 
         bbox = dict(facecolor = 'w', alpha = 0, edgecolor = 'w'),ha='center', va='center')         
"""
ax.set_ylim(bottom = 1, top = upperPlottingLim)
#"""
ax.set_ylabel(r"$|N_{\parallel}|$")
ax.set_xlabel(r"R$_{midplane}$ (m)")

ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=True,
                     bottom=True, top=False, left=True, right=True)

fig.tight_layout()
plt.show()



