"""
Plots the propagation domain for the LH waves.
I have found this not overly useful and the below code is rather deprecated
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
from matplotlib.collections import LineCollection
import netCDF4
import os, sys
from numpy import sqrt
import matplotlib

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

currentDir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentDir)

sys.path.append(parentdir)
import getGfileDict
import getInputFileDictionary
import helperFunctions as helper
import getTargetInfo

gfileDict = getGfileDict.getGfileDict()
inputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d')
genray_in = getInputFileDictionary.getInputFileDictionary('genray')
n_para_f = (genray_in["grill"]["anmax(1)"] + genray_in["grill"]["anmin(1)"])/2
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')

q = 1.602e-19    
m_e = 9.109e-31
m_D = 3.343e-27
eps_0 = 8.854e-12
c = 2.99792e8

w = cqlrf_nc.variables["freqcy"][:]*2*np.pi

tiscal = 1
enescal = 1
tescal = 1

plt.rc('xtick', labelsize = 17)
plt.rc('ytick', labelsize = 17)
plt.rc('axes', labelsize = 19)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 16)

rgrid = gfileDict["rgrid"]
zgrid = gfileDict["zgrid"]

def getBComponents(r, Z):
    B_zGrid = gfileDict["bzrz"]
    B_TGrid = gfileDict["btrz"]
    B_rGrid = gfileDict["brrz"]

    B_THFSMid = RectBivariateSpline(zgrid,rgrid,B_TGrid)(Z, r).flatten()
    B_zHFSMid = RectBivariateSpline(zgrid,rgrid,B_zGrid)(Z, r).flatten()
    B_rHFSMid = RectBivariateSpline(zgrid,rgrid,B_rGrid)(Z, r).flatten()

    #B_THFSMid = interp2d(rgrid,zgrid, B_TGrid, kind = 'linear')(r, Z)
    #B_zHFSMid = interp2d(rgrid,zgrid, B_zGrid, kind = 'linear')(r, Z)
    #B_rHFSMid = interp2d(rgrid,zgrid, B_rGrid, kind = 'linear')(r, Z)

    B_polMid = np.sqrt(B_zHFSMid**2 + B_rHFSMid**2)

    return B_THFSMid, B_polMid

def getPsiNormGrid():
    psirz = gfileDict["psirz"]
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]

    return (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

def getTeProf(rhos):
    global tescal
    ryain = inputFileDict['setup']['ryain']
    try:
        tescal = inputFileDict['setup']['tescal']
    except:
        pass

    T_e = inputFileDict['setup']['tein']*tescal
    
    return interp1d(ryain, T_e)(rhos)

def getTiProf(rhos):
    global tiscal
    ryain = inputFileDict['setup']['ryain']
    try:
        tiscal = inputFileDict['setup']['tiscal']
    except:
        pass

    T_i = inputFileDict['setup']['tiin']*tiscal
    
    return interp1d(ryain, T_i)(rhos)

def getDenProf(rhos):
    global enescal
    ryain = inputFileDict['setup']['ryain']
    try:
        enescal = inputFileDict['setup']['enescal']
    except:
        pass

    nProf = inputFileDict['setup']['enein(1,1)']*1e6*enescal
    
    return interp1d(ryain, nProf)(rhos)

def getVtheProfile_rho(TeProf):
    vthe = sqrt(2*TeProf*1.60218e-16)/sqrt(m_e)
    
    return vthe

#returns the plasma frequency for a given species and a given density
def getWp(species, denProf):
    stemValue = np.sqrt(q**2*denProf/eps_0)
        
    if species == 'e':
        return stemValue/np.sqrt(m_e)
    elif species == 'D':
        return stemValue/np.sqrt(m_D)
    else:
        raise Exception('Invalid species')

#returns the cyclotron frequency for a given species and given B
def getWc(species, Bprof):
    stemValue = q*Bprof
    
    if species == 'e':
        return stemValue*(-1/m_e)
    elif species == 'D':
        return stemValue*(1/m_D)
    else:
        raise Exception('Invalid species')
    
#returns S, D, R, L, and P for the input magnetic field and density
def getStixComponents(Bprof, denProf):
    w_ce = getWc('e', Bprof)
    w_cD = getWc('D', Bprof)
    
    w_pe = getWp('e', denProf)
    w_pD = getWp('D', denProf)

    R = 1- (w_pe**2/(w*(w + w_ce))) - (w_pD**2/(w*(w + w_cD)))
    L = 1- (w_pe**2/(w*(w - w_ce))) - (w_pD**2/(w*(w - w_cD)))
    P = 1 - (w_pe**2/w**2) - (w_pD**2/w**2)
    
    S = .5*(R+L); D = .5*(R - L)
    
    return S, D, R, L, P

def addRays(ax):
    delpwr= cqlrf_nc.variables["delpwr"][:][:lastForwardLobeRay,:] #power in the ray at each point
    spsi = cqlrf_nc.variables["spsi"][:][:lastForwardLobeRay,:] #radial like variable
    maxDelPwrPlot = .8 #what portion of ray power must have been damped before we stop plotting that ray

    norm = plt.Normalize(0, 1)
    #plot the ray using a LineCollection which allows the colormap to be applied to each ray
    for ray in range(len(wnpar)):
        delpwr[ray,:] = delpwr[ray,:]/delpwr[ray,0] #normalize the ray power to that ray's starting power
        mostPowerDep = helper.findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]) #find the index of the last ray point we want to plot

        points = np.array([spsi[ray][:mostPowerDep], np.abs(wnpar[ray][:mostPowerDep])]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Create a continuous norm to map from datline co points to colors
        lc = LineCollection(segments, norm = norm,cmap=plt.cm.jet, zorder = 5)
        # Set the values used for colormapping
        lc.set_array(delpwr[ray][:mostPowerDep])
        lc.set_linewidth(1)
        ax.add_collection(lc)


psirzNorm = getPsiNormGrid()

R_HFSMid = np.linspace(np.min(rgrid), gfileDict['rmaxis'], 1000)
psiNormHFSMid = interp2d(rgrid,zgrid, psirzNorm, kind = 'cubic')(R_HFSMid, 0)
mask = np.where((psiNormHFSMid <= 1)*(R_HFSMid > gfileDict['rmaxis']/2))

R_HFSMid = R_HFSMid[mask]

R_HFSMid = np.linspace(np.min(R_HFSMid), gfileDict['rmaxis'], 500)
psiNormHFSMid = interp2d(rgrid,zgrid, psirzNorm, kind = 'cubic')(R_HFSMid, 0)

psiNormHFSMid = psiNormHFSMid[psiNormHFSMid <= 1]
rhoHFS = np.sqrt(psiNormHFSMid)

B_Ts, B_pols = getBComponents(R_HFSMid, gfileDict['zmaxis'])
B_tots = np.sqrt(B_Ts**2 + B_pols**2)

gamma = B_pols/B_tots

denProf = getDenProf(rhoHFS)
TeProf = getTeProf(rhoHFS)
TiProf = getTiProf(rhoHFS)

vthe = getVtheProfile_rho(TeProf)

stix_S, stix_D, stix_R, stix_L, stix_P = getStixComponents(B_tots, denProf)

rya = cql_nc.variables["rya"][:]

wnpar = cqlrf_nc.variables['wnpar'][:]
startingNpara = wnpar[:,0]

lastForwardLobeRay = len(wnpar)
wnpar = wnpar[:lastForwardLobeRay,:]
Nphi = genray_nc.variables['wn_phi'][:][:lastForwardLobeRay,:]#np.copy(disp_cdf.variables['nphi'].data)
wr = cqlrf_nc.variables['wr'][:][:lastForwardLobeRay,:]*.01 #convert to m

avgnphi = np.mean((w/c)*Nphi[:,0]*wr[:,0])
R0 = cql_nc.variables["radmaj"][:]

trhs = R_HFSMid

#return propagation domain function
def f(npar,index):

    npar2 = npar*npar
    P0 = stix_P[index]*((npar2-stix_R[index])*(npar2-stix_L[index]))
    P2 = (npar2-stix_S[index])*(stix_P[index]+stix_S[index])+stix_D[index]**2
    P4 = stix_S[index]

    tmp = (npar*np.sqrt(1-gamma[index]**2) - (c/w)*avgnphi/R_HFSMid[index])**2 - (gamma[index]**2)*(-P2+np.sqrt(P2**2-4*P0*P4))/(2*P4)

    return tmp

#brain dead root solver
def getNparRoots(index, func):
    nparas = -1*np.linspace(0,10,1000)
    values = func(nparas, index)
    if np.nanmin(values) > 0 or np.nanmax(values) < 0:
        return []
    
    if np.isnan(values).any() > 0:
        lastNanIndex = np.where(np.isnan(values))[0][-1]

        values = values[lastNanIndex+1:]
        nparas = nparas[lastNanIndex+1:]

    roots = []

    if np.max(values) < 0:
        roots.append(np.max(nparas))
    else:
        lessThanZero = (values < 0)
        for j in range(len(values)-1):
            if lessThanZero[j] != lessThanZero[j+1]:
                roots.append(nparas[j])
        if lessThanZero[0] and lessThanZero[1]:
            roots.append(nparas[0])

    if len(roots) == 1:
        roots.append(-10)

    assert len(roots) <= 2

    return roots

def g(npar,index):
    npar2 = npar*npar
    P0 = stix_P[index]*((npar2-stix_R[index])*(npar2-stix_L[index]))
    P2 = (npar2-stix_S[index])*(stix_P[index]+stix_S[index])+stix_D[index]**2
    P4 = stix_S[index]
    tmp = P2**2-4*P0*P4
    return tmp

def n_estat(ipsi):
    nup = (avgnphi/R_HFSMid[ipsi])/(np.sqrt(1-gamma[ipsi]**2)-np.sqrt(-stix_P[ipsi]/stix_S[ipsi])*gamma[ipsi])
    nlo = (avgnphi/R_HFSMid[ipsi])/(np.sqrt(1-gamma[ipsi]**2)+np.sqrt(-stix_P[ipsi]/stix_S[ipsi])*gamma[ipsi])
    
    return [nup,nlo]

def getAccess():
    w_ce = getWc('e', B_tots)
    w_cD = getWc('D', B_tots)
    
    w2_pe = getWp('e', denProf)**2
    w2_pD = getWp('D', denProf)**2

    return np.sqrt(1 - w2_pD/w**2 + w2_pe/w_ce**2) + np.sqrt(w2_pe)/np.abs(w_ce)

eld = c/(3*vthe)

surfs = range(0, len(rhoHFS))

bounds=np.zeros([2,len(surfs)])
i=0

modeConversion = np.zeros(len(surfs))
nparas = np.linspace(0,10, 200)


PPDboundaryNpar = []
PPDboundaryRho = []
for i in range(len(rhoHFS)):
    accessRoots = getNparRoots(i, g)
    modeConversion[i] = np.min(accessRoots)
    """
    PPDroots = getNparRoots(i, f)
    PPDroots.sort()
    if len(PPDroots) > 0:
        PPDboundaryRho.append([rhoHFS[i] for root in PPDroots]) 
        PPDboundaryNpar.append(PPDroots)
    """
#PPDboundaryNpar = np.array(PPDboundaryNpar)
#PPDboundaryRho = np.array(PPDboundaryRho)

fig,ax = plt.subplots(figsize = (7.5,6))
ax.set_title(f'Discharge {shotNum}, $n_{{||,f}}$ = {n_para_f: 0.2f}\nenescal = {enescal}, Tescal = {tescal}, Tiscal = {tiscal}', y = 1.05)    
ax.plot(rhoHFS,eld,color='r',linestyle='--',label='Electron Landau Damping', linewidth = 2)  
ax.plot(rhoHFS,np.abs(modeConversion),color='b',linestyle='dotted',label='Mode Conversion', linewidth = 2)
#print(f'{PPDboundaryRho}')
#print(f'{PPDboundaryNpar}')

#ax.plot(PPDboundaryRho[:,0],np.abs(PPDboundaryNpar[:,0]),color='k', linewidth = 2)
#ax.plot(PPDboundaryRho[:,1],np.abs(PPDboundaryNpar[:,1]),color='k', linewidth = 2, label = 'Potential Propagation Domain')
#ax.plot([PPDboundaryRho[:,1][-1], PPDboundaryRho[:,1][-1]],np.abs([PPDboundaryNpar[:,0][-1], PPDboundaryNpar[:,1][-1]]),color='k', linewidth = 2)

cmap = matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(0,1),
         cmap = plt.get_cmap('jet'))
cmap.set_array([])
ticks = np.linspace(0,1,5)

cbar = fig.colorbar(cmap, ax = ax, shrink = .9, ticks = ticks, pad = .01)
cbar.set_label(r"Fractional power in ray")

ax.legend(loc='best',fontsize=14)
ax.set_ylim(0,6)
ax.set_xlabel(r'$\rho_{pol}$')
ax.set_ylabel('|-$N_\parallel|$')

addRays(ax)

fig.tight_layout()
plt.show()
