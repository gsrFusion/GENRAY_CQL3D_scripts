import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import numpy as np
import matplotlib.pyplot as plt
from omfit_classes import omfit_eqdsk, utils_fusion
from scipy.interpolate import interp1d
import netCDF4
import getTargetInfo
import shotToEqdsk
topmostShotDir = getTargetInfo.getTopmostShotDir()
targetDir = getTargetInfo.getTargetDir()
machine = getTargetInfo.getMachine()        
shot = getTargetInfo.getShotNum()
eqdskName = shotToEqdsk.getEqdskName(shot, machine = machine)

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.rc('axes', labelsize = 20)
plt.rc('figure', titlesize = 18)
OMFIT_nc_derived = netCDF4.Dataset(f'{topmostShotDir}/FIT.nc','r')

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
rya = np.ma.getdata(cql_nc.variables["rya"][:])

gfile = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')
omfit_rho_pol = np.sqrt(gfile['fluxSurfaces']['levels'])
avgR = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['R'])(rya)
avga = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['a'])(rya)

q_prof = np.abs(interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['q'])(rya))


rho_psi = np.sqrt(OMFIT_nc_derived.variables["psi_n"])
if len(rho_psi) ==0 or len(rho_psi) == 1:
    rho_psi = rho_psi[0]
n_prof = interp1d(rho_psi, np.copy(OMFIT_nc_derived.variables["n_e"][0]))(rya)
Te_prof = interp1d(rho_psi, np.copy(OMFIT_nc_derived.variables["T_e"][0]))(rya)
Zeff_prof = interp1d(rho_psi, np.copy(OMFIT_nc_derived.variables['Zeff'][0]))(rya)

eta = utils_fusion.eta_0(zeff = Zeff_prof, te = Te_prof, ne = n_prof)
spitzerFactor = utils_fusion.Spitzer_factor(zeff = Zeff_prof)

neoEta_prof = utils_fusion.eta_par_neo(zeff = Zeff_prof, te = Te_prof, 
            ne = n_prof, q = q_prof, R0 = avgR, r_minor = avga)

e = 1.6e-19
m_e = 9.109e-31
eps_0 = 8.85e-12
eVToJoule = 1.6e-19

lnLambda = 16
eta1 =  np.pi*e**2*np.sqrt(m_e)*lnLambda*Zeff_prof/(16*np.pi**2*eps_0**2*(eVToJoule*Te_prof)**1.5)
eta1 *= (1+1.198*Zeff_prof + 0.222*Zeff_prof**2)/(1+2.966*Zeff_prof + 0.753*Zeff_prof**2)

surfArea = gfile['fluxSurfaces']['geo']['surfArea']
dArea = np.insert(np.diff(surfArea),0,surfArea[0])
avgJt =  gfile['fluxSurfaces']['avg']['Jt']

darea_cql = np.ma.getdata(cql_nc.variables["darea"][:])/1e4#convert to m^2
area_cql = np.ma.getdata(cql_nc.variables["area"][:])/1e4#convert to m^2
curreq = np.ma.getdata(cql_nc.variables["curreq"][:])*1e4#convert to A/m^2
currentCumsum = np.cumsum(curreq**2*darea_cql)

dvol = np.ma.getdata(cql_nc.variables["darea"][:])/1e6#convert to m^3
Btot2 =  interp1d(omfit_rho_pol,gfile['fluxSurfaces']['avg']['Btot**2'])(rya)
Bcumsum = np.cumsum(dvol * Btot2)

mu = 4*np.pi*1e-7
L_i = (Bcumsum/mu)/(currentCumsum)

omfit_rho_pol = np.sqrt(gfile['fluxSurfaces']['levels'])
resistance = eta1 * 2* np.pi *avgR
resistance_neo = neoEta_prof * 2* np.pi *avgR
fig,ax = plt.subplots()
ax.plot(rya,L_i/resistance_neo, lw = 2)

ax2 = ax.twinx()
#ax2.plot(rya, resistance, lw = 2, color = 'r')
ax2.plot(rya, resistance_neo, lw = 2, color = 'g')
ax2.plot(rya, L_i, lw = 2, color = 'k')

plt.show()