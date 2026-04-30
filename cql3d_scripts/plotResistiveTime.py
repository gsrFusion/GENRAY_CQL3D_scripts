import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os, sys
from scipy.interpolate import interp1d
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import warnings
warnings.filterwarnings("ignore")
import getInputFileDictionary
import helperFunctions as helper
from omfit_classes import omfit_eqdsk
from omfit_classes import utils_fusion
import shotToEqdsk

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('figure', titlesize = 22)
plt.rc('legend', fontsize = 14)

import getInputFileDictionary
cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d')

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

eqdskName = shotToEqdsk.getEqdskName(shotNum, machine = machine)
geqdsk = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')

avgR = geqdsk['fluxSurfaces']['avg']['R']
Vloop = .5 #0.5 V for 203912 at 2700 ms
print(f'{Vloop/(2*np.pi*avgR)}')

rho_pol, n_e = helper.getCQLne(targetDir = targetDir)
rho_pol, T_e = helper.getCQLTe(targetDir = targetDir)
rho_pol, Zeff = helper.getCQLZeff(targetDir = targetDir)
rho_pol, n_D = helper.getCQLnD(targetDir = targetDir)
rho_pol, n_C = helper.getCQLnC(targetDir = targetDir)
rho_pol, T_D = helper.getCQLTD(targetDir = targetDir)

sigma_neo = utils_fusion.nclass_conductivity_from_gfile(
    psi_N=rho_pol**2,
    Te=np.array([T_e*1e3]),
    ne=np.array([n_e]),
    Ti=np.array([T_D*1e3]),
    gEQDSK=geqdsk,
    Zeff=np.array([Zeff]),
    nis = np.array([[n_D], [n_C]]),
    Zis=[1,6],
    Zdom=1,
    return_info_pack=False,
    plot_slice=None,  # Set to a time slice index to plot, set to None for no plot
    charge_number_to_use_in_ion_collisionality='Koh',
    charge_number_to_use_in_ion_lnLambda='Koh',
)[0]

sigma_neo = interp1d(geqdsk['fluxSurfaces']['levels']**.5, sigma_neo)(rho_pol)

eta_neo = 1/sigma_neo

e=1.602e-19
m_e = 9.109e-31
coLo = 16
eps_0 = 8.854e-12
eta_spi = (4*np.sqrt(2*np.pi)/3)*(Zeff*e**2 *np.sqrt(m_e)*coLo)/((4*np.pi*eps_0)**2 * (T_e*1.602e-16)**(3/2))

eta_spi *= (1+1.198*Zeff+0.222*Zeff**2)/(1+2.966*Zeff+.753*Zeff**2)

print(f'eta_spi[0], eta_neo[0]: {eta_spi[0], eta_neo[0]}')

mu_0 = 4*np.pi*1e-7
rs = geqdsk['fluxSurfaces']['avg']['a']
r_func = interp1d(geqdsk['fluxSurfaces']['levels'], rs)

L = rs[-1]
resistiveTime_neo = mu_0*L**2/(3.832**2*eta_neo)
resistiveTime_spi = mu_0*L**2/(3.832**2*eta_spi)
fig,ax = plt.subplots()
ax.plot(rho_pol, resistiveTime_neo, label = r'with $\eta_{neo}$', lw = 3)
ax.plot(rho_pol, resistiveTime_spi, label = r'with $\eta_{spitzer}$', lw = 3)
ax.set_ylabel(r'Resistive diffusion time $\mu_0 a^2/\eta$')
ax.set_xlabel(r'radial coordinate ($\rho_{pol}$)')
fig.tight_layout()
ax.legend()
plt.show()
