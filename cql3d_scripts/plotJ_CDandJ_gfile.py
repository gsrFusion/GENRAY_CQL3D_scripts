###
# Plots the total current density profile, the driven current density, and the bootstrap current density
###

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
import getTargetInfo
import getInputFileDictionary

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 18)
plt.rc('legend', fontsize = 14)

cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d')

targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

eqdskName = shotToEqdsk.getEqdskName(shotNum, machine = machine)
geqdsk = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
rya = np.ma.getdata(cql_nc.variables["rya"][:])

print(f'peak of JLh at rho_pol = {rya[np.argmax(curr)]}')

curr_gfile = geqdsk.surfAvg('Jpar', interp = 'cubic')/1e6#geqdsk["fluxSurfaces"]["avg"]["Jt"]/1e6#MA/m^2
psi_N_gfile = geqdsk["fluxSurfaces"]["levels"]
rya = np.ma.getdata(cql_nc.variables["rya"][:])
darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
totalCD = np.sum(curr*darea)
print(f'totalCD: {totalCD}')

#Needed for the bootstrap current calculation
rho_pol, n_e = helper.getCQLne(targetDir = targetDir,rho_pol=rya)
rho_pol, T_e = helper.getCQLTe(targetDir = targetDir,rho_pol=rya)
rho_pol, Zeff = helper.getCQLZeff(targetDir = targetDir,rho_pol=rya)
rho_pol, n_D = helper.getCQLnD(targetDir = targetDir,rho_pol=rya)
rho_pol, n_C = helper.getCQLnC(targetDir = targetDir,rho_pol=rya)
rho_pol, T_D = helper.getCQLTD(targetDir = targetDir,rho_pol=rya)

pressure = (n_e*T_e + n_D*T_D + n_C*T_D)*1.602e-16
J_BS_prof = utils_fusion.sauter_bootstrap(gEQDSKs = geqdsk, psi_N = rya**2, 
                Ti = np.array([T_D*1e3]), ne = np.array([n_e]), Te = np.array([T_e*1e3]),
                charge_number_to_use_in_ion_collisionality = 'Koh', charge_number_to_use_in_ion_lnLambda = 'Koh',
                Zis=[1,6], nis = np.array([[n_D], [n_C]]), R0 = 1.6955, p = np.array([pressure]), version = 'osborne')[0]

fig,ax = plt.subplots()
ax.plot(rya, curr, lw = 3, color ='k', label = r'$J_{LH}$')
ax.plot(psi_N_gfile**.5, np.abs(curr_gfile), lw = 2, label = r'$J_{||}$')
ax.plot(rya, np.abs(J_BS_prof/1e6), lw = 2, label = r'$J_{BS}$')
ax.set_xlabel(r'$\rho_{pol}$')
ax.set_ylabel(r'LH-driven current density (MA/$m^2$)')
ax.set_title(f'Shot {shotNum.split(".")[0]}, {shotNum.split(".")[1][1:]} ms', loc = 'right')

fig.tight_layout()
ax.legend(loc = 'best')

plt.show()
