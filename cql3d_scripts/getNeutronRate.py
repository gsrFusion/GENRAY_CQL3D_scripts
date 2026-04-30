import numpy as np
import matplotlib.pyplot as plt
import os,sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from scipy.interpolate import interp1d

import helperFunctions as helper
from omfit_classes import omfit_eqdsk

import shotToEqdsk
import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

eqdskName = shotToEqdsk.getEqdskName(shotNum, machine = machine)
geqdsk = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')

rho_psi, Ti = helper.getCQLTD()

rho_psi, n_D = helper.getCQLnD()

sigmav_DD = 2.33e-14*Ti**(-2/3)*np.exp(-18.76*Ti**(-1/3)) #cm^3/s

numReactions = n_D**2/2*sigmav_DD/1e6

rho_pols = np.sqrt(geqdsk['fluxSurfaces']['levels'])
vol = geqdsk['fluxSurfaces']['geo']['vol']
dvol = vol[1:]-vol[:-1]
rho_centers = (rho_pols[1:] + rho_pols[:-1])/2

surfArea = geqdsk['fluxSurfaces']['geo']['surfArea']

dvol_interp1d = interp1d(rho_centers, dvol, fill_value = "extrapolate", bounds_error = False)(rho_psi)
surfaceArea = interp1d(rho_pols, surfArea, fill_value = "extrapolate", bounds_error = False)(.99)

totalNeutrons = np.sum(dvol_interp1d * numReactions)
print(f'totalNeutrons: {totalNeutrons:.2e}')
print(f'approximate neutron flux: {totalNeutrons/surfaceArea:.2e} (neutrons/m^2)')