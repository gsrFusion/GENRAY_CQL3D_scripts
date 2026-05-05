###
# Outputs the DT fusion power assuming n_T = n_D
###
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__);dname = os.path.dirname(abspath);os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__));parentdir = os.path.dirname(currentdir);sys.path.append(parentdir)

import getTargetInfo
import getInputFileDictionary
import helperFunctions as helper

#####Setup and do plotting#####
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)
plt.rc('legend',fontsize=16)

targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
inputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d')
cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

dVol = cql_nc.variables["dvol"][:] *1e-6#convert from cm^3 to m^3
rya = np.ma.getdata(cql_nc.variables["rya"][:])#convert from cm^3 to m^3

rho_pol, n_D = helper.getCQLnD(rho_pol = rya)
rho_pol, T_D = helper.getCQLTD(rho_pol = rya)

#Bosh-Hale https://library.psfc.mit.edu/catalog/online_pubs/MFE_formulary_2014.pdf page 80
#T in keV
def Bosh_Hale_reactivity(T):
    B_G = 34.3827
    m_muc2 = 1124656
    C1 = 1.173023e-9
    C2 = 1.51361e-2
    C3 = 7.51886e-2
    C4 = 4.60643e-3
    C5 = 1.35e-2
    C6 = -1.0675e-4
    C7 = 1.366e-5

    theta = T/(1 - (T*(C2+T*(C4+T*C6)))/(1+T*(C3+T*(C5+T*C7))))

    xi = (B_G**2/(4*theta))**(1/3)
    react = C1*theta*np.sqrt(xi/(m_muc2*T**3))*np.exp(-3*xi)

    return react*1e-6#convert to m^3/s

reactivities = Bosh_Hale_reactivity(T_D)
reactionRate = reactivities*n_D**2 #assumes n_D = n_T
powerDensity = (2.819831e-12)*reactionRate
P_fusion = np.sum(powerDensity*dVol)
print(f'P_fusion: {P_fusion/1e6} MW')

fig,ax = plt.subplots()
ax.plot(rho_pol, powerDensity/1e6)
ax.grid()
ax.set_ylabel('DT fusion power density (MW/m^3)')
ax.set_xlabel(r'$\rho_{pol}$')
fig.tight_layout()
plt.show()