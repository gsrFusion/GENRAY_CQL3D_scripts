import numpy as np
import matplotlib.pyplot as plt
import os, sys
import netCDF4

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import helperFunctions as helper

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 13)

#returns the distribution function integrated over pitch angle
def getIntegratedOverPitchF(targetDir):
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

    rya = cql_nc.variables["rya"][:]#radial points at which f was solved

    #distribution function (lrz, jx, iy)
    #vnorm^3*s^3/cm^6
    f = cql_nc.variables["f"][:]

    #pitch angles mesh at which f is defined
    pitchAngleMesh = np.ma.getdata(cql_nc.variables["y"][:])

    #energy  = restmkev*(gamma-1)
    #energies corresponding to velocities jx
    enerkev = cql_nc.variables["enerkev"][:] 

    #minimum energy of particles to consider for plotting
    minEnergy = 50
    #index of that minimum energy in enerkev
    #this index is also the index for the corresponding velocity
    minEnergyIndex = 0#np.where(enerkev < minEnergy)[0][-1]
    #distribution function for energetic particles
    energeticF = f[:,:,:]


    energeticF_integOverPitch = np.zeros((len(rya), len(enerkev[minEnergyIndex:])))
    for rhoIndex in range(len(rya)):
        #this is the angular part of the spherical jacobian
        integOverPitch = 2*np.pi*np.trapz(energeticF[rhoIndex,:]*np.sin(pitchAngleMesh[rhoIndex]), pitchAngleMesh[rhoIndex], axis = 1)
        energeticF_integOverPitch[rhoIndex,:] = integOverPitch

    return energeticF_integOverPitch

def compareDistributionSlices():
    targetDirs = ['/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_203619.04120/DIIID_203619.04120_expSpectrum_2Zeff_second',
                  '/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_203619.04135/DIIID_203619.04135_expSpectrum_2Zeff',]

    rhosToPlot = [.1,.2,.3,.4,.5,.6,.7]

    fig, ax = plt.subplots(figsize=(8,5.5))

    titleString = ''

    integratedFs = [getIntegratedOverPitchF(targetDirs[0]), getIntegratedOverPitchF(targetDirs[1])]

    cql_nc = netCDF4.Dataset(f'{targetDirs[0]}/cql3d.nc','r')
    rya = cql_nc.variables["rya"][:]#radial points at which f was solved
    enerkev = cql_nc.variables["enerkev"][:] 

    for i in range(len(rhosToPlot)):
        rho_pol = rhosToPlot[i]
        rho_index = helper.findNearestIndex(rho_pol, rya)

        ax.plot(enerkev, integratedFs[0][rho_index,:] - integratedFs[1][rho_index,:], lw = 2, 
                label = r'$\rho_{pol} = $' + f'{rho_pol}')

    ax.legend(loc='best',ncol=2)
    ax.set_ylabel(r'Dist func (vnorm$^3$*sec$^3$/cm$^6$)')
    ax.set_title('4120 - 4135')
    ax.set_xlabel('Energy (keV)')
    ax.set_xlim([50,250])
    ax.set_ylim([-1e13,.25e13])
    fig.tight_layout()

    plt.show()


compareDistributionSlices()