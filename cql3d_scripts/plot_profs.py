###
# Plots the electron and ion densities and temperatures according to the cql3d input file
###

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__);dname = os.path.dirname(abspath);os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__));parentdir = os.path.dirname(currentdir);sys.path.append(parentdir)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

import getInputFileDictionary
inputFileDict = getInputFileDictionary.getInputFileDictionary('cql3d')

import numpy as np
import matplotlib.pyplot as plt
import helperFunctions as helper

style = 'onePlot'

if style == 'onePlot':
    nrows = 1
else:
    nrows = 2

rho_pol, n_e = helper.getCQLne()
rho_pol, T_e = helper.getCQLTe()
rho_pol, n_D = helper.getCQLnD()
rho_pol, T_D = helper.getCQLTD()
#rya, q = helper.getCQLq()

dischargeNumber = os.getcwd().split('/')[-1]
if len(dischargeNumber) != 6:
    dischargeNumber = os.getcwd().split('/')[-2]


#####Setup and do plotting#####
plt.rc('xtick', labelsize = 15)
plt.rc('ytick', labelsize = 15)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)
plt.rc('legend',fontsize=16)

if style == 'onePlot':
    figsize = (6.4,4.5)#(6.4,4.8)
else:
    figsize = (6.4,5.5)

fig,axes = plt.subplots(nrows = nrows, figsize=figsize)

if style == 'onePlot':
    axes.set_xlim([0,1])
    axes2 = axes.twinx()

    axes.plot(rho_pol, n_e/1e19, label = r'$n_e$', lw = 3, color = 'green')
    axes2.plot(rho_pol, T_e, label = r'$T_e$', lw = 3,linestyle = 'dashed', color = 'green')
    axes.plot([.2,.4],[-1,-2], label = r'$T_e$', lw = 3,linestyle = 'dashed', color = 'green')
    #q2Index = helper.findNearestIndex(2, np.abs(q))
    #axes.axvline(rya[q2Index],label = r'$q=2$', lw = 2,linestyle = 'dashdot', color = 'k')

    axes.plot(rho_pol, n_D/1e19, label = r'$n_D$', lw = 3, color = 'royalblue')
    axes2.plot(rho_pol, T_D, label = r'$T_D$', lw = 3,linestyle = 'dashed', color = 'royalblue')
    axes.plot([.2,.4],[-1,-2], label = r'$T_D$', lw = 3,linestyle = 'dashed', color = 'royalblue')


    
    #axes.plot(rya, q, label = r'$q$', lw = 3,linestyle = 'dashdot', color = 'k')
    """
    if shotNum == '203619.04130':
        import netCDF4
        sliceNC = netCDF4.Dataset(f'/home/grantr/codes/GENRAY_CQL3D_scripts/genray_batch/DIIID_shots/DIIID_203619.04130/TS_203619.04129.nc','r')
        
        print(f'uncertainty shape: {sliceNC.variables["n_e__uncertainty"][:,0].shape}')

        axes.errorbar(np.sqrt(sliceNC.variables['psi_n'][:]), sliceNC.variables['n_e'][:]/1e19, 
                      yerr = sliceNC.variables['n_e__uncertainty'][:,0]/1e19, marker = 'D', 
                      markersize = 5,zorder = 10, linestyle='none',
                      color='k', label = 'DIII-D Measurements')
        
        axes.errorbar(np.sqrt(sliceNC.variables['psi_n'][:]), sliceNC.variables['T_e'][:]/1e3, 
                      yerr = sliceNC.variables['T_e__uncertainty'][:,0]/1e3, marker = 'D', 
                      markersize = 5,zorder = 10, linestyle='none',
                      color='k', label = 'DIII-D Measurements')

        print(sliceNC.variables.keys())
    """



    axes.set_ylabel(r'Density ($10^{19}m^{-3}$)')
    axes2.set_ylabel(r'Temperature (keV)')
    axes.set_xlabel(r'$\rho_{pol}$')
    axes.grid()
    axes.legend(ncol = 2)
    #axes.legend(loc='lower center', bbox_to_anchor=(0.5,.99),ncol=3,labelspacing=0.3)
    axes.set_ylim([0,4])
    axes2.set_ylim([0,4])
    

else:

    axes[0].plot(rho_pol, n_e/1e19, lw = 3, color = 'tab:green')
    #axes[0].plot(rho_pol, n_D/1e19, label = 'D', lw = 3, color = 'tab:blue')

    axes[1].plot(rho_pol, T_e, lw = 3, color = 'tab:green')
    #axes[1].plot(rho_pol, T_D, label = 'D', lw = 3, color = 'tab:blue')

    if shotNum == '203619.04130':
        import netCDF4
        sliceNC = netCDF4.Dataset(f'/home/grantr/codes/GENRAY_CQL3D_scripts/genray_batch/DIIID_shots/DIIID_203619.04130/TS_203619.04129.nc','r')
        
        print(f'uncertainty shape: {sliceNC.variables["n_e__uncertainty"][:,0].shape}')

        axes[0].errorbar(np.sqrt(sliceNC.variables['psi_n'][:]), sliceNC.variables['n_e'][:]/1e19, 
                      yerr = sliceNC.variables['n_e__uncertainty'][:,0]/1e19, marker = 'd', 
                      markersize = 5,zorder = 10, linestyle='none',
                      mfc='none', mec='k',ecolor='k',label = 'TS data')
        
        axes[1].errorbar(np.sqrt(sliceNC.variables['psi_n'][:]), sliceNC.variables['T_e'][:]/1e3, 
                      yerr = sliceNC.variables['T_e__uncertainty'][:,0]/1e3, marker = 'd', 
                      markersize = 5,zorder = 10, linestyle='none',
                      mfc='none', mec='k',ecolor='k',label = 'TS data')

    axes[0].set_ylabel(r'n$_e$ ($10^{19}m^{-3}$)')
    axes[0].set_xlabel(r'$\rho_{pol}$')
    axes[0].grid()
    axes[0].set_xlim(right=1.025)
    #axes[0].set_yticks([0,.5,1,1.5,2,2.5])
    #axes[0].legend()

    axes[1].set_ylabel(r'T$_e$ (keV)')
    axes[1].set_xlabel(r'$\rho_{pol}$')
    axes[1].grid()
    #axes[1].set_yticks([0,.5,1,1.5,2,2.5])
    axes[1].set_xlim(right=1.05)
    #axes[1].legend()


#fig.suptitle(f'CQL3D {shotNum} Profiles\n enescal = {enescal}, Tescal = {tescal}, Tiscal = {tiscal}')

    axes[0].set_ylim([0,2.5])
    axes[1].set_ylim([0,2.5])

fig.tight_layout()
plt.savefig('206629_profs.jpeg',dpi=300)
plt.show()