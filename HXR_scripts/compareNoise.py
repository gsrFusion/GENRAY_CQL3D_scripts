import numpy as np
import matplotlib.pyplot as plt

import HXR_tomo

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.rc('axes', labelsize = 20)
plt.rc('figure', titlesize = 18)

import netCDF4
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

parentParentDir = os.path.dirname(parentdir)
sys.path.append(parentParentDir)
cql_nc = netCDF4.Dataset(f'{parentParentDir}/cql3d.nc','r')

rya = np.ma.getdata(cql_nc.variables["rya"][:])

_imageBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])
args0 = np.array([_imageBand, False, False, 50, 250, .62, True, False, .8, 0])
args1 = np.copy(args0); args1[-1] = .1
args2 = np.copy(args0); args2[-1] = .3
args3 = np.copy(args0); args3[-1] = .5

tomo0 = HXR_tomo.main([args0])
tomo1 = HXR_tomo.main([args1])
tomo2 = HXR_tomo.main([args2])
tomo3 = HXR_tomo.main([args3])

fig, ax = plt.subplots()
ax.plot(tomo0.grid,tomo0.solution/np.max(tomo0.solution), lw = 2,  label = r"No error")
ax.plot(tomo1.grid,tomo1.solution/np.max(tomo1.solution), lw = 2, label = r"0-10% error")
ax.plot(tomo3.grid,tomo2.solution/np.max(tomo2.solution), lw = 2, label = r"20-30% error")
ax.plot(tomo3.grid,tomo3.solution/np.max(tomo3.solution), lw = 2, label = r"40-50% error")

ax.plot(rya, tomo0.ne/np.max(tomo0.ne), label = r"$n_e(E \geq $" + f"${tomo0.E_pMin}$ keV$)$", color = 'k',
            linestyle = 'dashed', linewidth = 2)
ax.axhline(0,c='k', linestyle = 'dotted')
ax.set_ylim([-.05,1.1])
ax.set_xlabel(r"$\rho_{{pol}}$", fontsize = 26)
ax.set_ylabel("Normalized Units", fontsize = 18)
ax.set_xticks([-1,-.75,-.5,-.25,0,.25,.5,.75,1])     
ax.ticklabel_format(axis = 'y', scilimits = (0,0))
legend = ax.legend(fontsize = 17, loc = 'upper left')#bbox_to_anchor=(.4,1), loc='upper center') 

legend.get_frame().set_alpha(None)
legend.get_frame().set_facecolor((1, 1, 1, 0.1))  
     
fig.set_size_inches((7,6.5))#8,4.75
fig.tight_layout()
ax.set_xlim([0,1.01])

plt.show()
