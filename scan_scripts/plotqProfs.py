import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)

import numpy as np
import matplotlib.pyplot as plt
import netCDF4

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

machine = 'NTPT'

fig,ax = plt.subplots()


if machine == 'NTPT':
    fakeDevice = 'DIIID'
    
    stem = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeDevice}'

    if fakeDevice == 'DIIID':

        fakeShots = ['147634PT','147634NT','193765PT','193765NT']
        nparas = ['n2.5', 'n2.5', 'p2.5','p2.5']
        colors = ['mediumblue', 'crimson' ,'darkturquoise', 'orange']
        labels = ['PT 147634-like', 'NT 147634-like', 'PT 193765-like', 'NT 193765-like']
        power = 1
    elif fakeDevice == 'ARC':
        fakeShots = ['V3APT','V3ANT']
        nparas = ['n2.0', 'n2.0']
        colors = ['teal', 'goldenrod']
        labels = ['PT V3A-like', 'NT V3A-like']
        power = 10

    for i, fakeShot in enumerate(fakeShots):
        targetDir = f'{stem}.{fakeShot}/{machine}_{fakeDevice}.{fakeShot}_{nparas[i]}Npara_-0.5grillHeight_{power}MW'

        cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
        rya = cql_nc.variables["rya"][:]
        q_prof = cql_nc.variables["qsafety"][:]

        ax.plot(rya, np.abs(q_prof), lw = 3, color = colors[i], label = labels[i])


ax.set_xlabel(r'$\rho_{pol}$')
ax.set_ylabel(f'safety factor')
ax.grid()
ax.set_ylim([0,6])
ax.legend()
fig.tight_layout()
plt.savefig(f'toka_DIIID_q.jpg',dpi=300)

plt.show()