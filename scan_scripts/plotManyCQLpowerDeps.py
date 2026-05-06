###
# Plots the power deposition predicted by CQL3D for several simulations
###

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


import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
machine = getTargetInfo.getMachine()

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 14)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 12)

machine = 'DIIID'

if machine == 'DIIID':

    shotNum = '203912.02700'
    
    case = '203912 2700 E scan'

    if case == '203912 2700 E scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum'

        includeWallRef = True

        targetDirs = [
            f'{stem1}_first',
            f'{stem1}_p0.025Fld_first',
            f'{stem1}_p0.05Fld_first',
        ]
        labels = [
            r'E=0 V/m',
            r'E=0.025 V/m',
            r'E=0.05 V/m',
        ]
        colors = [ 'tab:blue', 'tab:green','tab:red', 'tab:purple', 'tab:grey']

fig,ax = plt.subplots()

#ax.set_title('prmt4 = 0.001, rksteps = 35000')

#adds the ray traces to ax
def plotPowerDep(targetDir, label = ''):
    print(targetDir)
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    powrft = cql_nc.variables["powrft"][:][-1,:]
    rya = cql_nc.variables["rya"][:]

    ax.plot(rya, powrft, lw = 3, label = label)


for i, targetDir in enumerate(targetDirs):
    plotPowerDep(targetDir,label = labels[i])

ax.set_ylabel("power (W/cm^3)")
ax.set_xlabel("rho_pol")
ax.legend()
fig.tight_layout()
plt.show()


    

