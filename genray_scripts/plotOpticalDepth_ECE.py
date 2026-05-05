###
# Plots the optical depth vs frequency
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
print(targetDir)
genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')#netCDF4.Dataset(f'{targetDir}/genray.nc','r')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 14)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

#adds the ray traces to ax
def main():
    wtau_em_nc = genray_ece_nc.variables['wtau_em_nc'][:,0]
    wfreq_nc = genray_ece_nc.variables['wfreq_nc'][:]
    print(wfreq_nc)

    fig,ax = plt.subplots()
    ax.plot(wfreq_nc, wtau_em_nc, lw = 3)
    ax.scatter(wfreq_nc, wtau_em_nc)
    ax.set_ylabel('Optical depth')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylim(bottom = 0)
    fig.tight_layout()

    plt.show()
    
main()
