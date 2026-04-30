"""
Plots the ray traces and the RF power deposition density
"""
import numpy as np
import matplotlib.pyplot as plt

import os, sys
from scipy.signal import find_peaks
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import netCDF4

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 14)
plt.rc('legend', fontsize = 12)

#plots either the toroidal and/or poloidal ray trajectories
def addECE(targetDir, ax, label = '', color = None, style = 'solid', includeWallRef = True, minFreq = 0):
    genray_ece_nc = netCDF4.Dataset(f'{targetDir}/genray_ece.nc','r')

    wtemp_rad_fr_nc = genray_ece_nc.variables['wtemp_rad_fr_nc'][:,0]
    wtemp_rad_fr_wall_nc = genray_ece_nc.variables['wtemp_rad_fr_wall_nc'][:,0]

    if includeWallRef:
        T_rad = wtemp_rad_fr_wall_nc
    else:
        T_rad = wtemp_rad_fr_nc

    freqs = genray_ece_nc.variables['wfreq_nc'][:]

    zorder = 1
    if color == 'tab:blue':
        zorder = 5
    
    ax.plot(freqs[freqs > minFreq], T_rad[freqs > minFreq], lw = 3, 
            label = label, zorder = zorder,color = color, linestyle = style)
    #ax.axvline(freqs[np.argmax(T_rad)], lw =2, color = 'k')
    #ax.scatter(freqs, T_rad,  color = color,)

def main():

    fig,ax = plt.subplots()


    shotNum = '203619.04130'
    stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum'
    
    case = '203619 4130 prmt6 ECE scan'

    if case == '203912 2700 on off':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_noLH'
        includeWallRef = True

        targetDirs = [
            [f'{stem1}_first', f'{stem1}_second', f'{stem1}_third', f'{stem1}_fourth', ],
            [f'{stem2}_first', f'{stem2}_second', f'{stem2}_third', f'{stem2}_fourth', ],
        ]
        labels = [
            r'With LH',
            r'Without LH',
        ]
        colors = [ 'tab:blue', 'tab:red','tab:red', 'tab:purple', 'tab:grey']
        ax.set_xlim([75,105])

    if case == '203912 2700 E scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum'

        includeWallRef = True

        targetDirs = [
            [f'{stem1}_first',f'{stem1}_second',f'{stem1}_third',f'{stem1}_fourth'],
            [f'{stem1}_p0.025Fld_first',f'{stem1}_p0.025Fld_second',f'{stem1}_p0.025Fld_third',f'{stem1}_p0.025Fld_fourth'],
            #[f'{stem1}_p0.025Fld_0.75pwrFactor_second'],
            [f'{stem1}_p0.05Fld_first',f'{stem1}_p0.05Fld_second',f'{stem1}_p0.05Fld_third',f'{stem1}_p0.05Fld_fourth'],
        ]
        labels = [
            r'E=0 V/m',
            r'E=0.025 V/m',
            #r'E=0.025 V/m, 0.75 pwrFactor',
            r'E=0.05 V/m',
        ]
        colors = [ 'tab:blue', 'tab:green','tab:red', 'tab:purple', 'tab:grey']
        ax.set_xlim([75,105])


    if case == '203619 4130 on off':
        ax.set_xlim([80,115])#115
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff_noLH/DIIID_{shotNum}_expSpectrum_2Zeff_noLH'
        includeWallRef = True

        targetDirs = [
            [f'{stem1}_first', f'{stem1}_second', f'{stem1}_third', f'{stem1}_fourth', ],
            #[f'{stem2}_first', f'{stem2}_second', f'{stem2}_third', f'{stem2}_fourth', ],
        ]
        labels = [
            r'With LH',
            r'Without LH',
        ]
        colors = [ 'tab:blue', 'tab:red','tab:red', 'tab:purple', 'tab:grey']

    if case == '203619 4130 port scan':
        ax.set_xlim([80,105])
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        includeWallRef = True

        targetDirs = [
            #[f'{stem1}_first', f'{stem1}_second'],
            [f'{stem1}_first_noPort', f'{stem1}_second_noPort'],
        ]
        labels = [
            r'With port',
            r'Without port',
        ]
        colors = [ 'tab:red', 'tab:red','tab:red', 'tab:purple', 'tab:grey']

    if case == '203619 4130 E resolution scan':
        stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}'
        includeWallRef = True

        targetDirs = [
            [f'{stem}/DIIID_{shotNum}_expSpectrum_2Zeff_p0.05Fld'],
            [f'{stem}/DIIID_{shotNum}_expSpectrum_2Zeff_p0.05Fld_1200jx'],
            [f'{stem}/DIIID_{shotNum}_expSpectrum_2Zeff_p0.05Fld_1200jx_6000enorm'],
        ]
        labels = [
            r'jx = 1000, enorm = 5000',
            r'jx = 1200, enorm = 5000',
            r'jx = 1200, enorm = 6000',
        ]
        colors = [ 'tab:blue', 'tab:red','tab:red', 'tab:purple', 'tab:grey']
        ax.set_xlim([80,115])

    if case == '203619 4130 id scan':
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff_id2/DIIID_{shotNum}_expSpectrum_2Zeff_id2'
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        includeWallRef = True

        targetDirs = [
            [f'{stem1}_first', f'{stem1}_second', f'{stem1}_third', f'{stem1}_fourth', ],
            [f'{stem2}_first', f'{stem2}_second', f'{stem2}_third', f'{stem2}_fourth', ],
        ]
        labels = [
            r'Thermal effects included',
            r'Cold plasma',
        ]
        colors = [ 'tab:blue', 'tab:red','tab:red', 'tab:purple', 'tab:grey']
        ax.set_xlim([80,115])

    if case == '203619 4130 Zeff scan':
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_1.6Zeff/DIIID_{shotNum}_expSpectrum_1.6Zeff'
        stem3 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2.4Zeff/DIIID_{shotNum}_expSpectrum_2.4Zeff'
        includeWallRef = True

        targetDirs = [
            #[f'{stem1}_first', f'{stem1}_second', f'{stem1}_third', f'{stem1}_fourth', ],
            [f'{stem2}_first', f'{stem2}_second', f'{stem2}_third', f'{stem2}_fourth', ],
            #[f'{stem3}_first', f'{stem3}_second', f'{stem3}_third', f'{stem3}_fourth', ],
        ]
        labels = [
            #r'Zeff = 1.6',
            r'With LH, Zeff = 2',
            #r'Zeff = 2.4',
        ]
        colors = [ 'tab:red', 'tab:blue','seagreen', 'tab:purple', 'tab:grey']
        colors = [ 'tab:blue', 'tab:blue','seagreen', 'tab:purple', 'tab:grey']
        ax.set_xlim([80,115])

    if case == '203619 4130 nfreq scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_5nfreq'],
            [f'{stem2}_10nfreq'],
            [f'{stem1}_second'],

        ]
        labels = [
            r'nfreq = 5',
            r'nfreq = 10',
            r'nfreq = 15',
        ]
        colors = [ 'tab:green', 'tab:red','tab:blue', 'tab:grey', 'tab:green']

        ax.set_xlim([80,105])
        plt.rc('legend', fontsize = 12)

    if case == '203619 4130 ireflm scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_ireflm1_lowerAcc'],
            [f'{stem2}_ireflm2_lowerAcc'],
            [f'{stem2}_ireflm3_lowerAcc'],

        ]
        labels = [
            r'ireflm = 1',
            r'ireflm = 2',
            r'ireflm = 3',
        ]
        colors = [ 'tab:green', 'tab:blue','tab:red', 'tab:grey', 'tab:green']
        ax.set_xlim([80,105])

    if case == '203619 4130 prmt4 ECE scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_1e-3prmt4ECE'],
            [f'{stem2}_1e-4prmt4ECE'],
            [f'{stem2}_1e-5prmt4ECE'],
            [f'{stem2}_1e-6prmt4ECE'],
            [f'{stem1}_second'],

        ]
        labels = [
            r'prmt4ECE = 1e-3',
            r'prmt4ECE = 1e-4',
            r'prmt4ECE = 1e-5',
            r'prmt4ECE = 1e-6',
            r'prmt4ECE = 1e-7',
        ]
        colors = [ 'tab:blue', 'tab:red','tab:purple', 'tab:grey', 'tab:green']

        ax.set_xlim([80,105])
        plt.rc('legend', fontsize = 10)


    if case == '203619 4130 nnkpar scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_20nnkpar'],
            [f'{stem2}_30nnkpar'],
            [f'{stem2}_40nnkpar'],
            [f'{stem2}_50nnkpar'],
            [f'{stem1}_second'],

        ]
        labels = [
            r'nnkpar = 20',
            r'nnkpar = 30',
            r'nnkpar = 40',
            r'nnkpar = 50',
            r'nnkpar = 60',
        ]
        colors = [ 'tab:orange', 'tab:red','tab:purple', 'tab:grey','tab:green', 'tab:blue']
        plt.rc('legend', fontsize = 10)
        ax.set_xlim([80,105])

    if case == '203619 4130 enorm scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_5500enorm'],
            [f'{stem2}_3500enorm'],
            [f'{stem2}_2500enorm'],
            [f'{stem1}_second'],

        ]
        labels = [
            r'enorm = 5500',
            r'enorm = 3500',
            r'enorm = 2500',
            r'enorm = 1500',
        ]
        colors = [ 'tab:red', 'tab:green','tab:purple', 'tab:blue','tab:green', 'tab:blue']
        plt.rc('legend', fontsize = 10)
        ax.set_xlim([80,105])

    if case == '203619 4130 xfac scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_0.5xfac'],
            [f'{stem2}_0.25xfac'],
            [f'{stem1}_second'],

        ]
        labels = [
            r'xfac = 0.5',
            r'xfac = 0.25',
            r'xfac = 0.1',
        ]
        colors = [ 'tab:red', 'tab:green','tab:blue', 'tab:blue','tab:green', 'tab:blue']
        plt.rc('legend', fontsize = 11)
        ax.set_xlim([80,105])

    if case == '203619 4130 jx scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_375jx'],
            [f'{stem2}_500jx'],
            [f'{stem2}_625jx'],
            [f'{stem1}_second'],

        ]
        labels = [
            r'jx = 375',
            r'jx = 500',
            r'jx = 625',
            r'jx = 750',
        ]
        colors = [ 'tab:purple', 'tab:red','tab:green', 'tab:blue','tab:green', 'tab:blue']
        plt.rc('legend', fontsize = 11)
        ax.set_xlim([80,105])

    if case == '203619 4130 prmt4 LH scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_1e-4prmt4LH'],
            [f'{stem2}_1e-5prmt4LH'],
            [f'{stem2}_1e-6prmt4LH'],
            [f'{stem2}_1e-7prmt4LH'],
            [f'{stem1}_second'],

        ]
        labels = [
            r'prmt4LH = 1e-4',
            r'prmt4LH = 1e-5',
            r'prmt4LH = 1e-6',
            r'prmt4LH = 1e-7',
            r'prmt4LH = 1e-8',
        ]
        colors = [ 'tab:green', 'tab:grey','tab:red', 'tab:purple', 'tab:blue']
        plt.rc('legend', fontsize = 10)
        ax.set_xlim([80,105])

    if case == '203619 4130 prmt6 LH scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_1e-2prmt6LH'],
            [f'{stem1}_second'],
            [f'{stem2}_2.5e-3prmt6LH'],


        ]
        labels = [
            r'prmt6LH = 1e-2',
            r'prmt6LH = 5e-3',
            r'prmt6LH = 2.5e-3',
        ]
        colors = [ 'tab:green', 'tab:blue','tab:red', 'tab:purple', 'tab:blue']
        plt.rc('legend', fontsize = 10)
        ax.set_xlim([80,105])

    if case == '203619 4130 prmt6 ECE scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_5e-6prmt4ECE_5e-3prmt6ECE'],
            [f'{stem2}_5e-6prmt4ECE_2.5e-3prmt6ECE'],
            [f'{stem2}_5e-6prmt4ECE_2e-3prmt6ECE'],
            [f'{stem2}_5e-6prmt4ECE_1.5e-3prmt6ECE'],
            [f'{stem2}_5e-6prmt4ECE_1e-3prmt6ECE'],
            [f'{stem2}_5e-6prmt4ECE_9e-4prmt6ECE'],

        ]
        labels = [
            r'prmt6LH = 5e-3',
            r'prmt6LH = 2.5e-3',
            r'prmt6LH = 2e-3',
            r'prmt6LH = 1.5e-3',
            r'prmt6LH = 1e-3',
            r'prmt6LH = 9e-4',
        ]

        #"""
        targetDirs = [
            [f'{stem2}_5e-3prmt6ECE_lowerAcc'],
            [f'{stem2}_2.5e-3prmt6ECE_lowerAcc'],
            [f'{stem2}_1.75e-3prmt6ECE_lowerAcc'],
            [f'{stem2}_ireflm2_lowerAcc'],
            [f'{stem2}_7.5e-4prmt6ECE_lowerAcc'],

        ]
        labels = [
            r'prmt6LH = 5e-3',
            r'prmt6LH = 2.5e-3',
            r'prmt6LH = 1.75e-3',
            r'prmt6LH = 1e-3',
            r'prmt6LH = 7.5e-4',
        ]
        #"""

        colors = [ 'tab:green', 'tab:blue','tab:red', 'tab:purple', 'tab:grey', 'tab:orange']
        plt.rc('legend', fontsize = 10)
        ax.set_xlim([80,105])

    if case == '203619 4130 pwr scan':
        stem1 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/DIIID_{shotNum}_expSpectrum_2Zeff'
        stem2 = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/DIIID_{shotNum}_expSpectrum_2Zeff/scans/DIIID_{shotNum}_expSpectrum_2Zeff_second'
        includeWallRef = True

        targetDirs = [
            [f'{stem2}_0.75pwrFactor'],
            [f'{stem1}_second'],
            [f'{stem2}_1.25pwrFactor'],

        ]
        labels = [
            r'$P_{LH, tot} = 59$ kW',
            r'$P_{LH, tot} = 78$ kW',
            r'$P_{LH, tot} = 97.5$ kW',
        ]
        colors = [ 'tab:grey', 'tab:blue','tab:red', 'tab:purple', 'tab:grey']
        plt.rc('legend', fontsize = 10)
        ax.set_xlim([80,105])

    if includeWallRef:
        #ax.set_ylabel(r'$T_{rad}$ with ECE reflections (keV)')
        ax.set_ylabel(r'Radiation temperature T$_{rad}$ (keV)')
    else:
        ax.set_ylabel(r'$T_{rad}$ without ECE reflections (keV)')
    ax.set_xlabel(r'Frequency (GHz)')

    if shotNum == '203619.04130':
        minFreq = 80
    if shotNum == '203912.02700':
        minFreq = 75
    print(targetDirs)
    for i,dirs in enumerate(targetDirs):
        print(f'dirs: {dirs}')
        for j,targetDir in enumerate(dirs):
            print(f'target: {targetDir}')
            if j == 0:
                addECE(targetDir, ax, label = f'{labels[i]}', color = colors[i], includeWallRef = includeWallRef, minFreq = minFreq)#, style = styles[i])
            else:
                addECE(targetDir, ax, color = colors[i], includeWallRef = includeWallRef, minFreq = minFreq)


    DIIID_values = None
    sigma = None
    if shotNum == '203912.02700':
        DIIID_values = [2.5824008 , 1.39407992, 0.94645083, 0.81220359, 0.88051116, 0.89109635, 1.01204979, 1.18636823, 1.30664837, 1.44913375, 1.58996427, 1.7724458 , 1.95014727, 2.09026408, 2.23521471, 2.51402879, 2.65753579, 2.90507698, 3.10071349, 3.12329721, 3.1874361 , 3.14659691, 3.09538341, 3.06353974, 3.09011865, 3.22139716, 3.36753941, 3.46460414, 3.08602786, 3.25660014, 3.37270665, 3.22256875, 2.84639525, 3.09703612, 2.82084775, 2.99611259, 2.50551295, 2.79195881, 3.65390229, 3.07619691]
        sigma = [0.55028987, 0.28307492, 0.09646015, 0.03130944, 0.04082583, 0.05411086, 0.06373718, 0.07313851, 0.07793425, 0.08551253, 0.092725  , 0.09627454, 0.09089012, 0.07455167, 0.07081261, 0.0598761 , 0.0691742 , 0.05191937, 0.14060463, 0.24287122, 0.29219761, 0.25364348, 0.21264099, 0.20740688, 0.21601106, 0.23491451, 0.25725874, 0.26862234, 0.24190952, 0.25306058, 0.26157564, 0.25357226, 0.22021876, 0.2354593 , 0.21303467, 0.21872422, 0.1703857 , 0.20761023, 0.32986897, 0.21331525]
    if shotNum == '203619.04120':
        DIIID_values = [1.04244792, 1.40595269, 1.46621025, 1.53681779, 1.77099395, 1.69662797, 1.69903719, 1.80195522, 1.73125625, 1.06115818, 0.62294263, 0.37888351, 0.34158564, 0.34218457, 0.35121584, 0.42226839, 0.36570656, 0.39801067, 0.42528152, 0.51973337, 0.65875411, 0.80058628, 0.99458027, 1.18808925, 1.38460314, 1.59218562, 1.65122569, 1.56188107, 1.84427464, 1.81966686, 1.91930115, 1.86943126, 1.65424848, 1.78176737, 1.52368999, 1.54178166, 1.25577819, 1.39429832, 1.79460871, 1.63761687]
        sigma = [0.37833145, 0.47928944, 0.50052214, 0.52229393, 0.59195691, 0.58576584, 0.60000652, 0.62103432, 0.59744811, 0.46333417, 0.26068753, 0.10866313, 0.05944199, 0.04560371, 0.03528531, 0.04621909, 0.04186141, 0.04623201, 0.05026684, 0.05760876, 0.06987601, 0.08355413, 0.10504228, 0.12705854, 0.14640516, 0.17318632, 0.18412949, 0.1756306 , 0.20610285, 0.20266417, 0.21830384, 0.20829406, 0.19420874, 0.21549417, 0.18965581, 0.19358425, 0.15670006, 0.16962603, 0.21081436, 0.18824719]
    if shotNum == '203619.04130':
        DIIID_values = [1.54097414, 2.0470016 , 2.13048434, 2.21609354, 2.52367711, 2.44742584, 2.48486996, 2.61922908, 2.52230239, 1.8477025 , 1.0734036 , 0.52141368, 0.41335404, 0.39061141, 0.39694598, 0.50469536, 0.43668789, 0.47516528, 0.50995892, 0.59003371, 0.70966178, 0.8857165 , 1.11072564, 1.33231473, 1.55832803, 1.79727972, 1.87265325, 1.77506018, 2.08727193, 2.06216955, 2.18464136, 2.12412047, 1.90382946, 2.06328702, 1.77424574, 1.79718685, 1.46392119, 1.62084031, 2.07128787, 1.88276935]
        sigma = [0.21576479, 0.28587329, 0.2846531 , 0.2692441 , 0.28571698, 0.28720346, 0.30950746, 0.33040395, 0.32098365, 0.42575818, 0.26297823, 0.0771246 , 0.03311502, 0.01597273, 0.02041748, 0.04939581, 0.04108725, 0.04425885, 0.05656646, 0.03094726, 0.0329524 , 0.03808838, 0.04697888, 0.05584624, 0.06511637, 0.0731259 , 0.08131821, 0.07991054, 0.08974364, 0.09070159, 0.10038162, 0.09688853, 0.09744002, 0.11312548, 0.10269894, 0.10453832, 0.08540171, 0.09263628, 0.11132904, 0.09691198]
    if shotNum == '203917.02800':
        DIIID_values = [1.80684423, 0.97713965, 0.74930722, 0.71002197, 0.80390221, 0.8277722 , 0.94484293, 1.10375786, 1.20774841, 1.33600509, 1.46999645, 1.65462768, 1.81962347, 1.96521151, 2.10276437, 2.34134698, 2.47660923, 2.70755672, 2.89824033, 2.9351964 , 3.05754352, 3.11625242, 3.13883495, 3.10678887, 3.12562847, 3.24273086, 3.36878586, 3.45581889, 3.09331155, 3.24725628, 3.37756133, 3.23020887]
    if shotNum == '203917.03700':
        DIIID_values = [2.17810798, 1.14457178, 0.82409877, 0.74516374, 0.82109463, 0.83510202, 0.94940329, 1.10898256, 1.21686029, 1.35108495, 1.48144066, 1.67197084, 1.83846533, 1.98413777, 2.12280178, 2.36333299, 2.50468421, 2.71795607, 2.87802553, 2.87901521, 2.9358139 , 2.91673779, 2.90014362, 2.87531114, 2.92114258, 3.05473685, 3.18457508, 3.26851344, 2.93988132, 3.0838697 , 3.21016407, 3.07113123]
    if shotNum == '203619.04135':
        DIIID_values = [1.60368681, 2.13407207, 2.21663046, 2.29522824, 2.60867739, 2.53293085, 2.57472825, 2.71596956, 2.61736989, 2.05226946, 1.19208932, 0.54341048, 0.41945848, 0.39594352, 0.40976572, 0.53768569, 0.46423119, 0.50477809, 0.54427397, 0.61066061, 0.69726187, 0.88158321, 1.11442661, 1.34308743, 1.57791173, 1.82937312, 1.92797005, 1.83195543, 2.15573716, 2.1278038 , 2.25577521, 2.19547153, 1.97241962, 2.1404593 , 1.84182107, 1.86504126, 1.51844096, 1.67978466, 2.14379025, 1.94669354]
        sigma = [0.1341951 , 0.18883961, 0.18574601, 0.17839438, 0.18083471, 0.18444228, 0.19793077, 0.2074492 , 0.20107335, 0.24571943, 0.14181228, 0.04638657, 0.02498391, 0.00758356, 0.01583842, 0.03902015, 0.03233019, 0.03512308, 0.03041063, 0.01309039, 0.04299508, 0.04064882, 0.03726995, 0.03546975, 0.03245882, 0.03304584, 0.04748873, 0.05136775, 0.06062879, 0.05750338, 0.06273083, 0.065131  , 0.06560726, 0.07531843, 0.06734589, 0.06698281, 0.05551658, 0.06307615, 0.08169595, 0.06898947]
    if shotNum == '203619.04160':
        DIIID_values = [1.0947361 , 1.51075852, 1.53742182, 1.60393035, 1.84599531, 1.81839633, 1.85638893, 1.99375379, 2.0087893 , 1.99359763, 1.52894008, 0.65624166, 0.42941415, 0.35647714, 0.37664655, 0.52303809, 0.44625595, 0.55302513, 0.57984692, 0.61192405, 0.65078831, 0.70703232, 0.88569111, 1.15079796, 1.307881  , 1.53021681, 1.74554634, 1.62745142, 1.99958134, 1.97456038, 2.08317494, 2.02504802, 1.8125037 , 1.96905077, 1.6932652 , 1.7107197 , 1.39077437, 1.53138399, 1.94379079, 1.76558971]
        sigma = [0.09508576, 0.13271788, 0.10178398, 0.1072572 , 0.1392663 , 0.15347809, 0.16369736, 0.18033884, 0.18601771, 0.23529281, 0.35801619, 0.13636345, 0.05226064, 0.0131841 , 0.0079182 , 0.02455934, 0.01823262, 0.02702142, 0.02513472, 0.02903367, 0.02742721, 0.0584963 , 0.05754063, 0.14372389, 0.07144881, 0.06958455, 0.13149452, 0.07596993, 0.12544164, 0.13220492, 0.14590059, 0.14193654, 0.11939251, 0.13524652, 0.12319853, 0.1291451 , 0.10660744, 0.11200608, 0.13487308, 0.11914928]
    DIIID_freqs = np.concatenate([
        np.arange(16, dtype=float) + 83.5,
        np.arange(16, dtype=float) + 98.5,
        2 * np.arange(8, dtype=float) + 115.5
    ])
    print(DIIID_freqs)
    if DIIID_values is not None and sigma is not None:
        ax.errorbar(DIIID_freqs[:len(DIIID_values)], DIIID_values, yerr = sigma[:len(DIIID_values)],
                     marker = 'd', markersize = 5,zorder = 10, linestyle='none',color='k', label = 'DIII-D Measurements')
    elif DIIID_values is not None:
        ax.scatter(DIIID_freqs[:len(DIIID_values)], DIIID_values, label = 'DIII-D Measurements',color = 'k',marker = 'd', zorder = 10)
    ax.set_title(f'Shot {shotNum.split(".")[0]}, {shotNum.split(".")[1][1:]} ms', loc = 'right')
    leg = ax.legend(framealpha = 1, loc = 'best',ncol=1)
    #leg.set_zorder(100)
    ax.set_ylim(bottom = 0)

    fig.tight_layout()
    plt.savefig('203619_prmt6ECEScan.jpeg',dpi=300)
    plt.show()


main()
