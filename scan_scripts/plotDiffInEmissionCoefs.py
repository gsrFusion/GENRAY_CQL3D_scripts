
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.collections import LineCollection

import matplotlib
import os, sys
from scipy.interpolate import interp1d
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
plt.rc('axes', titlesize = 14)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

#adds the ray traces to ax
def main():
    
    shotNum = '203619.04135'

    stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/redemption/DIIID_{shotNum}'

    t1 = f'{stem}_expSpectrum_2Zeff_0.005prmt6_1e-8prmt4_60nnkpar_prmt4ECE1e-7_n0.1Fld_2ireflm'
    t2 = f'{stem}_expSpectrum_2Zeff_0.005prmt6_1e-8prmt4_60nnkpar_prmt4ECE1e-7_n0.05Fld_2ireflm'

    ECE_nc_1 = netCDF4.Dataset(f'{t1}/genray_ece.nc','r')
    ECE_nc_2 = netCDF4.Dataset(f'{t2}/genray_ece.nc','r')

    wr_1 = ECE_nc_1.variables['wr_em_nc'][:]
    wr_2 = ECE_nc_2.variables['wr_em_nc'][:]

    ws_1 =  ECE_nc_1.variables['wsn_nc'][:]
    ws_2 =  ECE_nc_2.variables['wsn_nc'][:]

    emis_1 = ECE_nc_1.variables['wj_emis_nc'][:]
    emis_2 = ECE_nc_2.variables['wj_emis_nc'][:]

    freqs_1 = ECE_nc_1.variables['wfreq_nc'][:]
    freqs_2 = ECE_nc_2.variables['wfreq_nc'][:]

    assert np.equal(freqs_1, freqs_2).all()

    for i, freq in enumerate(freqs_1):

        fig,ax = plt.subplots()
        ax.plot(ws_1[i], emis_1[i])
        ax.plot(ws_2[i], emis_2[i])
        plt.show()
    


main()
