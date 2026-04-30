
"""
Plots n_para, its toroidal and poloidal components, and n_||,acc
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.optimize import fsolve

print('past imports')

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import helperFunctions as helper
import netCDF4

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')

def plotNEvolution(rays = None, followRay = True):

    Nparas = np.copy(genray_nc.variables["wnpar"]) #n_|| of the ray at each point along the ray trace
    Nperps = np.copy(genray_nc.variables["wnper"]) #n_perp of the ray at each point along the ray trace
    Nphis = np.copy(genray_nc.variables["wn_phi"])#toroidal index of refraction

    if rays == None:
        rays = range(len(Nparas))

    if followRay:
        sene = np.copy(genray_nc.variables["sene"])*1e6
        sbtot = np.copy(genray_nc.variables["sbtot"])/1e4
        sb_r = np.copy(genray_nc.variables["sb_r"])/1e4
        sb_z = np.copy(genray_nc.variables["sb_z"])/1e4
        sb_theta = np.sqrt(sb_r**2 + sb_z**2)
        sb_phi = np.copy(genray_nc.variables["sb_phi"])/1e4

        #dielectric tensor components along the ray
        cweps11 = genray_nc.variables["cweps11"][:]
        cweps11 = cweps11[0] + 1j*cweps11[1]
        cweps12 = genray_nc.variables["cweps12"][:]
        cweps12 = cweps12[0] + 1j*cweps12[1]
        cweps21 = genray_nc.variables["cweps21"][:]
        cweps21 = cweps21[0] + 1j*cweps21[1]
        cweps33 = genray_nc.variables["cweps33"][:]
        cweps33 = cweps33[0] + 1j*cweps33[1]

    delpwr= np.copy(cqlrf_nc.variables["delpwr"]) #power in the ray at each point
    radialVariable = (np.copy(genray_nc.variables["spsi"])) #rho_pol of the ray at each point along the ray trace
    wr = np.copy(genray_nc.variables["wr"])/100#major radius of ray
    w_theta_pol = np.copy(genray_nc.variables["w_theta_pol"])*np.pi/180#poloidal angle of ray in radians

    w = 2*np.pi*4.6e9#*u.rad/u.s #angular frequency of the wave

    fig,axes = plt.subplots()#,figsize = (7,7), sharex = True)
    maxDampingToPlot = .9 #when the ratio of ray power to ray starting power is below this number, the trace ends

    ryain, Tein = helper.getCQLTe()
    ryain, nein = helper.getCQLne()
    dampingCondition = getDampingCondition(Tein,nein)


    #axes.plot(ryain, np.sign(Nparas[0][0])*dampingCondition, label = r'$N_{ELD}$', lw = 3, color = 'r', linestyle = 'dashed')

    for j in range(len(rays)):
        rayNum = rays[j]

        delpwrRatios = delpwr[rayNum]/np.max(delpwr[rayNum])
        powerIndex = helper.findNearestIndex(1-maxDampingToPlot, delpwrRatios) #find the index of the last ray point we want to plot

        if followRay:
            S, D, P = np.real(cweps11[rayNum]), np.imag(cweps21[rayNum]), np.real(cweps33[rayNum])
            N_acc = np.sqrt((-D**2*(P+S) + 2*np.sqrt(D**2*P*S*(D**2-(P-S)**2))+S*(P-S)**2)/((P-S)**2))
            #N_acc_1 = np.sqrt((-D**2*(P+S) - 2*np.sqrt(D**2*P*S*(D**2-(P-S)**2))+S*(P-S)**2)/((P-S)**2))

            B_tot_abs = np.abs(sbtot[rayNum])
            B_tor = sb_phi[rayNum]
            #for normal Ip, the resulting B_pol is pointing backwards along the poloidal angle
            #so we need a minus sign to account for that
            B_pol = sb_theta[rayNum]

            Nthetas = (Nparas[rayNum]-Nphis[rayNum]*(B_tor/B_tot_abs))*(B_tot_abs/B_pol)

            points = np.array([radialVariable[rayNum][:powerIndex],Nparas[rayNum][:powerIndex]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            norm = plt.Normalize(0, 1)

            # Create a continuous norm to map from datline co points to colors
            lc = LineCollection(segments, norm = norm,cmap=plt.cm.jet, zorder = 5)
            # Set the values used for colormapping
            lc.set_array(delpwrRatios[:powerIndex])
            lc.set_linewidth(3)
            #axes.add_collection(lc)

            poloidalContribution = Nparas[rayNum][:powerIndex] - (Nphis[rayNum][:powerIndex] *B_tor[:powerIndex]/ np.abs(B_tot_abs[:powerIndex]))

            axes.plot(radialVariable[rayNum][:powerIndex],Nparas[rayNum][:powerIndex], label = r'$N_{||}$', lw = 3, color = 'k')


            axes.plot(radialVariable[rayNum][:powerIndex],poloidalContribution, label = r'$N_\theta \frac{B_\theta}{|B|}$', lw = 3, color = 'darkturquoise', linestyle = 'dashed')
            #axes.axhline(0,color = 'grey',lw = 1.5, linestyle = 'dashdot')
            axes.plot(radialVariable[rayNum][:powerIndex],Nphis[rayNum][:powerIndex] *B_tor[:powerIndex]/ np.abs(B_tot_abs[:powerIndex]), label = r'$N_\phi \frac{B_\phi}{|B|}$', lw = 3, color ='orange', linestyle = 'dashed')
            axes.plot(radialVariable[rayNum][:powerIndex],np.sign(Nparas[rayNum][0])*N_acc[:powerIndex], label = r'$N_{acc}$', lw = 3, color ='g', linestyle = 'dotted')
            
            
            
            #axes.plot(radialVariable[rayNum][:powerIndex],np.sign(Nparas[rayNum][0])*N_acc_1[:powerIndex], label = r'$N_{||, acc}$', lw = 3, color ='red', linestyle = 'dotted')


        else:
            #"""
            points = np.array([radialVariable[rayNum][:powerIndex],Nparas[rayNum][:powerIndex]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            norm = plt.Normalize(0, 1)

            # Create a continuous norm to map from datline co points to colors
            lc = LineCollection(segments, norm = norm,cmap=plt.cm.jet, zorder = 5)
            # Set the values used for colormapping
            lc.set_array(delpwrRatios[:powerIndex])
            lc.set_linewidth(3)
            axes.add_collection(lc)
            #axes.plot(radialVariable[rayNum],Nparas[rayNum])
            #"""
            #axes.scatter(np.sin(w_theta_pol[rayNum][:powerIndex]),Nparas[rayNum][:powerIndex], c = delpwrRatios[:powerIndex]*100, cmap = plt.cm.jet)

        

        #axes.set_ylim([-5,.5])
        axes.set_xlim([.76,1])
        axes.set_ylim([-4.65,0])
        #axes.set_ylim([np.nanmax(Nparas[:,0])*1.5, np.nanmin(Nparas[:,0])*.5])
        #axes.set_xlim([0,2*np.pi])

    axes.text(.792,-3.67, 'Ray damps', rotation = 40, fontsize = 14, ha='center', va='center')
    axes.text(.9714,-2.88, 'HFS launch', rotation = -5, fontsize = 14, ha='center', va='center')
    
    axes.annotate(
                "",
                xy=(.9467,-2.98),
                xytext=(.9948,-3.1),
                arrowprops=dict( linewidth = 2, 
                                linestyle = '-', 
                                color = 'k',
                                arrowstyle = 'simple',
                                joinstyle='miter',   # sharp corners
                                capstyle='butt'  
                                )
                )

    axes.set_xlabel(r"$\rho_{pol}$")#
    axes.set_ylabel('Index of refraction')
    axes.legend(loc = 'lower right', ncol = 2)
    fig.tight_layout()
    plt.savefig('180403_Npara_evo.jpeg',dpi=300)

    plt.show()

def getDampingCondition(Te, ne):
    e = 1.602e-19
    m_e = 9.109e-31
    eps_0 = 8.854e-12
    c=2.99e8
    dx = 0.01
    w_pes = np.sqrt(ne*e**2/(eps_0*m_e))

    condition = np.zeros(len(Te))

    v_ths = np.sqrt(2*Te*1.602e-16/m_e)


    def damping(Npara, w_pe,v_th):

        x = c/(Npara*v_th)
        return 2*np.sqrt(np.pi)*w_pe*dx*Npara*x**3*np.exp(-x**2)/c - 1
    """
    for i in range(len(condition)):
        from scipy.optimize import brentq
        # Use brentq (Brent's method) to find the root within an interval [a, b] where the function changes sign
        condition[i] = brentq(lambda x: damping(x, w_pes[i], v_ths[i]), .1, 20)
    """
    return 6.4/np.sqrt(Te)#condition



plotNEvolution(followRay=True, rays = [19])
#plotNEvolution(followRay=False)
