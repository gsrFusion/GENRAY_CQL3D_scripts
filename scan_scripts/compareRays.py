
"""
Plots n_para, its toroidal and poloidal components, and n_||,acc
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.collections import LineCollection
from scipy.optimize import fsolve

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 13)

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
from scipy.interpolate import interp1d
import shotToEqdsk
#from omfit_classes import omfit_eqdsk

def plotNEvolution(targetDir, ax, machine, rays = None, color = None, labelSuffix = ''):
    shotNum = targetDir.split('/')[6].split('_')[-1]
    eqdskName = shotToEqdsk.getEqdskName(shotNum, machine = machine)
    #geqdsk = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')

    #delta_gfile = geqdsk['fluxSurfaces']['geo']['delta']
    #rho_pol_gfile = geqdsk['fluxSurfaces']['levels']
    #delta_func = interp1d(rho_pol_gfile, delta_gfile, bounds_error=False, )

    cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
    #cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')   

    #rya = cql_nc.variables["rya"][:]

    Nparas = np.copy(genray_nc.variables["wnpar"]) #n_|| of the ray at each point along the ray trace
    Nr = np.copy(genray_nc.variables["wn_r"]) #n_|| of the ray at each point along the ray trace
    Nz = np.copy(genray_nc.variables["wn_z"]) #n_|| of the ray at each point along the ray trace
    Nperps = np.copy(genray_nc.variables["wnper"]) #n_per of the ray at each point along the ray trace
    Nphis = np.copy(genray_nc.variables["wn_phi"])#toroidal index of refraction

    if rays == None:
        rays = range(len(Nparas))

    sene = np.copy(genray_nc.variables["sene"])*1e6
    ste = np.copy(genray_nc.variables["ste"])
    szeff = np.copy(genray_nc.variables["szeff"])
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
    cweps13 = genray_nc.variables["cweps13"][:]
    cweps13 = cweps13[0] + 1j*cweps13[1]
    cweps21 = genray_nc.variables["cweps21"][:]
    cweps21 = cweps21[0] + 1j*cweps21[1]
    cweps22 = genray_nc.variables["cweps22"][:]
    cweps22 = cweps22[0] + 1j*cweps22[1]
    cweps23 = genray_nc.variables["cweps23"][:]
    cweps23 = cweps23[0] + 1j*cweps23[1]
    cweps33 = genray_nc.variables["cweps33"][:]
    cweps33 = cweps33[0] + 1j*cweps33[1]

    delpwr= np.copy(cqlrf_nc.variables["delpwr"])#power in the ray at each point
    radialVariable = (np.copy(genray_nc.variables["spsi"])) #rho_pol of the ray at each point along the ray trace
    wr = np.copy(genray_nc.variables["wr"])/100#major radius of ray
    thetas = np.copy(genray_nc.variables["w_theta_pol"])*np.pi/180#poloidal angle of ray in radians
    ws= np.copy(genray_nc.variables["ws"])/100 #poloidal length along ray
    wphi= np.copy(genray_nc.variables["wphi"]) #toroidal angle
    fluxn= np.copy(genray_nc.variables["fluxn"]) #toroidal angle
    salphal= np.copy(genray_nc.variables["salphal"]) #Linear damping wavenumber
    salphac= np.copy(genray_nc.variables["salphac"]) #collisional damping wavenumber
    seikon= np.copy(genray_nc.variables["seikon"]) #eikonal
    w = 2*np.pi*4.6e9#*u.rad/u.s #angular frequency of the wave
   
    maxDampingToPlot = .9 #when the ratio of ray power to ray starting power is below this number, the trace ends

    m_e = 9.109e-31
    m_D = 3.343e-27
    e = 1.602e-19
    eps_0 = 8.854e-12
    c = 2.99e8

    for j in range(len(rays)):
        rayNum = rays[j]
        deltas = 1#delta_func(radialVariable[rayNum])
        deltaHat = np.arcsin(deltas)
        Gamma = np.sin(deltaHat*np.sin(thetas[rayNum]) + thetas[rayNum])*(deltaHat*np.cos(thetas[rayNum])+1)

        delpwrRatios = delpwr[rayNum]/np.max(delpwr[rayNum])
        powerIndex = helper.findNearestIndex(1-maxDampingToPlot, delpwrRatios) #find the index of the last ray point we want to plot
        bounceIndex = helper.findBounceIndex(radialVariable[rayNum],bounceToFind = 3)

        endIndex = powerIndex#min(bounceIndex, powerIndex)

        if 'PT' in labelSuffix:
            for p in range(len(radialVariable[rayNum])):
                if radialVariable[rayNum][p] > 1.02 and radialVariable[rayNum][p+1] < radialVariable[rayNum][p]:
                    endIndex = p
                    break

        rhos = radialVariable[rayNum][:endIndex]

        S, D, P = np.real(cweps11[rayNum]), np.imag(cweps21[rayNum]), np.real(cweps33[rayNum])
        N_acc = np.sqrt((-D**2*(P+S) + 2*np.sqrt(D**2*P*S*(D**2-(P-S)**2))+S*(P-S)**2)/((P-S)**2))

        B_tot_abs = np.abs(sbtot[rayNum])
        B_tor = sb_phi[rayNum]
        B_pol = -sb_theta[rayNum]

        """
        increasingThetas = np.copy(thetas[rayNum][:endIndex])
        poloidalDist = ws[rayNum][:endIndex]

        def advanceTheta(t):
            for p in range(len(t)):
                if p == 0:
                    continue
                else:
                    if t[p-1] > 6.2 and t[p] < .01:
                        t[p:] += 2*np.pi
                
        advanceTheta(increasingThetas)

        dB_dl = (B_pol[:endIndex][1:] - B_pol[:endIndex][:-1])/(poloidalDist[1:] - poloidalDist[:-1])
        dtheta_dl = (increasingThetas[:endIndex][1:] - increasingThetas[:endIndex][:-1])/(poloidalDist[1:] - poloidalDist[:-1])

        db_dtheta = dB_dl/dtheta_dl

        thetaCenters = (increasingThetas[1:] + increasingThetas[:-1])/2

        B_pol_derivFunc = interp1d(thetaCenters,db_dtheta, bounds_error=False, fill_value = 'extrapolate')
        B_pol_deriv = B_pol_derivFunc(increasingThetas)
        """
        """

        w_pes = np.sqrt(sene[rayNum]*e**2/(eps_0*m_e))
        w_pDs = np.sqrt(sene[rayNum]*e**2/(eps_0*m_D))
        w_ces = e*B_tot_abs/m_e

        freqRatio = w_pes/w_ces

        maxEpsElement = np.zeros(endIndex)
        for l in range(len(maxEpsElement)):
            maxEpsElement[l] = np.max(np.abs([cweps11[rayNum][l], cweps12[rayNum][l], cweps13[rayNum][l], cweps22[rayNum][l], cweps23[rayNum][l], cweps33[rayNum][l]]))

        electrostaticRatio = (Nparas[rayNum]**2 + Nperps[rayNum]**2)[:endIndex] / maxEpsElement
        k_paras = w*Nparas[rayNum]/c
        k_perps = w*Nperps[rayNum]/c
        eps_para = -(w_pes**2/w**2)
        eps_perp = 1 + (w_pes/w_ces)**2 - (w_pDs/w)**2
        k_phis = Nphis[rayNum]*(w/c)
        Nthetas = (Nparas[rayNum]-Nphis[rayNum]*(B_tor/B_tot_abs))*(B_tot_abs/B_pol)
        k_thetas = Nthetas*(w/c)

        denom = (2*w_pes**2*k_paras**2/w**3)*(1 + w_pDs**2/(w**2*eps_perp))

        D = eps_para*k_paras**2 + eps_perp*k_perps**2
        """
        """
        kpara_gamma_component = 2*k_paras*eps_para*((Gamma/A) * (2*B_tor*k_phis/B_tot_abs - k_thetas*B_pol*B_tor**2/(B_tot_abs**3) - k_phis*B_tor**3/(B_tot_abs**3)))
        kpara_deriv_component = 2*k_paras*eps_para*(1/B_tot_abs)*(k_thetas - k_thetas *B_pol**2/B_tot_abs**2 - k_phis*B_tor*B_pol/B_tot_abs**2)
        kpara_gamma_component = kpara_gamma_component[:powerIndex]
        kpara_deriv_component = kpara_deriv_component[:powerIndex]*B_pol_deriv

        epsperp_gamma_component = (eps_para*k_paras**2/eps_perp)*(2*w_pes**2/w_ces**2)*B_tor**2*Gamma/(A*B_tot_abs**2)
        epsperp_deriv_component = (eps_para*k_paras**2/eps_perp)*(2*w_pes**2/w_ces**2)*B_pol/(B_tot_abs**2)
        epsperp_gamma_component = epsperp_gamma_component[:powerIndex]
        epsperp_deriv_component = epsperp_deriv_component[:powerIndex]*B_pol_deriv

        kperp_gamma_component = eps_perp*2*k_phis*(k_phis*B_pol/B_tot_abs - k_thetas)*((B_pol/B_tot_abs)*(Gamma/A) - B_pol*B_tor**2*Gamma/(A*B_tot_abs**3))
        kperp_deriv_component = eps_perp*2*k_phis*(k_phis*B_pol/B_tot_abs - k_thetas)*((1/B_tot_abs) - B_pol**2/B_tot_abs**3)
        kperp_gamma_component = kperp_gamma_component[:powerIndex]
        kperp_deriv_component = kperp_deriv_component[:powerIndex]*B_pol_deriv

        total = kpara_gamma_component + kpara_deriv_component + epsperp_gamma_component + epsperp_deriv_component + kperp_gamma_component + kperp_deriv_component
        """

        #dmdt = total/denom[:endIndex]

        thetaComponent = (Nparas[rayNum]-Nphis[rayNum]*(B_tor/B_tot_abs))
        phiComponent = Nphis[rayNum]*(B_tor/B_tot_abs)
        #ax.plot(ws[rayNum][:endIndex],(Nparas)[rayNum][:endIndex], lw = 3,color = color, label = f'{labelSuffix}')
        #ax.set_ylabel(r'N||')
        #ax.plot(ws[rayNum][:endIndex], Nparas[rayNum][:endIndex], lw = 3,color = color, label = f'{labelSuffix}')
        ax.plot(rhos[:endIndex], Nparas[rayNum][:endIndex], lw = 3, color = color, label = rf'$N_{{||}}$, {labelSuffix}')
        
        #"""
        threshold = .05
        radialDistance = threshold

        arrowColor = ''
        if color == 'mediumblue':
            arrowColor = 'royalblue'
        if color == 'crimson':
            arrowColor = 'orangered'
        if color == 'darkturquoise':
            arrowColor = 'teal'
        if color == 'orange':
            arrowColor = 'goldenrod'


        from matplotlib.patches import FancyArrowPatch
        for i in range(len(rhos)-1):
            if radialDistance >= threshold:
                dx = rhos[i+1] - rhos[i]
                dy = Nparas[rayNum][:endIndex][i+1] - Nparas[rayNum][:endIndex][i]
                arrow = FancyArrowPatch(
                    (rhos[i], Nparas[rayNum][:endIndex][i]), (rhos[i]+dx, Nparas[rayNum][:endIndex][i]+dy),
                    arrowstyle='->', color=arrowColor,
                    mutation_scale=15, lw=1.5,zorder=10
                )
                ax.add_patch(arrow)

                dy1 = thetaComponent[:endIndex][i+1] - thetaComponent[:endIndex][i]
                arrow1 = FancyArrowPatch(
                    (rhos[i], thetaComponent[:endIndex][i]), (rhos[i]+dx, thetaComponent[:endIndex][i]+dy1),
                    arrowstyle='-|>', color=arrowColor,
                    mutation_scale=20, lw=0,zorder=10
                )
                ax.add_patch(arrow1)
            
                radialDistance = 0

            else:
                dist = np.sqrt((rhos[i]- rhos[i-1])**2)
                radialDistance += dist
        #"""
        
        ax.plot(rhos, thetaComponent[:endIndex], lw = 3, linestyle = 'dashed', color = color, label = rf'$N_{{\theta}}B_{{\theta}}/|B|$, {labelSuffix}')
        #ax.plot(radialVariable[rayNum][:endIndex], phiComponent[:endIndex], lw = 3, linestyle = 'dotted', color = color, label = rf'$N_{{\phi}}B_{{\phi}}/|B|$, {labelSuffix}')
        
        #ax2.set_ylabel(r'$dm/dt$')
        

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

    for i in range(len(condition)):
        from scipy.optimize import brentq
        # Use brentq (Brent's method) to find the root within an interval [a, b] where the function changes sign
        condition[i] = brentq(lambda x: damping(x, w_pes[i], v_ths[i]), .1, 20)

    return condition


def makePlot():
    fig,ax = plt.subplots(figsize = (6.4,5.25))#figsize=(6.5,10))#, nrows = 3)
    ax.set_ylim([-4,0.5])
    #ax.set_ylim([-.5,5])
    machine = 'NTPT'


    if machine == 'NTPT':
        fakeMachine = 'DIIID'

        if fakeMachine == 'MANTA':

            power =10
            Npara = -2.1
            thgrill = 140

            targetDirs = [
                        f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeMachine}.PT05/{machine}_{fakeMachine}.PT05_n{np.abs(Npara)}Npara_{thgrill}thgrill_{power}MW',
                        #f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeMachine}.NT05/{machine}_{fakeMachine}.NT05_n{np.abs(Npara)}Npara_{thgrill}thgrill_{power}MW',
                        ] 
            labels = ['PT', 'NT']

            colors = ['tab:blue', 'tab:orange']

        if fakeMachine == 'DIIID':

            power =1
            #"""
            Npara = -2.8
            grillHeight = -.25

            targetDirs = [
                        f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeMachine}.147634PT/{machine}_{fakeMachine}.147634PT_n{np.abs(Npara)}Npara_{grillHeight}grillHeight_{power}MW',
                        f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeMachine}.147634NT/{machine}_{fakeMachine}.147634NT_n{np.abs(Npara)}Npara_{grillHeight}grillHeight_{power}MW',
                        ] 
            labels = ['PT', 'NT']
            
            colors = ['mediumblue', 'crimson']


            """
            Npara = 2.5
            grillHeight = 0.4
            shot = '193765'

            targetDirs = [
                        f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeMachine}.{shot}PT/{machine}_{fakeMachine}.{shot}PT_p{np.abs(Npara)}Npara_{grillHeight}grillHeight_{power}MW',
                        f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{fakeMachine}.{shot}NT/{machine}_{fakeMachine}.{shot}NT_p{np.abs(Npara)}Npara_{grillHeight}grillHeight_{power}MW',
                        ] 
            labels = ['PT', 'NT']
            
            colors = ['darkturquoise', 'orange']
            """


    if machine == 'DIIID':
        shot = 203619
        time = '.04135'
        Npara = -2.7
        shotNum = f'{shot}{time}'
        stem = f'/home/grantr/symlinks/genray_batch/DIIID_shots/DIIID_{shotNum}/numRaysTest/DIIID_{shotNum}'
        targetDirs =[
            f'{stem}_expSpectrum_id2_1Zeff_10000nrelt_0.005prmt6_0.000002prmt4_oneRay',
            #f'{stem}_expSpectrum_id2_1Zeff_10000nrelt_0.005prmt6_0.000002prmt4_4nthin_30nnkpar',
            #f'{stem}_expSpectrum_1Zeff_10000nrelt_0.005prmt6_0.000002prmt4_4nthin_30nnkpar',
            #f'{stem}_expSpectrum_129Eqdsk_1Zeff_10000nrelt_0.005prmt6_0.000002prmt4_4nthin_30nnkpar',
            ]
    
        labels = [
                'id = 2,', 
                #'id = 16_SAM', 
                #'id = 16_SAM, 129x129 eqdsk', 
                ]

        colors = ['tab:blue', 'tab:red', 'tab:green','tab:purple','tab:pink','tab:brown','tab:grey']



    for i in range(len(targetDirs)):
        targetDir = targetDirs[i]
        #ray  20
        #90 for 147, 10 for 193
        plotNEvolution(targetDir, ax, machine, color = colors[i], rays = [90],labelSuffix = labels[i])

    ax.text(1.01, -1.2, "LCFS", size=16, rotation=-90.,
         ha="center", va="center"
         )
    ax.text(.99, -4.6, "147634-like", size=16,
         ha="center", va="center"
         )
    ax.set_xlabel(r"$\rho_{pol}$")
    #ax.set_xlabel(r"$\rho_{pol}$")
    ax.set_ylabel(r'Refractive Index')
    #ax.set_ylabel(r'delpwr')
    ax.set_xlim(left= 0.55)
    ax.grid()
    ax.axvline(1, color ='k', lw = 2, linestyle= 'dashed')
    ax.legend(loc='lower center', bbox_to_anchor=(0.5,.99),ncol=2,labelspacing=0.3)
    
    fig.tight_layout()
    plt.savefig(f'toka_147634_evoCompare_n2.8_-0.25Height.jpeg',dpi=300)

    plt.show()


makePlot()
