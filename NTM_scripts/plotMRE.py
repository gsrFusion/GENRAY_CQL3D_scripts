#analysis based off of https://iopscience.iop.org/article/10.1088/0029-5515/46/4/006/pdf


import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import numpy as np
import matplotlib.pyplot as plt
from omfit_classes import omfit_eqdsk, utils_fusion
import netCDF4
import shotToEqdsk
import helperFunctions as helper
import getInputFileDictionary
from scipy.interpolate import interp1d
import xarray as xr
import pickle
from scipy.interpolate import RectBivariateSpline as RBS
from scipy.interpolate import interp1d

plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.rc('axes', labelsize = 18)
plt.rc('axes', titlesize = 14)
plt.rc('legend', fontsize = 13)

#relative alignment is in units of DeltaR/deltaLH
def generateGrowthRates(info=None, ws=None, ax_w=None, K_1Dict = None,
                        plotLabel=None, GENCQLResultDict = None, relAlign = 0, waveType='None',index=0):
    
    q2Index = helper.findNearestIndex(2, info.q_prof)
    rho_pol = info.rho_pol
    q2Rho = rho_pol[q2Index]
    #print(f'q = 2 happens at rho_pol = {q2Rho}')
    #q2FluxWidths_psiN = GENCQLResultDict['q2FluxWidths_psiN']

    if waveType == 'LH':
        JPeakCenters = GENCQLResultDict['JPeakCenters_rhop']
        JPeakFWHMs = GENCQLResultDict['JPeakFWHMfits_m']
        JPeakMags = GENCQLResultDict['JPeakMagsfits_Aperm2']
        JPeakNparas = GENCQLResultDict['JPeakNparas']
        FWHM_func = interp1d(JPeakCenters, JPeakFWHMs, kind = 'linear', bounds_error = False, fill_value = np.inf)
        mag_func = interp1d(JPeakCenters, JPeakMags, kind = 'linear', bounds_error = False, fill_value = 0)

        Npara_func = interp1d(JPeakCenters, JPeakNparas, kind = 'linear', bounds_error = False, fill_value = 0)

        R_LFSmidplane = info.R_LFSmidplane
        R_q2 = R_LFSmidplane[q2Index]
        R_offset = R_LFSmidplane - R_q2
        alignments = R_offset / FWHM_func(rho_pol)
        #mask for where is there actually possible current drive
        mask = np.where((rho_pol > np.min(JPeakCenters))*(rho_pol < np.max(JPeakCenters)))
        masked_alignments = alignments[mask]

        if np.max(masked_alignments) < relAlign or relAlign < np.min(masked_alignments):
            #pass
            print(f'desired alignment is not possible given input CD fits, returning nan')
            return np.nan

        dampingLoc_R = R_LFSmidplane[mask][helper.findNearestIndex(relAlign, masked_alignments)]
        indexOfDampingLoc = helper.findNearestIndex(dampingLoc_R, R_LFSmidplane)

        damping_rho = rho_pol[indexOfDampingLoc]
        """
        fig,ax = plt.subplots()
        ax.plot(rho_pol, alignments)
        ax.axvline(damping_rho, color = 'k')
        ax.axvline(q2Rho, color = 'r')
        plt.show()
        """

        del_CD = FWHM_func(damping_rho)
        J_CD_peak = mag_func(damping_rho)
        Npara = Npara_func(damping_rho)

        print(f'for relative alignment {relAlign}, need N|| = {Npara}')
            
    elif waveType == 'EC':
        del_CD = 6.15/100
        J_CD_peak = 8.065e4
        plotLabel = rf'1 MW EC, $\Delta R/\delta_{{EC}} = 0$'
    else:
        J_CD_peak = 0
        del_CD = np.inf


    a_2 = 3.6

    F = 1 - 2.43*relAlign + 1.4*relAlign**2 - .23*relAlign**3
    delDeltar = -(5/32)*np.pi**(-3/2)*a_2*(info.L_q_prof[q2Index]/del_CD)*(J_CD_peak/info.J_para_prof[q2Index])*F #function of r
    dwdt = np.zeros(len(ws))
    q2Index = helper.findNearestIndex(2, info.q_prof)

    #Delta_prime_re = StuartXarray.Re_DeltaPrime.sel(m=2,n=1).item()
    #Delta_prime_im = StuartXarray.Im_DeltaPrime.sel(m=2,n=1).item()

    #Di = StuartXarray.Di.sel(m=2,n=1).item()

    widthPoints = K_1Dict['widthPoints']
    alignmentPoints = K_1Dict['alignmentPoints']
    K_1_matrix = K_1Dict['K_1_width_alignment']
    K_1_func = RBS(widthPoints, alignmentPoints, K_1_matrix)

    CD_JRatio_ofInterest = J_CD_peak/info.J_BS_prof[q2Index]
    para_JRatio_ofInterest = (info.J_BS_prof/info.J_para_prof)[q2Index]
    w_marg_ofInterest = info.w_marg[q2Index]
    L_q_ofInterest =info. L_q_prof[q2Index]

    us = ws/del_CD
    K_1s = np.zeros(len(us))
    if relAlign is not None:
        K_1s = K_1_func(us, relAlign).ravel()

    for i in range(len(ws)):
        w = ws[i]
        #fluxWidth = q2FluxWidths_psiN[i]
        Deltar = 0#2*(q2Rho**2)*Delta_prime_re * (fluxWidth/2)**(-1+np.sqrt(-4*Di))*np.sqrt(-4*Di)
        #Deltar_old = Delta_prime_re * (fluxWidth/2)**(np.sqrt(-4*Di))

        #print(f'new - old: {Deltar-Deltar_old}')
        #K_1 = 0.007 + 0.863*us -0.376*us**2+0.0464*us**3
        #this is really tau_r/ r * dw/dt
        dwdt[i] = (Deltar + delDeltar + a_2*para_JRatio_ofInterest*(L_q_ofInterest/w)*(1 - w_marg_ofInterest**2/(3*w**2) - K_1s[i]*CD_JRatio_ofInterest))

    if len(ws) > 1:
        if relAlign < 0 :
            ax_w.plot(ws, dwdt, label = f'{plotLabel}', lw = 3,color = 'tab:olive')
        else:
            ax_w.plot(ws, dwdt, label = f'{plotLabel}', lw = 3)

    else:
        return dwdt

class plasmaInfo:
    def __init__(self,targetDir, gfile,rho_pol):
        #cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
        #rya = np.ma.getdata(cql_nc.variables["rya"][:])

        self.rho_pol = rho_pol
        omfit_rho_pol = np.sqrt(gfile['fluxSurfaces']['levels'])

        self.q_prof = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['q'])(rho_pol)

        ryain, Tein = helper.getCQLTe(targetDir)
        ryain, TDin = helper.getCQLTD(targetDir)
        ryain, nein = helper.getCQLne(targetDir)
        ryain, nDin = helper.getCQLnD(targetDir)
        ryain, nCin = helper.getCQLnC(targetDir)
        ryain, Zeffin = helper.getCQLZeff(targetDir)

        Te_prof = interp1d(ryain, Tein)(rho_pol)
        TD_prof = interp1d(ryain, TDin)(rho_pol)
        ne_prof = interp1d(ryain, nein)(rho_pol)
        nD_prof = interp1d(ryain, nDin)(rho_pol)
        nC_prof = interp1d(ryain, nCin)(rho_pol)

        pressure = (ne_prof*Te_prof + nD_prof*TD_prof + nC_prof*TD_prof)*1.602e-16

        self.R_LFSmidplane = helper.convertRhopolToRmidplane(rho_pol, targetDir = targetDir, side = 'LFS')
        #R_LFSmidplane = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['a'])(rho_pol)

        eps = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['geo']['eps'])(rho_pol)
        avgBp = interp1d(omfit_rho_pol, gfile['fluxSurfaces']['avg']['Bp'])(rho_pol)

        self.J_BS_prof = utils_fusion.sauter_bootstrap(gEQDSKs = gfile, psi_N = self.rho_pol**2, 
                        Ti = np.array([TD_prof*1e3]), ne = np.array([ne_prof]), Te = np.array([Te_prof*1e3]),
                        charge_number_to_use_in_ion_collisionality = 'Koh', charge_number_to_use_in_ion_lnLambda = 'Koh',
                        Zis=[1,6], nis = np.array([[nD_prof], [nC_prof]]), R0 = 1.6955, p = np.array([pressure]), version = 'osborne')[0]#A/m^2, function of r
        self.J_para_prof = interp1d(omfit_rho_pol, gfile.surfAvg('Jpar', interp = 'cubic'))(self.rho_pol)#A/m^2, function of r

        m_D = 3.3e-27 #kg
        q_D = 1.6e-19
        keV_to_Joule = 1.60218e-16

        v_th = np.sqrt(2*keV_to_Joule*Te_prof/m_D) #including factor of two since we care of perp velocity
        poloidal_larmor_radius_D = m_D*v_th/(q_D*avgBp)
        #marginal island size
        self.w_marg = 2*np.sqrt(eps)*(poloidal_larmor_radius_D)#function of r
        
        dq_dr = np.gradient(np.asarray(self.q_prof), self.R_LFSmidplane)
        self.L_q_prof = self.q_prof / dq_dr #function of r

def generateComparisonPlot():
    #fig_w,ax_w = plt.subplots(figsize = (7,5.8))#7,5.8
    fig_w,ax_w = plt.subplots(figsize = (7,5))#7,5.8

    machine = 'DIIID'
    shot = '180403'
    time = '.04400'

    targetSuffix = '_1MW'#'_0.15Width_0.2deltaT'
    storageStem = '/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage'

    targetDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{time}/{machine}_{shot}{time}_n3.0Npara{targetSuffix}'
    gfile = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/g{shot}{time}')
    rho_pol = np.load(f'{storageStem}/NTM_rho_pol.npy')
    ws=np.load(f'{storageStem}/NTM_ws_m.npy')

    info = plasmaInfo(targetDir, gfile,rho_pol)

    with open(f'{storageStem}/K_1Dict.pkl', 'rb') as f:
        K_1Dict = pickle.load(f)

    

    #StuartXarray = xr.open_dataset('/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/g180403_65_04400_xr',engine = 'scipy')

    if shot == '184281':
        relativeAlignments = [0,0,.25,.5,0]#[0, 0, 0.25, 0.5]
        waveTypes = ['None', 'LH', 'LH', 'LH', 'LH']# = [False, True, True,True,True]#[False, True, True, True]
        LHpowers = [0,2,2,2,1]
        plotLabels = ['No LH']
    if shot == '180403':
        relativeAlignments = [0,0,-1]#.25,.5,]
        waveTypes = ['None', 'LH', 'EC']
        LHpowers = [0,1,1]#[0,1,1,1,]
        plotLabels = ['No LH']    

        """
        relativeAlignments = [0,0,.25,.5]#.25,.5,]
        waveTypes = ['None', 'LH', 'LH', 'LH']
        LHpowers = [0,1,1,1]#[0,1,1,1,]
        plotLabels = ['No LH'] 
        """

    if shot == '199605':
        relativeAlignments = [0,0,.25,.5,0]
        waveTypes = ['None', 'LH', 'LH', 'LH', 'LH']
        LHpowers = [0,2,2,2,1]
        plotLabels = ['No LH']   

        relativeAlignments = [0,0,.25,.5,]
        waveTypes = ['None', 'LH', 'LH', 'LH']
        LHpowers = [0,2,2,2,]
        plotLabels = ['No LH'] 

    #"""
    for i in range(len(waveTypes)):
        CDresultDict = None
        

        if waveTypes[i] == 'LH':
            plotLabels.append(rf'{LHpowers[i]} MW LH, $\Delta R/\delta_{{LH}} = $' + f'{relativeAlignments[i]}')
            targetSuffix = f'_{LHpowers[i]}MW'
            with open(f'{storageStem}/{machine}_{shot}{time}/{machine}_{shot}{time}{targetSuffix}_resultDict.pkl', 'rb') as f:
                CDresultDict = pickle.load(f)

        if waveTypes[i] == 'EC':
            plotLabels.append(rf'{LHpowers[i]} MW EC, $\Delta R/\delta_{{LH}} = $' + f'{relativeAlignments[i]}')

        plotLabel = plotLabels[i]
        relAlign = relativeAlignments[i]

        generateGrowthRates(info=info, ws=ws, ax_w=ax_w,K_1Dict=K_1Dict,
                        plotLabel=plotLabel, GENCQLResultDict = CDresultDict, relAlign = relAlign, waveType=waveTypes[i], index=i)
    

    """ 
    islandWs = [0.05,.075,.1,.15,.2]
    deltaT=0.1 
    if shot == '199605':
        alignments = [0,.25,.5,0]
        powers = [2,2,2,1]
        colors = ['tab:orange','tab:green','tab:red','tab:purple']
    
    if shot == '184281':
        alignments = [0,.25,.5,0]
        powers = [2,2,2,1]
        colors = ['tab:orange','tab:green','tab:red','tab:purple']
    
    if shot == '180403':
        alignments = [0,.25,.5]#[0,.25,.5]
        powers = [1,1,1]#[1,1,1]
        colors = ['tab:orange','tab:green','tab:red']
    for p in range(len(alignments)):
        islandGrowthRates = np.zeros(len(islandWs))
        for k in range(len(islandWs)):
            islandSuffix = f'_{powers[p]}MW_{islandWs[k]}Width_{deltaT}deltaT'#islandFileSuffixes[k]
            islandw = islandWs[k]

            with open(f'{storageStem}/{machine}_{shot}{time}/{machine}_{shot}{time}{islandSuffix}_resultDict.pkl', 'rb') as f:
                CDresultDict = pickle.load(f)

            plotLabel = ''
            islandGrowthRates[k] = generateGrowthRates(info=info, ws=[islandw], ax_w=ax_w, K_1Dict=K_1Dict,
                            plotLabel=plotLabel, GENCQLResultDict = CDresultDict, relAlign = alignments[p], waveType = 'LH')
        ax_w.scatter(islandWs, islandGrowthRates,color = colors[p])
        ax_w.plot(islandWs, islandGrowthRates,linestyle = 'dashed',color = colors[p],lw=2)
    ax_w.axvline(-1,lw=2,linestyle='dashed', color = 'k', label = fr'With $\delta T_e = {deltaT}$')
    
    """
    ax_w.set_ylabel(r'$\frac{\tau_r}{r}\frac{dw}{dt}(w)$ at q=2 surface')
    ax_w.set_xlabel(r'island width w (m)')
    ax_w.axhline(0, lw = 2, linestyle = 'dotted', color = 'k')
    ax_w.text(.16+.0085, -5.25,'Stable',fontsize = 22)
    ax_w.text(.16, 2,'Unstable',fontsize = 22)
    #ax_w.text(.16, 4.8,'Unstable',fontsize = 22)
    """
    ax_w.plot([-1,-.5],[-1,-.5], label = r'1 MW LH, $\Delta R/\delta_{LH} = 0$', lw = 2, color = 'tab:orange')
    ax_w.plot([-1,-.5],[-1,-.5], label = r'1 MW LH, $\Delta R/\delta_{LH} = 0.25$', lw = 2, color = 'tab:green')
    ax_w.plot([-1,-.5],[-1,-.5], label = r'1 MW LH, $\Delta R/\delta_{LH} = 0.5$', lw = 2, color = 'tab:red')
    ax_w.plot([-1,-.5],[-1,-.5], label = r'With $\delta T_e = 0.1$', lw = 2, color = 'k', linestyle = 'dashed')
    """
    #ax_w.legend(loc='lower center', bbox_to_anchor=(0.5,.99),ncol=2,labelspacing=0.3)
    ax_w.legend(loc='best',labelspacing=0.3)
    ax_w.set_ylim([-6,6])
    ax_w.set_xlim([0,.22])
    fig_w.tight_layout()
    plt.savefig('NTM_EC.jpg', dpi=300)

    plt.show()

generateComparisonPlot()