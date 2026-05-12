###
# Let me start by saying this script is a bit of an abomination, and I'm sorry for that
# That being said, I think it is possible to readily become friends with the abomination
#
# This script makes a new directory for the genray/cql runs if needed, moves the eqdsk and namelist file templates there,
# then uses the eqdsk and profile files to fill in those templates as desired
#
# To write to the namelists, it loops through the files and looks for the desired variables and then sets their values
# There is a python module that does this in a much cleaner way, but I had issues with some variables, so I wrote this instead
#

GENRAY_CQL3D_scriptsLocation = f'/home/grantr/codes/GENRAY_CQL3D_scripts'

import numpy as np
import matplotlib.pyplot as plt
import shotToEqdsk
from scipy.interpolate import interp1d
import netCDF4
import os
import helperFunctions as helperFuncs
import getGfileDict
import pickle 
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 20)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 20)

class InputFileHelper:

    """
    targetDir
    waveType -> LH, or EC
    makeDir -> if the target dir doesn't exist, make it
    overwrite -> if the target dir exist and there are already input files in it, overwrite them
    doPlot -> plot profiles to see if things are working well
    nScale, TScale , ZeffScale  -> scale factors for n, T, Zeff
    numCQLToFokkerPlanck -> number of surface to run cql at
    rya -> flux surfaces to Fokker Planck
    ndens, njene -> number of grid points when inputting the density and temperature profiles into GENRAY,CQL3D, respectively
    isScoping -> if it's just a scoping run, I reduce the resolution
    eqsym -> How CQL3D should treat the eqdsk, see cql3d helpfile
    OMFIT_nc_derived ->
    thgrill -> poloidal location of LH or FW antenna
    powerInLobes -> array (or int) of power in lobe(s)
    N_para_edges -> used for anmin, anmax in genray grill setup
    pwrScale -> cql3d's power scale, see cql3d helpfile
    N_para_peaks -> center of each LH or FW N|| lobe
    ecconeParamDict -> if waveType = 'EC', this is used to fill out the eccone namelist section. Barely used so needs to be checked
    makeIsland -> add a magnetic island 
    islandParamDict -> dictionary of parameters about the desired island
    """
    def __init__(self, targetDir,  
                 waveType = 'LH',
                 makeDir = True, overwrite = True, doPlot = True,
                 nScale = 1, TScale = 1, ZeffScale = 1,  
                 numCQLToFokkerPlanck = 50, ndens = 51, njene= 51, rya = None,
                 isScoping = False, eqsym = 'none', OMFIT_nc_derived = None,
                 thgrill=None, powerInLobes = None,  N_para_edges = None, pwrScale = 1, N_para_peaks = -2.7,
                 frqncy = 4.6e9,
                 ecconeParamDict = None, 
                 makeIsland = False, islandParamDict = None
                ):

        self.targetDir = targetDir
        targetSplit = targetDir.split('/')
        self.parentShotDir = targetSplit[6]
        self.shot = self.parentShotDir.split('_')[1]
        self.machine = self.parentShotDir.split('_')[0]
        self.makeDir = makeDir
        self.overwrite = overwrite
       
        self.nScale = nScale
        self.TScale = TScale
        self.ZeffScale = ZeffScale
        self.pwrScale = pwrScale
        
        self.numCQLToFokkerPlanck = numCQLToFokkerPlanck
        self.ndens = ndens
        self.njene = njene
        self.rya = rya

        self.isScoping = isScoping
        self.doPlot = doPlot
        
        self.eqsym = eqsym
        
        self.OMFIT_nc_derived = OMFIT_nc_derived
        self.waveType = waveType
        self.frqncy = frqncy

        self.islandParamDict = islandParamDict
        self.makeIsland = makeIsland

        self.numSpecies = -1

        if self.waveType == 'LH':
            self.thgrill = thgrill
            self.N_para_peaks = N_para_peaks
            self.N_para_edges = N_para_edges
            self.powerInLobes = powerInLobes

        if self.waveType == 'EC':
            if ecconeParamDict == None:
                raise Exception('you need to provide eccone parameters')
            
            self.zst = ecconeParamDict['zst']
            self.rst = ecconeParamDict['rst']
            self.alfast = ecconeParamDict['alfast']
            self.betast = ecconeParamDict['betast']
            self.alpha1 = ecconeParamDict['alpha1']
            self.powtot = ecconeParamDict['powtot']
            self.isX = ecconeParamDict['isX'] #X mode or O mode

            self.ioxm = -1
            if not self.isX:
                self.ioxm = 1

        if self.nScale != 1:
            print(f'ALERT: nScale IS NOT UNITY. ARE YOU SURE?')
        if self.TScale != 1:
            print(f'ALERT: TScale IS NOT UNITY. ARE YOU SURE?')
        if self.ZeffScale != 1:
            print(f'ALERT: ZeffScale IS NOT UNITY. ARE YOU SURE?')

    
    #copies the input file templates into the target directory
    #if the target directory does not exist, a new directory will be made if makeDir = true
    #if overwrite = True, any existing input files in the target directory will be overwritten
    def copyInputFileTemplates(self):
        
        #if the target dir doesn't exist, either make one or break
        if not os.path.exists(self.targetDir):
            if self.makeDir:
                print(f'making dir {self.targetDir}')
                os.system(f'mkdir -p {self.targetDir}')
            else:
                raise Exception(f"target directory {self.targetDir} does not exist") 
        
        dir_list = os.listdir(self.targetDir)
        if ('cqlinput' in dir_list or 'genray.in' in dir_list) and self.overwrite == False:
            raise Exception('input files already present. Refusing to overwrite')
        
        else:
            """
            For DIII-D cases, since we know what the predicted launched spectrum of the full launcher is for a given peak N|| value,
            there's the option to just specify the peak N|| value and it will pull the template that has the spectrum already baked in
            However, this option is somewhat deprecated
            Suggest always using generateNparaSpectrum.py to create a new spectrum each time
            This will be required for matching experiment
            """
            if self.machine == 'DIIID':
                os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/cqlinput {self.targetDir}/cqlinput')
                os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/genray_ece.in {self.targetDir}/genray_ece.in')
                if isinstance(self.N_para_peaks, float):
                    factor = -1*np.sign(self.N_para_peaks)
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/genray_{self.N_para_peaks*factor}.in {self.targetDir}/genray.in')
                else:
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/genray_{len(self.N_para_peaks)}Lobes.in {self.targetDir}/genray.in')
            
            elif self.machine == 'KSTAR':
                os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/cqlinput {self.targetDir}/cqlinput')
                os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/genray_{len(self.N_para_peaks)}Lobes.in {self.targetDir}/genray.in')

            elif self.machine == 'WEST':
                os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/cqlinput {self.targetDir}/cqlinput')
                if isinstance(self.N_para_peaks, float):
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}s/templates/{self.machine}_templates/genray_{1}Lobe.in {self.targetDir}/genray.in')
                else:
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/genray_{len(self.N_para_peaks)}Lobe.in {self.targetDir}/genray.in')
            
            elif self.machine == 'FENIX':
                os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/cqlinput {self.targetDir}/cqlinput')
                if self.waveType == 'LH':
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/genray.in {self.targetDir}/genray.in')
                if self.waveType == 'EC':
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/{self.machine}_templates/genray_EC.in {self.targetDir}/genray.in')
            
            elif self.machine == 'NTPT':
                if 'DIIID' in self.shot:
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/DIIID_templates/cqlinput {self.targetDir}/cqlinput')
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/DIIID_templates/genray_1Lobes.in {self.targetDir}/genray.in')
                elif 'ARC' in self.shot:
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/NTPT_templates/cqlinput_ARC {self.targetDir}/cqlinput')
                    os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/templates/NTPT_templates/genray_1Lobes_ARC.in {self.targetDir}/genray.in')

            #copy in the files to run GENRAY/CQL3D
            os.system(f'cp ~/codes/genr_yuri.pbs {self.targetDir}')
            os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/runGENThenCQL.sh {self.targetDir}')
            os.system(f'cp {GENRAY_CQL3D_scriptsLocation}/runECE.sh {self.targetDir}')
            #if you want to optionally have two sets of scripts for whether a run will be computationally intensive
            #if you don't want this option, just have one set of these copying commands
            if self.isScoping == False:
                os.system(f'cp ~/codes/cql.pbs {self.targetDir}/cql.pbs')
                os.system(f'cp ~/codes/genr_sam_beefy.pbs {self.targetDir}/genr_sam.pbs')
            else:
                os.system(f'cp ~/codes/genr_sam.pbs {self.targetDir}/genr_sam.pbs')
                os.system(f'cp ~/codes/cql_scoping.pbs {self.targetDir}/cql.pbs')

    #populates the empty input file templates
    """
    The basic idea is to first get the profile fit file and make a set of functions that will then be used to produce the profiles on the input grids
    Each machine has a different type of profile file, so there's a messy set of if statements for each machine type
    Additionally, not each machine has the same number of species. 
        For example DIII-D -> electrons, deuterons, carbon
        ARC -> D, T, e, impurity
        self.numSpecies takes care of this
    """
    def populateInputFiles(self):
        if self.machine == 'DIIID':
            self.numSpecies = 3
            if self.OMFIT_nc_derived == None:
                self.OMFIT_nc_derived = netCDF4.Dataset(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/FIT.nc','r')
            else:
                print(f'fit file provided')
                
            rho_psi = np.sqrt(self.OMFIT_nc_derived.variables["psi_n"])
            if len(rho_psi) ==0 or len(rho_psi) == 1:
                rho_psi = rho_psi[0]
            n_e = np.copy(self.OMFIT_nc_derived.variables["n_e"][0])
            T_e = np.copy(self.OMFIT_nc_derived.variables["T_e"][0])/1e3
            
            try:
                n_i = np.copy(self.OMFIT_nc_derived.variables['n_2H1'][0])
                n_12C6 = np.copy(self.OMFIT_nc_derived.variables["n_12C6"][0])
                T_i = np.copy(self.OMFIT_nc_derived.variables["T_12C6"][0]/1e3)
                Zeff = np.copy(self.OMFIT_nc_derived.variables['Zeff'][0])
            except:
                print(f'couldnt find everytging in .nc file. Copying electron info into ions')
                Zeff = np.ones(len(n_e))*self.ZeffScale
                
                n_i = np.copy(n_e)*(6-Zeff)/5
                n_12C6 = np.copy(n_e)*(Zeff-1)/30

                T_i = np.copy(T_e)

            neFunc = interp1d(rho_psi, n_e, kind = 'quadratic')
            niFunc = interp1d(rho_psi, n_i, kind = 'quadratic')
            n12C6Func = interp1d(rho_psi, n_12C6, kind = 'quadratic')
            TeFunc = interp1d(rho_psi, T_e, kind = 'quadratic')
            TiFunc = interp1d(rho_psi, T_i, kind = 'quadratic')

            Zeff[Zeff<1] = 1
            self.ZeffFunc = interp1d(rho_psi, Zeff)

            self.speciesLabels = ['e', 'D', '12C6']
            self.nFunctions = [neFunc, niFunc, n12C6Func]
            self.TFunctions = [TeFunc, TiFunc, TiFunc]
        elif self.machine == 'NTPT':
            if 'DIIID' in self.shot:
                self.numSpecies = 3
                if self.OMFIT_nc_derived == None:
                    self.OMFIT_nc_derived = netCDF4.Dataset(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/FIT.nc','r')
                else:
                    print(f'fit file provided')
                    
                rho_psi = np.sqrt(self.OMFIT_nc_derived.variables["psi_n"])
                if len(rho_psi) ==0 or len(rho_psi) == 1:
                    rho_psi = rho_psi[0]
                n_e = np.copy(self.OMFIT_nc_derived.variables["n_e"][0])
                T_e = np.copy(self.OMFIT_nc_derived.variables["T_e"][0])/1e3
                
                try:
                    n_i = np.copy(self.OMFIT_nc_derived.variables['n_2H1'][0])
                    n_12C6 = np.copy(self.OMFIT_nc_derived.variables["n_12C6"][0])
                    T_i = np.copy(self.OMFIT_nc_derived.variables["T_12C6"][0]/1e3)
                    Zeff = np.copy(self.OMFIT_nc_derived.variables['Zeff'][0])
                except:
                    print(f'couldnt find everything in .nc file. Copying electron info into ions')
                    Zeff = np.ones(len(n_e))*self.ZeffScale
                    
                    n_i = np.copy(n_e)*(6-Zeff)/5
                    n_12C6 = np.copy(n_e)*(Zeff-1)/30

                    T_i = np.copy(T_e)

                neFunc = interp1d(rho_psi, n_e)
                niFunc = interp1d(rho_psi, n_i)
                n12C6Func = interp1d(rho_psi, n_12C6)
                TeFunc = interp1d(rho_psi, T_e, kind = 'quadratic')
                TiFunc = interp1d(rho_psi, T_i)

                Zeff[Zeff<1] = 1
                self.ZeffFunc = interp1d(rho_psi, Zeff)

                self.speciesLabels = ['e', 'D', '12C6']
                self.nFunctions = [neFunc, niFunc, n12C6Func]
                self.TFunctions = [TeFunc, TiFunc, TiFunc]
            elif 'ARC' in self.shot:
                self.numSpecies = 4
                with open(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/ARCprofs.pkl', 'rb') as f:
                    ARC_profs = pickle.load(f)
                    
                rho_psi = np.sqrt(ARC_profs['polflux'])
                n_e = ARC_profs['ne']*1e19
                n_i = ARC_profs['ni']*1e19
                T_e = ARC_profs['te']
                T_i = ARC_profs['ti']

                n_D = n_i[:,0]
                n_T = n_i[:,1]
                n_imp = n_i[:,2]
                T_D = T_i[:,0]
                T_T = T_i[:,1]
                T_imp = T_i[:,2]                    

                neFunc = interp1d(rho_psi, n_e)
                ndFunc = interp1d(rho_psi, n_D)
                ntFunc = interp1d(rho_psi, n_T)
                nimpFunc = interp1d(rho_psi, n_imp)
                TeFunc = interp1d(rho_psi, T_e)
                TdFunc = interp1d(rho_psi, T_D)
                TtFunc = interp1d(rho_psi, T_T)
                TimpFunc = interp1d(rho_psi, T_imp)

                Zeff = (n_D + n_T + n_imp*(4.3537**2))/n_e

                Zeff[Zeff<1] = 1
                self.ZeffFunc = interp1d(rho_psi, Zeff)

                self.speciesLabels = ['e', 'D', 'T', 'imp']
                self.nFunctions = [neFunc, ndFunc, ntFunc, nimpFunc]
                self.TFunctions = [TeFunc, TdFunc, TtFunc, TimpFunc]
        
        elif self.machine == 'KSTAR':
            self.numSpecies = 3

            psi_n = []
            ne = []
            Te = []
            ni = []
            Ti = []

            with open(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/NE_fit.dat') as ne_file:
                lines = ne_file.readlines()
                for o,line in enumerate(lines):
                    if o < 2:
                        continue
                    else:
                        splitted = list(map(float, line.split()))
                        ne.append(splitted[2]*1e18)
                        psi_n.append(splitted[1])
            with open(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/TE_fit.dat') as Te_file:
                lines = Te_file.readlines()
                for o,line in enumerate(lines):
                    if o < 2:
                        continue
                    else:
                        splitted = list(map(float, line.split()))
                        Te.append(splitted[2]/1e3)
            with open(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/TI_fit.dat') as Ti_file:
                lines = Ti_file.readlines()
                for o,line in enumerate(lines):
                    if o < 2:
                        continue
                    else:
                        splitted = list(map(float, line.split()))
                        Ti.append(splitted[2]/1e3)

            Zeff = np.ones(len(ne))*2
                
            ni = np.copy(ne)*(6-Zeff)/5
            n12C6 = np.copy(ne)*(Zeff-1)/30

            rho_psi = np.sqrt(np.array(psi_n))
            print(f'roh_psi, ne, ni: {len(rho_psi), len(ne), len(ni)}')
            ne = np.array(ne)
            ni = np.array(ni)

            neFunc = interp1d(rho_psi, ne)
            niFunc = interp1d(rho_psi, ni)
            n12C6Func = interp1d(rho_psi, n12C6)
            TeFunc = interp1d(rho_psi, Te, kind = 'quadratic')
            TiFunc = interp1d(rho_psi, Ti)

            Zeff[Zeff<1] = 1
            self.ZeffFunc = interp1d(rho_psi, Zeff)

            self.speciesLabels = ['e', 'D', '12C6']
            self.nFunctions = [neFunc, niFunc, n12C6Func]
            self.TFunctions = [TeFunc, TiFunc, TiFunc]

        elif self.machine == 'WEST':
            self.numSpecies = 3
            rho_psi = np.load(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/rho_pol.npy','r')
            n_e = np.load(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/ne.npy','r')
            T_e = np.load(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/Te_keV.npy','r')
            n_12C6 = np.zeros(len(n_e))
            Exception('this probably is breaking')
            n_i = n_e
            T_i = T_e
            Zeff = np.ones(len(n_e))
            neFunc = interp1d(rho_psi, n_e, bounds_error = False, fill_value=(np.max(n_e), np.min(n_e)))
            niFunc = interp1d(rho_psi, n_i, bounds_error = False, fill_value=(np.max(n_i), np.min(n_i)))
            n12C6Func = interp1d(rho_psi, n_12C6, bounds_error = False, fill_value=(0,0))
            TeFunc = interp1d(rho_psi, T_e, bounds_error = False, fill_value=(np.max(T_e), np.min(T_e)))
            TiFunc = interp1d(rho_psi, T_i, bounds_error = False, fill_value=(np.max(T_i), np.min(T_i)))
            self.ZeffFunc = interp1d(rho_psi, Zeff, bounds_error = False, fill_value=(1,1))

            self.speciesLabels = ['e', 'D', '12C6']
            self.nFunctions = [neFunc, niFunc, n12C6Func]
            self.TFunctions = [TeFunc, TiFunc, TiFunc]

        elif self.machine == 'FENIX':
            self.numSpecies = 2
            rho_psi = np.load(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/rho_pol.npy','r')
            n_e = np.load(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/ne_FENIX.npy','r')
            T_e = np.load(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/Te_FENIX.npy','r')

            n_i = n_e
            T_i = T_e
            Zeff = np.ones(len(n_e))*1


            neFunc = interp1d(rho_psi, n_e, bounds_error = False, fill_value=(np.max(n_e), np.min(n_e)))
            niFunc = interp1d(rho_psi, n_i, bounds_error = False, fill_value=(np.max(n_i), np.min(n_i)))
            TeFunc = interp1d(rho_psi, T_e, bounds_error = False, fill_value=(np.max(T_e), np.min(T_e)))
            TiFunc = interp1d(rho_psi, T_i, bounds_error = False, fill_value=(np.max(T_i), np.min(T_i)))
            self.ZeffFunc = interp1d(rho_psi, Zeff, bounds_error = False, fill_value=(1,1))

            self.speciesLabels = ['e', 'D']
            self.nFunctions = [neFunc, niFunc]
            self.TFunctions = [TeFunc, TiFunc]

        #if there's an island, modify the temperature profile
        if self.makeIsland:
            width = self.islandParamDict['width'] #m
            deltaT = self.islandParamDict['deltaT']
            islandq = self.islandParamDict['islandq']
            Te_withIsland = np.copy(T_e)
            gfileDict = getGfileDict.getGfileDict(targetDir = self.targetDir)

            qpsi = gfileDict['qpsi']#safety factor for equally spaced psi points
            psi = np.linspace(0,1,len(qpsi))
            rho_p_q2 = np.sqrt(psi)[helperFuncs.findNearestIndex(islandq, qpsi)]

            R_LFS = helperFuncs.convertRhopolToRmidplane(rho_psi, self.targetDir, side = 'LFS')
            targetSurface_R = R_LFS[helperFuncs.findNearestIndex(rho_p_q2, rho_psi)]
            coord = R_LFS -targetSurface_R
            #go from rho_pol to R_LFS
            #center island at q=2 surface

            mask1 = np.where((-width/2 <= coord)*(coord < 0))
            mask2 = np.where((width/2 >= coord)*(coord >= 0))
           
            Te_withIsland[mask1] = T_e[mask1]*(1 + deltaT*(1-(coord[mask1]/(width/2))**2))*(1+deltaT*np.sin(np.pi*coord[mask1]/(width/2)))
            Te_withIsland[mask2] = T_e[mask2]*(1 + deltaT*(1-(coord[mask2]/(width/2))**2))
            
            Te_5 = np.copy(T_e)
            Te_10 = np.copy(T_e)
            Te_15 = np.copy(T_e)

            width = 0.05
            mask1 = ((-width/2 <= coord)&(coord < 0))
            mask2 = ((width/2 >= coord)&(coord >= 0))
           
            mask3 = mask1 | mask2
            Te_5[mask1] = T_e[mask1]*(1 + deltaT*(1-(coord[mask1]/(width/2))**2))*(1+deltaT*np.sin(np.pi*coord[mask1]/(width/2)))
            Te_5[mask2] = T_e[mask2]*(1 + deltaT*(1-(coord[mask2]/(width/2))**2))
            
            Te_5 = Te_5[mask3]
            rho_psi_5 = np.copy(rho_psi)[mask3]


            width = .15
            mask1 = ((-width/2 <= coord)&(coord < 0))
            mask2 = ((width/2 >= coord)&(coord >= 0))
           
            Te_15[mask1] = T_e[mask1]*(1 + deltaT*(1-(coord[mask1]/(width/2))**2))*(1+deltaT*np.sin(np.pi*coord[mask1]/(width/2)))
            Te_15[mask2] = T_e[mask2]*(1 + deltaT*(1-(coord[mask2]/(width/2))**2))

            mask3 = mask1 | mask2
            Te_15 = Te_15[mask3]
            rho_psi_15 = np.copy(rho_psi)[mask3]

            self.TeFunc = interp1d(rho_psi, Te_withIsland, bounds_error = False, fill_value=(np.max(T_e), np.min(T_e)))
            #"""
            if self.doPlot:
                fig,ax = plt.subplots(figsize = (6.4, 4.2))
                ax.axvline(np.sqrt(psi)[helperFuncs.findNearestIndex(islandq, qpsi)], color = 'k', lw = 2, linestyle = 'dashdot', label = r'$q = 2$')
                ax.plot(rho_psi_5, Te_5, linestyle = 'dashed', lw = 3, color = 'crimson', label = r'$\delta T_e = 0.1$' '\n'+ r'$w = 5$ cm')
                ax.plot(rho_psi, T_e, lw = 3, color = 'forestgreen', label = r'$\delta T_e = 0$', zorder = -1)
                
                #ax.plot(rho_psi, Te_10, linestyle = 'dashed', lw = 3, color = 'blue', label = r'$\delta T_e = 0.1, w = 10$ cm')
                ax.plot(rho_psi_15, Te_15, linestyle = 'dashed', lw = 3, color = 'orange', label = r'$\delta T_e = 0.1$' '\n'+ r'$w = 15$ cm')
                #ax.plot(rho_psi, Te_withIsland, linestyle = 'dashed', lw = 3, color = 'crimson', label = r'$\delta T_e = 0.1$')
                
                
                ax.set_ylabel('Electron temperature (keV)', fontsize = 18)
                ax.set_xlabel(r'$\rho_{pol}$', fontsize = 18)
                ax.set_xlim([.55,.9])
                ax.set_ylim([1.3,2.5])
                ax.legend(fontsize = 14, ncol=2,framealpha = 1,)
                fig.tight_layout()
                plt.show()
            #"""
            T_e = Te_withIsland

        if self.doPlot:
            #plot profiles
            fig, ax = plt.subplots()
            ax2 = ax.twinx()

            for p in range(self.numSpecies):
                ax.plot(rho_psi, self.nFunctions[p](rho_psi), label = self.speciesLabels[p])
                ax.scatter(rho_psi, n_e)#, label = self.speciesLabels[p])
                ax2.plot(rho_psi, self.TFunctions[p](rho_psi), label = self.speciesLabels[p], linestyle = 'dashed')

            ax.legend(loc = 'best')
            ax.set_xlim([0,1])

            #"""
            plt.show()
            #plot zeff
            fig, ax = plt.subplots()
            ax.plot(rho_psi, Zeff, color = 'k')
            ax.set_xlim([0,1])
            ax.set_ylabel('Zeff')
            #"""
            plt.show()

        self.writeGENRAY()
        if self.machine == 'DIIID':
            self.writeGENRAY_ECE()
        self.writeCQL()

    #write all the relevant variables for the cql3d namelist
    #eqdsk name is used in the deprecated electric field portion
    #for DIII-D shots, N_para is a float, the peak of the forward spectrum
    #for WEST scoping shots, N_para is a tuple, each a peak of the lobe
    def writeCQL(self):
        cqlinput = open(f'{self.targetDir}/cqlinput','r+')
        cqlinput.seek(0)

        if self.rya is None:
            self.rya = np.round(np.linspace(0.01,0.99,self.numCQLToFokkerPlanck),4)
        else:
            print(f'len(self.rya): {len(self.rya)}, self.numCQLToFokkerPlanck: {self.numCQLToFokkerPlanck}')
            assert len(self.rya) == self.numCQLToFokkerPlanck

        ryain = np.linspace(0,1,self.njene)

        #convert to 1/cm^3
        cql3dDensities = [None]*self.numSpecies
        cql3dTemperatures = [None]*self.numSpecies
        for i in range(self.numSpecies):
            cql3dDensities[i] = self.nFunctions[i](ryain)/1e6#convert to 1/cm^3
            cql3dTemperatures[i] = self.TFunctions[i](ryain)

        zeffin = self.ZeffFunc(ryain)

        E_cm = np.zeros(len(ryain))

        #convert these variables to string in preparation for adding them to the input files
        ryain_str = str((np.round(ryain,6)).tolist())[1:-1].replace(',','')
        zeffin_str = str((np.round(zeffin,6)).tolist())[1:-1].replace(',','')
        E_cm_str = str((np.round(E_cm,6)).tolist())[1:-1].replace(',','')

        ## write to cqlinput:
        cqlinput.seek(0)
        cqlData = cqlinput.readlines()

        #WEST and DIII-D have Ip in opposite directions for co-current drive
        #since I have bsign as a function of the sign of Npara, we need this factor
        if self.machine == 'WEST' or self.machine == 'FENIX':
            bsign = 1
        if self.machine == 'NTPT':
            Npara_sign = np.sign(self.N_para_peaks[0])
            if Npara_sign == -1:
                bsign = 1
            elif Npara_sign == 1:
                bsign = -1
        if self.machine == 'MANTA':
            bsign = 1
        if self.machine == 'DIIID':
            if isinstance(self.N_para_peaks, float):
                Npara_sign = np.sign(self.N_para_peaks)
            else:
                Npara_sign = np.sign(self.N_para_peaks[0])
            bsign = Npara_sign*-1
        if self.machine == 'KSTAR':#needs to be checked
            if isinstance(self.N_para_peaks, float):
                Npara_sign = np.sign(self.N_para_peaks)
            else:
                Npara_sign = np.sign(self.N_para_peaks[0])
            bsign = Npara_sign


        for i in range(len(cqlData)):
            if 'eqsym' in cqlData[i]:
                cqlData[i] = f" eqsym = '{self.eqsym}'\n"
            if 'eqdskin' in cqlData[i]:
                print(f'eqdskName:{self.eqdskName}')
                cqlData[i] = f" eqdskin =  '{self.eqdskName}'\n"
            if 'lrz =' in cqlData[i]:
                cqlData[i] = f' lrz = {self.numCQLToFokkerPlanck}\n'
            if 'rya(1)' in cqlData[i]:
                cqlData[i] = f"rya(1) = {str(self.rya.tolist())[1:-1].replace(',','')}\n"

            if 'enein(1' in cqlData[i]:
                speciesNumber = int(cqlData[i].split('=')[0][-3])
                if speciesNumber == self.numSpecies+1:
                    speciesIndex = 0
                else:
                    speciesIndex = speciesNumber - 1
                enein_str = str((np.round(cql3dDensities[speciesIndex],6)).tolist())[1:-1].replace(',','')
                cqlData[i] = f' enein(1,{speciesNumber}) = {enein_str}\n'

            if 'ryain' in cqlData[i]:
                cqlData[i] = f' ryain = {ryain_str}\n'
            if 'tein' in cqlData[i]:
                tein_str = str((np.round(cql3dTemperatures[0],6)).tolist())[1:-1].replace(',','')
                cqlData[i] = f' tein = {tein_str}\n'
            if 'tiin' in cqlData[i]:
                tiin_str = str((np.round(cql3dTemperatures[1],6)).tolist())[1:-1].replace(',','')
                cqlData[i] = f' tiin = {tiin_str}\n'
            if 'elecin' in cqlData[i]:
                cqlData[i] = f' elecin = {E_cm_str}\n'
            if 'zeffin' in cqlData[i]:
                cqlData[i] = f' zeffin = {zeffin_str}\n'
            if 'enescal' in cqlData[i]:
                cqlData[i] = f' enescal = {self.nScale}\n'
            if 'tescal' in cqlData[i]:
                cqlData[i] = f' tescal = {self.TScale}\n'
            if 'tiscal' in cqlData[i]:
                cqlData[i] = f' tiscal = {self.TScale}\n'
            if 'pwrscale(1)' in cqlData[i]:
                cqlData[i] = f' pwrscale(1) = {self.pwrScale}\n'
            if 'njene' in cqlData[i]:
                cqlData[i] = f' njene = {self.njene}\n'
            if 'bsign' in cqlData[i]:
                cqlData[i] = f' bsign = {bsign}\n'

            #if it's scoping, reduce the resolution
            if self.isScoping and self.machine != 'FENIX':
                if 'jx' in cqlData[i]:
                    cqlData[i] = f' jx =  750\n'
                if 'lz' in cqlData[i]:
                    cqlData[i] = f' lz =  50\n'
                if 'iy' in cqlData[i]:
                    cqlData[i] = f' iy =  120\n'
                if 'enorm' in cqlData[i] and not('kenorm' in cqlData[i]):
                    cqlData[i] = f' enorm =  1500\n'
                if 'nstop' in cqlData[i]:
                    cqlData[i] = f' nstop =  15\n'
                if 'nplot' in cqlData[i]:
                    cqlData[i] = f' nplot =  15\n'
            
        cqlinput.seek(0)
        cqlinput.writelines(cqlData)
        cqlinput.seek(0)
        cqlinput.writelines(cqlData)

    def writeGENRAY_ECE(self): 
        genray_in = open(f'{self.targetDir}/genray_ece.in','r+')
        genray_in.seek(0)
 
        ## write to genray:
        genray_in.seek(0)
        genrayData = genray_in.readlines()
        i = 0

        R_wall = self.gfileDict['xlim']
        Z_wall = self.gfileDict['ylim']
        n_wall = len(R_wall)
        assert len(R_wall) == len(Z_wall)

        while i < len(genrayData):
            if 'n_wall' in genrayData[i]:
                genrayData[i]= f" n_wall =  {n_wall}\n"
            if 'r_wall' in genrayData[i]:
                genrayData[i]= f' r_wall =  {" ".join(map(str, R_wall))}\n'
            if 'z_wall' in genrayData[i]:
                genrayData[i]= f' z_wall =  {" ".join(map(str, Z_wall))}\n'

            if 'ndens' in genrayData[i]:
                genrayData[i]= f" ndens =  {self.ndens}\n"

            if 'eqdskin' in genrayData[i]:
                genrayData[i] = f' eqdskin= "{self.eqdskName}"\n'
            
            if 'dentab\n' in genrayData[i]:
                for j in range(self.ndens):
                    dentabSingleLine = f''

                    for p in range(self.numSpecies):
                        dentabSingleLine += f' {self.dentab[self.numSpecies*j + p]}'

                    if j == 0:
                        genrayData[i+j+1] = f" prof ={dentabSingleLine}\n"
                    else:
                        if '&' in genrayData[i+j+1]:
                            genrayData = np.insert(genrayData, i+j+1, f"{dentabSingleLine}\n")
                        else:    
                            genrayData[i+j+1] = f"{dentabSingleLine}\n"
                i += self.ndens
            if 'temtab\n' in genrayData[i]:
                for j in range(self.ndens):

                    temtabSingleLine = f''

                    for p in range(self.numSpecies):
                        temtabSingleLine += f' {self.temtab[self.numSpecies*j + p]}'

                    if j == 0:
                        genrayData[i+j+1] = f" prof = {temtabSingleLine}\n"
                    else:
                        if '&' in genrayData[i+j+1]:
                            genrayData = np.insert(genrayData, i+j+1, f"{temtabSingleLine}\n")
                        else:    
                            genrayData[i+j+1] = f"{temtabSingleLine}\n"
                i += self.ndens
            #"""
            if 'zeff1' in genrayData[i]:
                for j in range(len(self.zeff1)):
                    if j == 0:
                        genrayData[i+j] = f" zeff1 = {self.zeff1[j]}\n"
                    if j < len(self.zeff1) - 1:
                        if '&' in genrayData[i+j+1]:
                            genrayData = np.insert(genrayData, i+j+1, f"{self.zeff1[j+1]}\n")
                        else:    
                            genrayData = np.insert(genrayData, i+j+1, f"{self.zeff1[j+1]}\n")#genrayData[i+j+1] = f"{zeff1[j]}\n"
                i += len(self.zeff1)

            i+=1
            
        #cqlinput.truncate(0)
        genray_in.seek(0)
        genray_in.writelines(genrayData)
        genray_in.seek(0)
        genray_in.writelines(genrayData) 

    #populate genray input file
    #for DIII-D shots, N_para is a float, the peak of the forward spectrum
    #for WEST scoping shots, N_para is a tuple, each a peak of the lobe
    def writeGENRAY(self):

        genray_in = open(f'{self.targetDir}/genray.in','r+')
        genray_in.seek(0)

        rho_gen = np.zeros(self.ndens)
        for i in range(len(rho_gen)):
            rho_gen[i] = (i)/(self.ndens-1)
            
        zeff_gen_values = np.round(self.ZeffFunc(rho_gen),6)
        zeff1 = [None]*int(np.ceil(len(zeff_gen_values)/5))
        
        for i in range(len(zeff1)):
            string = ''
            if i > 0:
                string = '    '
            string += (str(zeff_gen_values[i*5:i*5+5].tolist())[1:-1].replace(',',''))
            zeff1[i] = string 

        self.zeff1 = zeff1

        n_species = [None]*self.numSpecies
        T_species = [None]*self.numSpecies

        for i in range(self.numSpecies):
            n_species[i] = self.nFunctions[i](rho_gen)
            T_species[i] = self.TFunctions[i](rho_gen)

        #density/temperature for electrons, deuterium, and carbon
        self.dentab = np.zeros(self.numSpecies*self.ndens)
        self.temtab = np.zeros(self.numSpecies*self.ndens)

        for i in range(self.ndens):
            for l in range(self.numSpecies):
                self.dentab[i*self.numSpecies + l] = n_species[l][i]
                self.temtab[i*self.numSpecies + l] = T_species[l][i]
            
        R_wall = self.gfileDict['xlim']
        Z_wall = self.gfileDict['ylim']
        n_wall = len(R_wall)
        assert len(R_wall) == len(Z_wall)

        if self.machine == 'NTPT':
            n_wall = 5
            new_R = []
            new_Z = []
            new_R.append(np.min(R_wall))
            new_R.append(np.max(R_wall))
            new_R.append(np.max(R_wall))
            new_R.append(np.min(R_wall))
            new_R.append(np.min(R_wall))
            R_wall = new_R

            new_Z.append(np.max(Z_wall))
            new_Z.append(np.max(Z_wall))
            new_Z.append(np.min(Z_wall))
            new_Z.append(np.min(Z_wall))
            new_Z.append(np.max(Z_wall))
            Z_wall = new_Z

        ## write to genray:
        genray_in.seek(0)
        genrayData = genray_in.readlines()
        i = 0
        while i < len(genrayData):
            if 'n_wall' in genrayData[i]:
                genrayData[i]= f" n_wall =  {n_wall}\n"
            if 'r_wall' in genrayData[i]:
                genrayData[i]= f' r_wall =  {" ".join(map(str, R_wall))}\n'
            if 'z_wall' in genrayData[i]:
                genrayData[i]= f' z_wall =  {" ".join(map(str, Z_wall))}\n'

            if 'ndens' in genrayData[i]:
                genrayData[i]= f" ndens =  {self.ndens}\n"
            #dentab and temtab get 3 entries per line to keep it somewhat legible
            if 'dentab\n' in genrayData[i]:
                for j in range(self.ndens):
                    dentabSingleLine = f''

                    for p in range(self.numSpecies):
                        dentabSingleLine += f' {self.dentab[self.numSpecies*j + p]}'

                    if j == 0:
                        genrayData[i+j+1] = f" prof ={dentabSingleLine}\n"
                    else:
                        if '&' in genrayData[i+j+1]:
                            genrayData = np.insert(genrayData, i+j+1, f"{dentabSingleLine}\n")
                        else:    
                            genrayData[i+j+1] = f"{dentabSingleLine}\n"
                i += self.ndens
            if 'temtab\n' in genrayData[i]:
                for j in range(self.ndens):

                    temtabSingleLine = f''

                    for p in range(self.numSpecies):
                        temtabSingleLine += f' {self.temtab[self.numSpecies*j + p]}'

                    if j == 0:
                        genrayData[i+j+1] = f" prof = {temtabSingleLine}\n"
                    else:
                        if '&' in genrayData[i+j+1]:
                            genrayData = np.insert(genrayData, i+j+1, f"{temtabSingleLine}\n")
                        else:    
                            genrayData[i+j+1] = f"{temtabSingleLine}\n"
                i += self.ndens
            #"""
            if 'zeff1' in genrayData[i]:
                for j in range(len(self.zeff1)):
                    if j == 0:
                        genrayData[i+j] = f" zeff1 = {self.zeff1[j]}\n"
                    if j < len(self.zeff1) - 1:
                        if '&' in genrayData[i+j+1]:
                            genrayData = np.insert(genrayData, i+j+1, f"{self.zeff1[j+1]}\n")
                        else:    
                            genrayData = np.insert(genrayData, i+j+1, f"{self.zeff1[j+1]}\n")#genrayData[i+j+1] = f"{zeff1[j]}\n"
                i += len(self.zeff1)
                
            if 'temp_scale' in genrayData[i]:
                variableName = genrayData[i].split('=')[0]
                genrayData[i] = f'{variableName}= {self.TScale}\n'
            if 'den_scale' in genrayData[i]:
                variableName = genrayData[i].split('=')[0]
                genrayData[i] = f'{variableName}= {self.nScale}\n'
            if 'eqdskin' in genrayData[i]:
                genrayData[i] = f' eqdskin= "{self.eqdskName}"\n'
            if 'frqncy' in genrayData[i]:
                genrayData[i] = f' frqncy= {self.frqncy}\n'

            if self.waveType == 'LH':
                #for normal Bt DIIID shots
                if isinstance(self.N_para_peaks, float) and np.sign(self.N_para_peaks) > 0 and (self.machine == 'DIIID' or self.machine == 'NTPT'):
                    if 'anmax' in genrayData[i]:
                        splitInfo = genrayData[i].split('=')
                        varName = splitInfo[0]
                        maxNpara = float(splitInfo[1].strip())
                        newMax = np.round(-1*maxNpara +0.4,decimals = 3)
                        genrayData[i] = f'{varName}= {newMax}\n'
                    if 'anmin' in genrayData[i]:
                        splitInfo = genrayData[i].split('=')
                        varName = splitInfo[0]
                        minNpara = float(splitInfo[1].strip())
                        newMin = np.round(-1*minNpara -0.4,decimals = 3)
                        genrayData[i] = f'{varName}= {newMin}\n'

                if self.machine == 'WEST' or self.machine == 'FENIX' or self.machine == 'MANTA':
                    if isinstance(self.N_para_peaks, float):
                        self.N_para_peaks = np.array([self.N_para_peaks])

                if not (self.N_para_edges is None):
                    if 'anmax' in genrayData[i]:
                        splitInfo = genrayData[i].split('=')
                        varName = splitInfo[0]
                        index = int(varName.split('(')[1][0])
                        genrayData[i] = f'{varName}= {self.N_para_edges[index-1,1]:.4f}\n'
                    if 'anmin' in genrayData[i]:
                        splitInfo = genrayData[i].split('=')
                        varName = splitInfo[0]
                        index = int(varName.split('(')[1][0])
                        genrayData[i] = f'{varName}= {self.N_para_edges[index-1,0]:.4f}\n'
                        
                if 'powers' in genrayData[i] and self.powerInLobes is not None:
                    splitInfo = genrayData[i].split('=')
                    varName = splitInfo[0]
                    index = int(varName.split('(')[1][0])
                    genrayData[i] = f'{varName}= {self.powerInLobes[index-1]:.1f}\n'

                if 'thgrill' in genrayData[i] and self.thgrill is not None:
                    splitInfo = genrayData[i].split('=')
                    varName = splitInfo[0]
                    if isinstance(self.thgrill, tuple):
                        index = int(varName.split('(')[1][0])
                        genrayData[i] = f'{varName}= {self.thgrill[index-1]}\n'
                    else:
                        genrayData[i] = f'{varName}= {self.thgrill}\n'

            elif self.waveType == 'EC':
                if 'zst' in genrayData[i]:
                    genrayData[i] = f' zst = {self.zst}\n'
                if 'rst' in genrayData[i]:
                    genrayData[i] = f' rst = {self.rst}\n'
                if 'alfast' in genrayData[i]:
                    genrayData[i] = f' alfast = {self.alfast}\n'
                if 'betast' in genrayData[i]:
                    genrayData[i] = f' betast = {self.betast}\n'
                if 'alpha1' in genrayData[i]:
                    genrayData[i] = f' alpha1 = {self.alpha1}\n'
                if 'powtot' in genrayData[i]:
                    genrayData[i] = f' powtot = {self.powtot}\n'
                if 'ioxm' in genrayData[i]:
                    genrayData[i] = f' ioxm = {self.ioxm}\n'


            #reduce resolution for scoping
            if self.isScoping:
                if self.waveType == 'LH':
                    if 'nthin'in genrayData[i]:
                        variableName = genrayData[i].split('=')[0]
                        genrayData[i] = f'{variableName}= 4\n'
                    if 'nnkpar' in genrayData[i]:
                        variableName = genrayData[i].split('=')[0]
                        genrayData[i] = f'{variableName}= 25\n'
                #
                if 'nrelt' in genrayData[i]:
                    genrayData[i] = f' nrelt =  6000\n'
                if ' prmt6' in genrayData[i].split(' = '):
                    genrayData[i] = f' prmt6 =  0.005\n'
                if 'prmt4' in genrayData[i]:
                    genrayData[i] = f' prmt4 =  1e-07\n'
                #"""

            i += 1
            
        #cqlinput.truncate(0)
        genray_in.seek(0)
        genray_in.writelines(genrayData)
        genray_in.seek(0)
        genray_in.writelines(genrayData)  

    #clean the eqdsk sitting in the parent shot dir and copy it to the target directory
    #an older version of EFIT in OMFIT produced a header with weird characters, so this cleaning fixes that if present
    def cleanAndCopyEQDSK(self):
        shotNumber = self.shot.split('.')[0]

        self.eqdskName = shotToEqdsk.getEqdskName(self.shot, self.machine)
        if self.machine == 'DIIID':
            eqdskTime = self.eqdskName.split('.')[1][1:]
            eqdskFile = open(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/{self.eqdskName}','r+')
            eqdskData = eqdskFile.readlines()

            if 'b\'' in eqdskData[0]: #some DIII-D equilibria had messed up headers
                original = eqdskData[0]
                segments = original.split()
                newData = f'  {segments[1][:-2]}     xx/yy/zzzz    #{shotNumber}  {eqdskTime}ms           {segments[-3]}  {segments[-2]}  {segments[-1]}\n'
                eqdskData[0] = newData
                print(newData)
                eqdskFile.seek(0)
                eqdskFile.writelines(eqdskData)

        os.system(f'cp /home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}/{self.eqdskName} {self.targetDir}/{self.eqdskName}')
        self.gfileDict = getGfileDict.getGfileDict(f'/home/grantr/symlinks/genray_batch/{self.machine}_shots/{self.parentShotDir}')

    def copySetupAndClean(self):

        self.copyInputFileTemplates()
        print('input files copied')
        self.cleanAndCopyEQDSK()
        print('eqdsk copied')
        print(f'thgrill: {self.thgrill}, {self.thgrill == np.inf}')
        if self.thgrill == np.inf:
            print('calculating what thgrill should be')
            self.thgrill = helperFuncs.getThgrill(self.targetDir)

        self.populateInputFiles()
        print('input files populated')