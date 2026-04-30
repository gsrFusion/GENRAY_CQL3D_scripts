import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from matplotlib.widgets import Slider,Button
from IPython import embed
from  scipy.linalg import eigh, solve_banded
from  scipy.interpolate import interp1d, interp2d

import netCDF4
import BuildLmatrix

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import AttenMat
import CountMatrix

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
import getGfileDict
gfileDict = getGfileDict.getGfileDict()

import getInputFileDictionary
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d')

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.rc('axes', labelsize = 20)
plt.rc('figure', titlesize = 18)

#distribution function
f = np.ma.getdata(cql_nc.variables["f"])
#pitch angles mesh at which f is defined
pitchAngleMesh = np.ma.getdata(cql_nc.variables["y"][:])

c = 299792458
normalizedVel = np.ma.getdata(cqlrf_nc.variables["x"][:])
vnorm = np.ma.getdata(cqlrf_nc.variables["vnorm"][:])
#see cql3d manual for how these energies are obtained from the normalized velocity
cql3dEnergies = (6.242e15)*(-1 + np.sqrt(1 + np.square(normalizedVel*vnorm/100)/c**2))*(9.109e-31*c**2)

#cql3d rho bins
rya = np.ma.getdata(cql_nc.variables["rya"][:])

nv = cqlinput['setup']['nv']#number of sight lines
nen = cqlinput['setup']['nen']#number of energy bins in CQL3D's synthetic x-ray diagnostic
en_ = np.ma.getdata(cql_nc.variables["en_"][:])#energy bins used by CQL3D for the XR detector
emin = cqlinput['setup']['enmin']#min photon energy we are looking for
emax = cqlinput['setup']['enmax']#max photon energy we are looking for

def update_fill_between(fill,x,y_low,y_up,min,max ):
    paths, = fill.get_paths()
    nx = len(x)
    
    y_low = np.maximum(y_low, min)
    y_low[y_low==max] = min
    y_up = np.minimum(y_up,max)
    
    vertices = paths.vertices.T
    vertices[:,1:nx+1] = x,y_up
    vertices[:,nx+1] =  x[-1],y_up[-1]
    vertices[:,nx+2:-1] = x[::-1],y_low[::-1]
    vertices[:, 0] = x[0],y_up[0]
    vertices[:,-1] = x[0],y_up[0]
    
def update_errorbar(err_plot, x,y,yerr):

    plotline, caplines, barlinecols = err_plot

    # Replot the data first
    plotline.set_data(x,y)

    # Find the ending points of the errorbars
    error_positions = (x,y-yerr), (x,y+yerr)

    # Update the caplines
    if len(caplines) > 0:
        for j,pos in enumerate(error_positions):
            caplines[j].set_data(pos)

    # Update the error bars
    barlinecols[0].set_segments(list(zip(list(zip(x,y-yerr)), list(zip(x,y+yerr))))) 
        

def triband_transpose(A):
    #transpose the tridiagonal band matrix A
    AT = np.copy(A) 
    AT[0,1:],AT[2,:-1] = A[2,:-1],A[0,1:]
    return AT
            
def triband_diag_multi(A,D):
    #multiply tridiagonal band matrix a diagonal matrix
    
    AD = np.copy(A)
    AD[0,1:]*=D[:-1]
    AD[1]*=D
    AD[2,:-1]*=D[1:]
    return AD

class HXR_tomography():
    reg_level_guess = .6
    reg_level_min = .4

    def __init__(self, chordsToUse = np.arange(1,nv+1,1), DRF = False,
            attenuate = True, ignoreThermal = False, E_pMin = emin, E_pMax = emax, regParam = .55,
            particularInversion = True, showPlot = True, n = 0.8, noiseLevel = 0, 
            brightness = None, emissivityProfile = None):
        self.name = "147634_tomo"
        BuildLmatrix.main(chordsToUse,n)
        
        #load data from L matrix
        projection_matrix = np.load(f'{targetDir}/Lmat{shotNum}_{n}n.npz')
        self.L = projection_matrix['L']   #m^-1
 
        self.grid = projection_matrix['rho_center']
        self.chords = projection_matrix['chords']
        
        #rhos at which the sightlines are most tangent to B
        self.rho_tg  = projection_matrix['rho_tg']
        #Rs at which the sightlines are most tangent to B
        self.R_tg = projection_matrix['R_tg']
        #Zs at which the sightlines are most tangent to B
        self.Z_tg = projection_matrix['Z_tg']
        #distances along the sightlines at which the sightlines are most tangent to B
        self.L_tg = projection_matrix['L_tg']
        #max of cos(theta) for each sightline, where theta is the angle between the sightline and B
        self.maxParallelity = projection_matrix['maxParallelity']
        
        #slice data in camera rows
        ind = np.where(np.diff(self.R_tg) > 0.1)[0]+1
        self.cam_row_ind = [slice(i,j) for i,j in zip(np.r_[0,ind], np.r_[ind,len(self.R_tg)])]

        #the regularization parameter to be used if particularInversion = True
        self.regParam = regParam
        #if true, do the inversion for a particular regularization parameter
        self.particularInversion = particularInversion
        #if true, show the plot
        self.showPlot = showPlot
        #minimum energy photons considered for the inversion
        self.E_pMin = E_pMin
        #maximum energy photons considered for the inversion
        self.E_pMax = E_pMax
        #if true, attenuate the brightness given by cql3d according to the attenuation curves of the sapphire window and SS plate
        self.attenuate = attenuate
        #if true, the thermal counts are ignored in the inversion
        self.ignoreThermal = ignoreThermal
        #convolve the detector response functions with the brightness
        self.DRF = DRF
        #percentage of the brightness by which to varying the brightness
        self.noiseLevel = noiseLevel

        #simulated fast electron density profile
        if type(emissivityProfile) is np.ndarray:
            self.ne = emissivityProfile
        else:
            self.ne = self.getTargetNe()
        #simulated brightness
        if type(brightness) is np.ndarray:
            self.brightness = brightness
        else:
            self.brightness = self.getBrightness()
            

        grid_edges = np.linspace(0,1,len(self.ne)+1)
        self.grid_ne = (grid_edges[1:]+grid_edges[:-1])/2
        
        #guess of uncertainty
        self.err = self.brightness*0.05+self.brightness.max()*0.01  #assume 5% error
        self.scale = np.median(self.brightness) #just a normalization to avoid calculation with huge exponents
 
        self.nr = self.L.shape[1]+1
        self.nt = 1 
 

    def regul_matrix(self,bias_left=True, bias_right=True):
        #regularization band matrix, 2. order derivative, bias left or right side to zero 
        bias = .1
        D = np.ones((3,self.nr-1))
        D[1,:] *= -2
        D[0,1] = 0
        D[2,-2] = 0

        D[1,[0,-1]] = 1e-5 #just to make D invertible
        
        if bias_left:
            D[1,0] = bias

        if bias_right:
            D[1,-1] = bias
        
        return D
    
    def PRESS(self,g, prod,S,U):
        #predictive sum of squares        
        w = 1./(1.+np.exp(g)/S**2)
        ndets = len(prod)
        return np.sum((np.dot(U, (1-w)*prod)/np.einsum('ij,ij,j->i', U,U, 1-w))**2)/ndets
    
        
    def GCV(self,g, prod,S,U):
        #generalized crossvalidation        
        w = 1./(1.+np.exp(g)/S**2)
        ndets = len(prod)
        return (np.sum((((w-1)*prod))**2)+1)/ndets/(1-np.mean(w))**2
    
    
    def FindMin(self,F, x0,dx0,prod,S,U,tol=0.01):
        #stupid but robust minimum searching algorithm.

        fg = F(x0, prod, S,U)
        while abs(dx0) > tol:
            fg2 = F(x0+dx0, prod,S,U)
                                
            if fg2 < fg:
                fg = fg2
                x0 += dx0                
                continue
            else:
                dx0/=-2.
                
        return x0, np.log(fg2)
       

    def calc_tomo(self, reg_level = 0.8, nfisher = 3, eps = 1e-2, optim_regul=False):
        #calculate tomography using optimised minimum Fisher regularisation
        #Odstrcil, T., et al. "Optimized tomography methods for plasma 
        #emissivity reconstruction at the ASDEX  Upgrade tokamak.
        #Review of Scientific Instruments 87.12 (2016): 123505.
 
        #prepare regularisation operator - contains all prior information
        #biased_edges - assume zero value at the boundaries of the grid.
        D = self.regul_matrix( np.any(self.grid < 0),True)

        #weight the contribution matrix and data by the uncertainty
        T = self.L/self.err[:,None]*self.scale
        d = self.brightness/self.err

        #flat initial estimate of the weight matrix W
        W = np.ones(self.nr-1)        
  
        #iterative calculation of minimum Fisher regularisation
        for ifisher in range(nfisher):
            #multiply tridiagonal regularisation operator by a diagonal weight matrix W
            WD = triband_diag_multi(D, W**0.5)
            
            #transpose the band matrix WD
            DTW = triband_transpose(WD) 
            
            #####    solve Tikhonov regularization (optimised for speed)
            H = solve_banded((1,1),DTW,T.T, overwrite_ab=True,check_finite=False)
            
            valid = H.sum(0) != 0
            if not np.all(valid):
                print('Warning - some LOS are not linearly independent')
            
            #fast method to calculate U,S,V = svd(H.T) of rectangular matrix 
            LL = np.dot(H.T, H)
            S2,U = eigh(LL,overwrite_a=True, check_finite=True,lower=True)  
            S = np.maximum(S2,S2[-1]*1e-20)**.5 #singular values S can be negative due to numerical uncertainty 

            #projection of the data on the U base
            p = np.dot(d,U)
            
            #calculate regularisation parameter
            g0 = np.interp(reg_level, np.linspace(0,1,len(S)), 2*np.log(S))
   
            #filtering factors attenuating high frequency eigenvectors
            w = 1./(1.+np.exp(g0)/S**2)
 
            #solution
            y = np.dot(H,np.dot(U/S2,w*p))
            #final inversion of  solution, reconstruction
            y = solve_banded((1,1),WD,y, overwrite_ab=False,overwrite_b=True,check_finite=False) 

            if ifisher < nfisher-1:
                #weight matrix for the next iteration
                W = 1/np.maximum(y,eps)#**.5
                   
        #calculate the basis in the reconstruction space
        V = np.dot(H,U/S)  
        V = solve_banded((1,1),WD,V, overwrite_ab=True,overwrite_b=True,check_finite=False)
        
        #estimate optimal regularisation level for initial guess
        if optim_regul:
            g0, log_fg2 = self.FindMin(self.PRESS, g0,1,p,S,U)
            reg_level_guess = np.interp(g0,  np.log(S2),
                                             np.linspace(0,1,len(S)))
            self.reg_level_guess = max(self.reg_level_min, reg_level_guess)**2

        self.backprojection_int = np.dot(self.L,y)*self.scale        
        self.chi2 = np.sum((d -np.dot(p*w,U.T))**2)/np.size(d)
        self.solution = y #this is the inversion
        self.reg_level = reg_level

        self.solution_err = np.sqrt(np.dot(V**2,(w/S)**2)*np.maximum(self.chi2,1))

    def do_inversion(self):
        if self.particularInversion:
            self.calc_tomo(reg_level = self.regParam**.5,  nfisher = 3, eps = 1e-5)
            if self.showPlot:
                self.plotParticular()
        else:
            self.calc_tomo(reg_level = self.reg_level_guess**.5,  nfisher = 3, eps = 1e-5)
            if self.showPlot:
                self.plotVariableInversion()

    #plot a particular inversion
    #optionally also plot the reintegrated brightness
    def plotParticular(self):
        f1 = self.plotParticularReconstruction()
        plt.show()
    
        #f2 = self.plotParticularReintB()
        #f2.show()
    
    #plot the inversion for the specified regularization parameter
    def plotParticularReconstruction(self):
        fig, ax = plt.subplots(dpi = 100)
        maxSol = np.max(self.solution)

        ax.plot(self.grid,self.solution/maxSol, label = "Predicted\nemissivity", linewidth = 3, color = 'mediumblue')
        #tomo_var = ax.fill_between(self.grid, (self.solution-self.solution_err)/maxSol,
        #            (self.solution+self.solution_err)/maxSol,alpha=.5,facecolor='royalblue',
        #            edgecolor='None', label = "statistical\n uncertainty")
        #ax.plot(rya, self.ne/np.max(self.ne), label = r"$n_e(E \geq $" + f"${self.E_pMin}$ keV$)$", color = 'k',
        ax.plot(rya, self.ne/np.max(self.ne), label = r"$n_{e,\mathrm{ fast}}$", color = 'k',
            linestyle = 'dashed', linewidth = 3) 


        ax.axhline(0,c='k', linestyle = 'dotted')
                
        ax.set_ylim([-.05,1.1])
        ax.set_xlabel(r"$\rho_{{pol}}$", fontsize = 26)
        ax.set_ylabel("Predicted emissivity", fontsize = 18)
        #ax.set_ylabel("Normalized Units", fontsize = 18)
        ax.set_xticks([-1,-.75,-.5,-.25,0,.25,.5,.75,1])     
        ax.ticklabel_format(axis = 'y', scilimits = (0,0))
        
        self.plotAbsorption(ax)
        
        legend = ax.legend(fontsize = 17,  loc='best',ncol=1)#,bbox_to_anchor=(.44,1),) 

        legend.get_frame().set_alpha(None)
        legend.get_frame().set_facecolor((1, 1, 1, 0.1))  
     
        fig.set_size_inches((7,7))#8,4.75
        fig.tight_layout()
        ax.set_xlim([0,1.01])

        return fig
    
    def plotAbsorption(self,ax):  
        urfpwrl = cqlrf_nc.variables["urfpwrl"][:]
        sdpwr = cqlrf_nc.variables["sdpwr"][:]
        spsi = cqlrf_nc.variables["spsi"][:]
        delpwr = cqlrf_nc.variables["delpwr"][:]*1e-7*1e-6#convert to MW
    
        dvol = cql_nc.variables["dvol"]
        rya = np.copy(cql_nc.variables["rya"])
        powrft = cql_nc.variables["powrft"][-1]
        ionDep = powrft*dvol/1e6

        electronPower = urfpwrl * delpwr
        rhos = rya
        rhoBinEdges = (rya[1:] + rya[:-1])/2
        frontEdge = rya[0] - (rya[1] - rya[0])/2
        backEdge = rya[-1] + (rya[-1] - rya[-2])/2

        rhoBinEdges = np.concatenate(([frontEdge],rhoBinEdges,[backEdge]))

        linDep = np.zeros(len(rhos))

        ePowerTotal = 0
        for i in range(len(delpwr)):
            ePowerInLCFS = electronPower[i][spsi[i] <= rhoBinEdges[-1]]
            spsiInLCFS = spsi[i][spsi[i] <= rhoBinEdges[-1]]
            ePowerTotal += np.sum(ePowerInLCFS)
            indices = np.digitize(spsiInLCFS, rhoBinEdges, right = False)
            np.add.at(linDep, indices-1, ePowerInLCFS)

        totalDep = linDep + ionDep
        normedTotal = np.copy(totalDep/np.max(totalDep))
        #ax.plot(rya, normedTotal, lw=2, linestyle=':', label = 'RF power\ndeposition', color = 'red')
    
    #plot the reintegrated emissivity for the inversion of a specified regularization parameter
    def plotParticularReintB(self):        
        figRe, axRe = plt.subplots()
        axRe.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
        axRe.plot(self.rho_tg[self.rho_tg.argsort()], self.backprojection_int[self.rho_tg.argsort()], label = f"Reintegrated Emissivity", lw=1.5)
        axRe.plot(self.rho_tg[self.rho_tg.argsort()], self.brightness[self.rho_tg.argsort()], label = "Brightness", lw=1.5)
    
        axRe.legend(fontsize = 16, loc = 'upper left')
        axRe.set_xlim([0,1.05])
        axRe.set_ylim([-.1,max(max(self.backprojection_int), max(self.brightness))*1.05])
        axRe.set_ylabel("Brightness (counts/s)")
        axRe.set_xlabel(r'tangential $\rho_{pol}$')
        figRe.set_size_inches((7,7))
        figRe.tight_layout(rect=[0, 0.0, 1, 1])
        return figRe

    #allow the user to vary the regularization parameter and plot the resulting inversion and reintegrated brightness
    def plotVariableInversion(self):
        f,ax = plt.subplots(1,2, sharex='col', figsize=(13,7), dpi = 100, num = 100)
        ax_time = plt.axes([0.2, 0.05, 0.65, 0.03], facecolor='y')
        slide_reg = Slider(ax_time, 'Regularisation:', 0.0, 1,
                    valinit=self.reg_level_guess, valfmt='%1.3f')

        f.subplots_adjust(bottom=.2)
 
        tomo_var = ax[0].fill_between(self.grid, self.grid*0, self.grid*0,
                                      alpha=.5, facecolor='b', edgecolor='None')
        tomo_mean, = ax[0].plot([],[], lw=2)
         
        retrofits = []
        retrofits_inter = []

        #ground true
        ne = np.interp(np.abs(self.grid), self.grid_ne, self.ne) 
        ax[0].plot(self.grid, ne/np.max(ne),'--', linewidth = 2)
                
        reIntCurve, = ax[1].plot([], [], label = "Reintegrated emissivity")
        brightnessCurve = ax[1].plot(self.rho_tg[self.rho_tg.argsort()], self.brightness[self.rho_tg.argsort()], label = "brightness")
  
        ax[0].axvline(0,c='k',ls='--')
        ax[1].axvline(0,c='k',ls='--')
        ax[0].axhline(0,c='k')
        ax[1].axhline(0,c='k')
        ax[0].set_xlim(self.grid[0]-0.01, self.grid[-1]+0.01)
        ax[1].set_xlim(self.rho_tg.min()-.01, self.rho_tg.max()+.01)
        ax[0].set_ylim(-.05, 1.1)#-self.ne.max()*0.02, self.ne.max()*1.2)
        ax[1].set_ylim(0, self.brightness.max()*1.2)
        ax[0].set_xlabel(r'$\rho_{pol}$')
        ax[1].set_xlabel(r'tangential $\rho_{pol}$')
        ax[1].set_ylabel('Brightness [W/m$^2$s]')
        ax[0].set_ylabel('Density [a.u.]')
        ax[1].ticklabel_format(axis='y', scilimits=[-3, 3])

        title = f.suptitle('')

        #recalculates the inversion when a new regularization parameter is chosen
        def update(val):
            
            self.calc_tomo(reg_level = val**.5,  nfisher = 3, eps = 1e-5)
            scale = np.max(self.solution)#np.sum(self.solution*ne)/np.sum(ne**2)

            update_fill_between(tomo_var,self.grid,(self.solution-self.solution_err)/scale,
                                (self.solution+self.solution_err)/scale,-np.inf,np.inf)
            tomo_mean.set_data(self.grid,self.solution/scale)
            reIntCurve.set_data(self.rho_tg[self.rho_tg.argsort()], self.backprojection_int[self.rho_tg.argsort()])

            scaledSolution = self.solution/np.max(self.solution)
            scaledNe = ne/np.max(ne)
            MSE = (1/len(ne))*np.sum((scaledNe - scaledSolution)**2)
            RMSE = np.sqrt(MSE)
            #RMSENormed = RMSE/np.mean(ne/np.max(ne))

            avgStatErr = np.sum(self.solution_err)/(np.max(self.solution)*len(self.solution))

            title.set_text('  $\chi^2/nDoF$ = %.1f'%( self.chi2) + f"  Scaled RMSE: {RMSE: .3f}" + \
                f" <Stat Err> = {avgStatErr: .3f}" + f"   {len(self.chords)} sightlines")
            f.canvas.draw_idle()


        ax[1].legend()
        slide_reg.on_changed(update)
        update(slide_reg.valinit)
        
        axbutton = plt.axes([0.85, 0.1, 0.1, 0.05])
        self.finalize_button = Button(axbutton, 'Finalize Plot')
        self.finalize_button.on_clicked(self.finalizePlot)
        plt.show()
    #this function was originally going to try to close the variable inversion plot
    #but I cannot figure out a clean way to do that without closing piscope as well
    def finalizePlot(self, event):
        self.plotParticular()

    """    
    def save(self,args):
        np.savez_compressed('Reconstruction',backprojection= self.backprojection_int[::self.interleave],
                            reg_param =self.reg_level,  grid=self.grid,chi2=self.chi2,
                            density=self.solution ,density_err = self.solution_err )
        
        print('Saved to Reconstruction_'+self.name+'.npz file')
    """

    #returns the electron density profile of all electrons a certain energy
    #minEnergy is in keV
    def getTargetNe(self):
        ne = np.zeros(len(rya))
        indices =  []#np.where(cql3dEnergies < self.E_pMin)[0]
        if len(indices) == 0:
            indices = np.where(cql3dEnergies < 50)[0]
        minCQL3DEnergyIndex = indices[-1]
        fRelevant = f[:, minCQL3DEnergyIndex:, :]

        if self.attenuate:
            attenMat = AttenMat.getAttenFunc(self.chords, cql3dEnergies[minCQL3DEnergyIndex:])
            
            megaMat = np.zeros(fRelevant.shape)
            for i in range(len(attenMat)):
                megaMat += fRelevant * attenMat[i][None, :, None]
            megaMat /= len(self.chords)

            #fRelevant = megaMat#fRelevant# * attenuation[None, :, None]

        for i in range(0, len(rya)):
            integFOverVel = np.ma.getdata(np.trapz(fRelevant[i,:,:]*normalizedVel[minCQL3DEnergyIndex:, None]**2, 
                normalizedVel[minCQL3DEnergyIndex:, None], axis = 0))
            #print(f"{type(integFOverVel), type(pitchAngleMesh)}")
            ne[i] = 2*np.pi*np.trapz(integFOverVel*np.sin(pitchAngleMesh[i]), pitchAngleMesh[i], axis = 0)

        return ne
    
    def getBrightness(self):
        #get the counts per chord per energy bin
        _,countMatrix = CountMatrix.getCountMatrix(self.chords, attenuate = self.attenuate, 
            includeResponseFunc = self.DRF, E_pMin=self.E_pMin,E_pMax=self.E_pMax)

        print(f'countMatrix.shape: {countMatrix.shape}')
        countsPerChord = np.sum(countMatrix, axis = 2)
        print(f'max of countsPerchord: {np.max(countsPerChord)}')

        #adds the thermal and nonthermal contributions since in experiment we can't tell them apart
        countMatrixThermal = countMatrix[0]
        countMatrixNonThermal = countMatrix[1]

        #sum over energy bin to get the counts per chord 
        thermalCountsPerChord = np.sum(countMatrixThermal, axis = -1)
        nonthermalCountsPerChord = np.sum(countMatrixNonThermal, axis = -1)
        print(f"max of count matrix: {np.max(thermalCountsPerChord + nonthermalCountsPerChord)}")      
        """
        from numpy.random import poisson, uniform
        for i in range(len(countsOfInterest)):
            cleanCount = countsOfInterest[i]
            #noisyCount = poisson(cleanCount)
            if self.noiseLevel != 0:
                adjustment = uniform(self.noiseLevel-.1, self.noiseLevel,1)
                factor = 1
                negPos = uniform(-1,1, 1)
                if negPos < 0:
                    factor = factor - adjustment
                else:  
                    factor = factor + adjustment
                #print(factor)
                noisyCount = cleanCount * factor
            
                countsOfInterest[i] = noisyCount
        """
        print(f'attenuate = {self.attenuate}')
        print(f"min, max, av thermal counts/s per chord: {np.min(thermalCountsPerChord) : .3e}, {np.max(thermalCountsPerChord) : .3e}, {np.sum(thermalCountsPerChord)/len(self.chords) : .3e} between {self.E_pMin} and {self.E_pMax} keV")
        print(f"min, max, av nonthermal counts/s per chord: {np.min(nonthermalCountsPerChord) : .3e}, {np.max(nonthermalCountsPerChord) : .3e}, {np.sum(nonthermalCountsPerChord)/len(self.chords) : .3e} between {self.E_pMin} and {self.E_pMax} keV")
        
        minTimeResolutions = 1/(nonthermalCountsPerChord*.01)
        print(f'min time resolutions: {minTimeResolutions}')
        print(f'median min time resolution: {np.median(minTimeResolutions)}')
        
        print(f'nonthermal counts: {nonthermalCountsPerChord}')
        if self.ignoreThermal:
            return nonthermalCountsPerChord
        else:
            return thermalCountsPerChord + nonthermalCountsPerChord


def gaussian(mu, std, xs):
    return 1/np.sqrt(2*np.pi)*np.exp(-(xs-mu)**2/std**2/2)

def fakeHeavi(rho):
    if rho > .5:
        return 1
    else:
        return 0

def getFakeEmissFunc():
    fakeRhos = np.linspace(0,1,1000)
    gaussian1 = gaussian(0,.025, fakeRhos)
    gaussian2 = gaussian(.2,.025, fakeRhos)
    gaussian3 = gaussian(.4,.025, fakeRhos)
    gaussian4 = gaussian(.6,.025, fakeRhos)
    gaussian5 = gaussian(.8,.03,fakeRhos)
    fakeEmiss = gaussian1 + gaussian2 + gaussian3 + gaussian4 + gaussian5

    fakeEmissivityFunc = interp1d(fakeRhos, fakeEmiss)

    return fakeEmissivityFunc

def getFakeBrightness(chords, fakeEmissivityFunc):
    brightness = np.zeros(len(chords))
    #polar thetas as measured from the vertical axis
    thet1 = cqlinput['setup']['thet1']*np.pi/180.
    #toroidal thetas as measured from the x axis
    thet2 = cqlinput['setup']['thet2']*np.pi/180.
    for i in range(len(chords)):
        chordNum = chords[i]
        torTheta = thet2[chordNum-1]
        polarTheta = thet1[chordNum-1]
        losDir = np.array([np.cos(torTheta)*np.sin(polarTheta), np.sin(torTheta)*np.sin(polarTheta), np.cos(polarTheta)])

        brightness[i] = integrateFakeEmissivity(fakeEmissivityFunc, losDir)

    #fig,ax = plt.subplots()
    #ax.plot(fakeRhos, fakeEmiss)
    #plt.show()
    return brightness
#returns the rho at a given position
def getRhoFromRZ(r,z):
    rgrid = gfileDict["rgrid"]
    zgrid = gfileDict["zgrid"]

    #relevant variables to find the normalized poloidal flux
    psirz = gfileDict["psirz"]
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)
    #interpolated function for poloidal flux
    psirzNormFunc = interp2d(rgrid, zgrid[9:-10], psirzNorm[9:-10, :])

    return np.sqrt(psirzNormFunc(r,z))

def integrateFakeEmissivity(fakeEmissivityFunc, losDir):
    brightness = 0

    x_sxr = cqlinput['setup']['x_sxr']/100.  # [m]
    y_sxr = 0#by convetion of CQL3D
    z_sxr = cqlinput['setup']['z_sxr']/100.  # [m]

    #current point being looked at
    currentX = x_sxr
    currentY = y_sxr
    currentZ = z_sxr

    dx = .67/100
    dl = np.sqrt(dx**2 + (losDir[1]*(-dx/losDir[0]))**2 + (losDir[2]*(-dx/losDir[0]))**2)
    
    #while we are looking at a point inside a reasonable volume
    maxR = np.sqrt(x_sxr**2 + y_sxr**2)
    while(np.abs(currentZ) < 1.12 and 1 < np.sqrt(currentX**2 + currentY**2) <= maxR):
        currentR = np.sqrt(currentX**2 + currentY**2)
        rho = getRhoFromRZ(currentR, currentZ)[0]
        if rho <= rya[-1]:
            brightness += fakeEmissivityFunc(rho)*dl

        #get the next point to look at based on the direction of the line of sight
        newX = -dx+currentX
        newY = losDir[1]*(-dx/losDir[0]) + currentY
        newZ = losDir[2]*(-dx/losDir[0]) + currentZ
        
        currentX = newX
        currentY = newY
        currentZ = newZ
    return brightness

def main(args):
    tomo = None
    if len(args) == 0:
        #all sightlines
        _all = np.arange(1,nv+1,1)
        #four and half rows, our typical choice
        _imageBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])
        _28 = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65])
        _upperBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])+1
        _HFS = np.array([42,52,62, 43,53,63,  44,54,64, 45,55,65, 41])
        _LFS = np.array([2,12,22,32, 3,13,23,33,  4,14,24,34, 5,15,25,35, 11,21,31,])
        _highCounts = np.array([29,40,50,60,61,51,12,62,63,64,54,55,46,47,37,27,17,5,4,30,20,21,36,33,34])
        
        #fakeEmissFunc = getFakeEmissFunc()
        #fakeEmissivity = fakeEmissFunc(rya)
        #fakeBrightness = getFakeBrightness(_imageBand, fakeEmissFunc)
    
        #regParam =  .6 for 2.3, .62 for 2.7, .63 for 3.1
	# = .63 for .25 MW, .64 for .5 MW
        tomo = HXR_tomography(chordsToUse = _28, DRF = False, attenuate = True, ignoreThermal = True,
                E_pMin = 50, E_pMax = 1e3, regParam = .6, particularInversion = True,
                showPlot = True, n = .8, noiseLevel = 0)#brightness = fakeBrightness, emissivityProfile = fakeEmissivity)
    else:
        car = los
        tomo = HXR_tomography(chordsToUse = args[0], DRF = args[1], attenuate = args[2],
                E_pMin = args[3], E_pMax = args[4], regParam = args[5], particularInversion = args[6],
                showPlot = args[7], n = args[8], noiseLevel = args[9], brightness = None, emissivityProfile = None)

    tomo.calc_tomo(optim_regul=True)
    tomo.do_inversion()

    return tomo
    
main([])





