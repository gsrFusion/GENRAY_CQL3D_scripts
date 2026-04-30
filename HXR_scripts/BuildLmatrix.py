#from map_equ_eqdsk import equ_map
import numpy as np
from IPython import embed
import matplotlib.pylab as plt 
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
from time import time as T
from omfit_classes.omfit_eqdsk import OMFITeqdsk #slow import
from map_equ_eqdsk import equ_map

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import getInputFileDictionary
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d')
import shotToEqdsk

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

#Load eqdsk file

coord_out = 'rho_pol'
eqdskName = shotToEqdsk.getEqdskName(shotNum,machine)
eq_time = float(eqdskName.split('.')[1])/1000 #s
print(f'eqdskName: {eqdskName}, time: {eq_time}')
eqdsk = OMFITeqdsk(f'{targetDir}/{eqdskName}')
#prepare EQ mapping object

eqm = equ_map({eq_time*1e3:eqdsk, (eq_time*1e3+1):eqdsk})
eqm.Open()
 

##### Parameters
#number of radial grid binds
nr = 200
#number of LOS
nv = 66
#relativistic beaming correction
n = -1#.8
#discretisation of each LOS
los_len = 10001#1001
#number of virtual LOS accounting for finite divergence of observation cones
nvirt = 1#200
#divergence of observation cones
divergence = 0#*.02275

#radial grid
rho_grid_edges = np.linspace(0,1, nr+1)
#rho_grid_edges = np.linspace(0,1, nr+1)

rho_grid_centers = (rho_grid_edges[1:]+rho_grid_edges[:-1])/2
  
#Diagnostic geometry
#location of XR detector in CQL3D coordinates
x_hxr = cqlinput['setup']['x_sxr']/100.  # [m]
y_hxr = 0#by convetion of CQL3D
z_hxr = cqlinput['setup']['z_sxr']/100.  # [m]
HXR_xyz = (x_hxr, y_hxr, z_hxr)


def main(chords, _n, saveFilename = None):
    chords = chords - 1

    global n;n = _n

    if saveFilename == None:
        saveFilename = f'{targetDir}/Lmat{shotNum}_{n}n'

    #polar thetas of sightlines as measured from the z axis
    thet1s = cqlinput['setup']['thet1']*np.pi/180.
    #toroidal thetas of sightlines as measured from the x axis
    thet2s = cqlinput['setup']['thet2']*np.pi/180.
    #number of LOS
    nv = int(cqlinput['setup']['nv'])
    

    #splitting chors in the dfferent diagonals
    ch_diags = [np.r_[47,57],np.r_[28:nv:10],np.r_[9:nv:10]]\
            +[np.r_[i:nv:10] for i in range(7)]\
            +[np.r_[7:38:10],np.r_[8,18]]

    #################  prepare L matrix #####################


    #unit vector in the LOS direction
    losDir = np.array([np.cos(thet2s)*np.sin(thet1s), np.sin(thet2s)*np.sin(thet1s), np.cos(thet1s)])

    ###add virtual LOS with some cone divergence
    #n1 n2 are units vectors perpendicuar to losDir
    n1  = np.array([losDir[0]*0, -losDir[2],losDir[1]])
    n1 /= np.linalg.norm(n1,axis=0)

    n2 = np.array([losDir[1]**2+losDir[2]**2, -losDir[0]*losDir[1],-losDir[0]*losDir[2]])
    n2 /= np.linalg.norm(n2,axis=0)

    #add nvirt LOS around the center of LOS, very rought model, but better than nothing
    alpha = np.linspace(0,2*np.pi,nvirt-1, endpoint=False)[:,None,None]

    cone_losDir = np.sin(divergence)*(np.cos(alpha)*n1[None]+np.sin(alpha)*n2[None])

    #add them to losDir and the central LOS will be also just a virtual LOS
    losDir = np.tile(losDir,(nvirt,1,1))
    losDir[1:] += cone_losDir #tilt the LOS in perpendicular direction
    losDir[1:] *= np.cos(divergence) #make losDir a unit vector again



    ##calculate R tangential for each LOS, L_tg is distance to tangential point nearest tokamak axis
    L_tg = - np.dot(HXR_xyz[:2],losDir[:,:2])/(losDir[:,:2]**2).sum(1)
    R_tg = np.hypot(x_hxr+losDir[:,0]*L_tg, y_hxr+losDir[:,1]*L_tg)
    Z_tg = y_hxr+losDir[:,2]*L_tg
    rho_tg = eqm.rz2rho(R_tg, Z_tg, eq_time, coord_out = coord_out)[0]
    rho_tg = rho_tg.reshape(R_tg.shape).mean(0)

    #calculate just for the center of LOS, t = 0 is tangential point 
    t = np.linspace(0,2,los_len)
    #LOS_xyz has size (nvirt, ndim=3, nv,  LOS_xyz)
    LOS_xyz = np.array(HXR_xyz,ndmin=3).T[None]+(losDir*L_tg[:,None])[...,None]*t
    #R,Z along LOSs
    LOS_R =  np.hypot(LOS_xyz[:,0],LOS_xyz[:,1])
    LOS_Z = LOS_xyz[:,2]
    #toroidal angle along LOSs
    Phi = np.arctan2(LOS_xyz[0,1], LOS_xyz[0,0], ) 


    #lenght of each step along LOS
    dL = (t[-1]-t[0])/los_len*L_tg


    # normalised toroidal flux coordinate rho_tor alogn LOSs
    LOS_rho = eqm.rz2rho(LOS_R, LOS_Z, eq_time, coord_out = coord_out)[0]
    LOS_rho = LOS_rho.reshape(LOS_R.shape)

    Raxis = np.interp(eq_time, eqm.t_eq, eqm.ssq['Rmag'])
    #HFS will have negative sign of rho
    if any(rho_grid_edges < 0):
        LOS_rho[LOS_R < Raxis] *= -1 
        rho_tg[R_tg.mean(0) < Raxis] *= -1


    ##calculate angle between LOS and magnetic field along LOS
    #Br, Bt, Bz along LOSs, only for cental LOS
    Brzt = eqm.rz2brzt(LOS_R[0].flatten()[:,None], LOS_Z[0].flatten()[:,None], eq_time*np.ones(LOS_Z[0].size))
    Br,Bz,Bt = np.reshape(Brzt, (3,)+LOS_Z[0].shape)

    #orthogonal projections of Br and Bt to x,y directions
    #DIII-D standart Bt direction is in CW direction
    By =  Br*np.sin(Phi) - Bt*np.cos(Phi)
    Bx =  Br*np.cos(Phi) +Bt*np.sin(Phi)
 
    #unit magnetic vector along LOS
    unitBlos  = np.array((Bx,By,Bz))
    unitBlos /= np.linalg.norm(unitBlos,axis=0)
    parallelity = np.abs(np.sum(unitBlos*losDir[0,:,:,None],0)) #cos of angle
 
    maxParallelity = parallelity.max(axis=1)

    par_tg = parallelity[:,los_len//2]
 

    #weighted sum all dL values which fits in between bin edges
    #dLmat2 is calculate by summing all dL contributiions to each grid bin, used for benchmarking 
    dLmat2 = np.zeros((nv,nr))
    #calculate exactly length of chord in each grid bin - it is accurate but also more difficult 
    dLmat  = np.zeros((nv,nr))

    #iterate over chords and calculate contribution to each grid bin
    for ilos in range(nv):
        #weight is given by dL value and is is splitted equally between all nvirt virtual LOSs
        weights = np.tile(dL[:,ilos], (los_len,1))/nvirt
        weights *= parallelity[ilos,:,None]**n#multiply by cos of angle between LOS and mag. field line
        dLmat2[ilos],edges = np.histogram(LOS_rho[:,ilos].T.flatten(), rho_grid_edges,
                                            weights=weights.flatten(), density=False)


        for iv in range(nvirt):
            #get first and last point just outside of lcfs
            ind = np.where(np.abs(LOS_rho[iv,ilos]) < 1)[0]
            if len(ind) == 0: #no crossection with grid
                continue
            ilcfs_in  = ind[0]-1
            ilcfs_out = ind[-1]+2
            rho_cut = LOS_rho[iv,ilos,ilcfs_in:ilcfs_out]

            #find local minima and maxima, rho is monotonous in between        
            imin = argrelextrema(rho_cut,np.less_equal)[0]
            imax = argrelextrema(rho_cut,np.greater_equal)[0]
            i_extrema = np.unique(np.r_[imin,imax])

            L = np.arange(ilcfs_out-ilcfs_in)*dL[iv,ilos]
            L_turn = 0

            for i in range(len(i_extrema)-1):
                #split LOS in regions with monotonously changing rho to make the inversion
                monotone_ind = slice(i_extrema[i],i_extrema[i+1]+1)

                if rho_cut[monotone_ind][0]== rho_cut[monotone_ind][-1]:
                    #special case, whole monotone_ind is constant
                    continue

                #input for interpolation needs to be monotonous
                Ledge = interp1d(rho_cut[monotone_ind], L[monotone_ind],
                                 bounds_error=False, fill_value = 0)(rho_grid_edges)
                par = interp1d(rho_cut[monotone_ind], 
                               parallelity[ilos,ilcfs_in:ilcfs_out][monotone_ind],
                                fill_value = 'extrapolate')(rho_grid_centers)            


                if not np.any(Ledge>0):
                    #special case, whole monotone_ind is within single grid cell
                    continue

                weight = np.abs(par)**n/nvirt

                #index of turning points
                i_tg_in, i_tg_out = np.where(Ledge>0)[0][[0,-1]]

                dLmat[ilos, i_tg_in:i_tg_out] += np.abs(Ledge[i_tg_in:i_tg_out]-Ledge[i_tg_in+1:i_tg_out+1])*weight[i_tg_in:i_tg_out]

                #add the turning point singularity
                if  rho_cut[monotone_ind][-1]-rho_cut[monotone_ind][0] > 0: # if L increases with rho
                    if i_tg_in > 0:
                        dLmat[ilos,i_tg_in-1] +=  (Ledge[i_tg_in]-L_turn)*weight[i_tg_in-1]
                    if i_tg_out < nr: #this should happen just in the last step
                        L_turn = Ledge[i_tg_out]
                        assert i < len(i_extrema)-1
                else:   # if L decreases with rho
                    if i_tg_out < nr: #this should happen just in the first step
                        dLmat[ilos,i_tg_out] +=  (Ledge[i_tg_out]-L_turn)*weight[i_tg_out-1]
                    L_turn = Ledge[i_tg_in]


    #assume that LOS which has not crossed the plasma contributes only to the outermost point
    los_inside_plasma = np.any(dLmat!=0,axis=1)
    #select LOS where the LOS is inside of the plasma
    los_inside_plasma = los_inside_plasma[chords]
    chords = chords[los_inside_plasma]

    np.savez(saveFilename,L = dLmat[chords, :], L2 = dLmat2[chords, :],
             rho_center=rho_grid_centers,
             rho_edge=rho_grid_edges,chords=chords+1,los_inside_plasma=los_inside_plasma,
             n=n,R_tg=R_tg[:,chords].mean(0), Z_tg = Z_tg[:,chords][0], L_tg = L_tg[:,chords][0],rho_tg=rho_tg[chords], maxParallelity = maxParallelity[chords])
