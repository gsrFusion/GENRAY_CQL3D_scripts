import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.ticker as tkr
from matplotlib.patches import Polygon
from scipy.interpolate import interp2d
import matplotlib.cm as cm
import matplotlib
import netCDF4
import os, sys

#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

import DetectorInformation
import CountMatrix
import constants
import BuildLmatrix

currentDir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentDir)

cql_nc = netCDF4.Dataset(f'{parentdir}/cql3d.nc','r')
cqlrf_nc = netCDF4.Dataset(f'{parentdir}/cql3d_krf001.nc','r')
sys.path.append(parentdir)
import getInputFileDictionary
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d',pathprefix=f'{parentdir}/')
import getGfileDict
gfileDict = getGfileDict.getGfileDict(pathprefix=f'{parentdir}/')

nv = cqlinput['setup']['nv']#number of detectors
emin = cqlinput['setup']['enmin']#min photon energy simulated
emax = cqlinput['setup']['enmax']#max photon energy simulated

en_ = np.ma.getdata(cql_nc.variables["en_"][:]) #energy bins used by CQL3D for the XR detector
eflux = cql_nc.variables["eflux"][:] #count rates of the chords. eflux[0] is thermal bremsstrahlung and eflux[1] is the nonthermal 
thet1 = np.array(cqlinput['setup']['thet1'])*np.pi/180 #polar thetas as measured from the vertical axis
thet2 = np.array(cqlinput['setup']['thet2'])*np.pi/180 #toroidal thetas as measured from the x axis

#location of XR detector
x_sxr = cqlinput['setup']['x_sxr']/100.  # [m]
y_sxr = 0#by convetion of CQL3D
z_sxr = cqlinput['setup']['z_sxr']/100.  # [m]

#Plots the ray traces and spot sizes of the detectors
#Color of the spot size corresponds to the number of counts it measures
def plotRaysAndCounts():
    fig,ax = plt.subplots(figsize = (6,9), dpi = 100)
    addRays(ax)
    drawFluxSurfaces(ax)
    addAllSightlines(ax,fig)

    #####Setup of plots#####
    ax.set_ylim([-1.4,1.4])
    ax.set_xlim([.95,2.4])
    ax.set_ylabel("z (m)")
    ax.set_xlabel("R (m)")
    ax.set_aspect(1)
    plt.rc('xtick', labelsize = 18)
    plt.rc('ytick', labelsize = 18)
    plt.rc('axes', labelsize = 20)
    plt.rc('figure', titlesize = 18)
    fig.tight_layout()
    plt.show()
    ########################

#adds the ray traces to ax
def addRays(ax):
    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    rbbbs = gfileDict["rbbbs"] #R points of the LCFS
    zbbbs = gfileDict["zbbbs"] # Z points of the LCFS
    
    wr  = cqlrf_nc.variables["wr"][:] #major radius of the ray at each point along the trace
    wz  = cqlrf_nc.variables["wz"][:] #height of the ray at each point along the trace
    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    wr *= .01; wz*=.01 #convert to m from cm
    
    maxDelPwrPlot = .9 #what portion of ray power must have been damped before we stop plotting that ray

    norm = plt.Normalize(0, 1)

    #plot the ray using a LineCollection which allows the colormap to be applied to each ray
    for ray in range(len(wr)):
        
        mostPowerDep = findNearestIndex(1 - maxDelPwrPlot, delpwr[ray]) #find the index of the last ray point we want to plot
        delpwr[ray,:] = delpwr[ray,:]/delpwr[ray,0] #normalize the ray power to that ray's starting power
        points = np.array([wr[ray][:mostPowerDep], wz[ray][:mostPowerDep]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Create a continuous norm to map from data points to colors
        lc = LineCollection(segments, norm = norm,cmap=plt.cm.jet)
        # Set the values used for colormapping
        lc.set_array(delpwr[ray][:mostPowerDep])
        lc.set_linewidth(1)
        ax.add_collection(lc)

    ax.plot(xlim, ylim, 'r', lw = 2)#plot wall
    ax.plot(rbbbs, zbbbs, 'k', lw = 1.5)#plot LCFS

#this adds all of the sightlines, specifically their spot sizes at their position of greatest tangency to the magnetic field
def addAllSightlines(ax, fig):
    #Choose which set of chords to plot
    #all of the chosen chords must interset the plasma in order for the L_matrix to contain them
    _all = np.arange(1,nv+1,1).astype(int)
    _imageBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])
    _HFS = np.array([42,52,62, 43,53,63,  44,54,64, 45,55,65, 41])
    _LFS = np.array([2,12,22,32, 3,13,23,33,  4,14,24,34, 5,15,25,35, 11,21,31,])
    _upperBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64, 5,15,25,35,45,55,65, 11,21,31,41])+1
    chordsToPlot = _imageBand
    
    #build a new length matrix to get all the chord geometry we need
    BuildLmatrix.main(chordsToPlot,.8, saveFilename = 'lengthMatrices/Lmat147634_plotRaysAndCounts')
    L_matrix = np.load(f'lengthMatrices/Lmat147634_plotRaysAndCounts.npz')

    chords = L_matrix['chords']
    R_tgs = L_matrix['R_tg'] #Rs at which the sightlines are most tangent to B
    Z_tgs = L_matrix['Z_tg'] #Zs at which the sightlines are most tangent to B
    L_tgs = L_matrix['L_tg'] #distances along the sightlines at which the sightlines are most tangent to B
    maxParallelities = L_matrix['maxParallelity'] #maximum value of the dot product between detector unit vector and magnetic field unit vector
    tangencyAngles = np.arccos(maxParallelities)*180/np.pi
    print(chords[np.argmax(maxParallelities)])
    leastTangent = np.max(tangencyAngles)#angle of the sightline of greatest tangency to B

    #add the sightlines one by one
    for i in range(len(chords)):
        addParticularSightline(ax, chords[i], tangencyAngles[i], L_tgs[i], R_tgs[i], Z_tgs[i], leastTangent)


    ###### Set up colorbar ######
    cmap = matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(0,leastTangent),
         cmap = plt.get_cmap('viridis'))
    cmap.set_array([])
    ticks = np.linspace(0,leastTangent,5)

    formatter = tkr.ScalarFormatter(useMathText=True)
    #formatter.set_powerlimits((0,0))

    cbar = fig.colorbar(cmap, ax = ax, shrink = .7, ticks = ticks, format = formatter)
    cbar.set_label(r"$\theta = $ acos($\hat{k} \cdot \hat{B}$) (degrees)")
    plt.rc('font', **{'size':'10'})
    #############################

#adds the toroidal and poloidal lines of sight of a particular detector
#the toroidal line of sight is simply a line, but projecting the line of sight onto the poloidal plane is more complicated
#It is done by starting at the detector and following the derivative of the line while recording the associated (r,z) point
"""Needs to be review - etendues and hence angles are wrong"""
def addParticularSightline(ax, chordNum, tangencyAngle, L_tg, R_tg, Z_tg, leastTangent):
    cornerAnglePairs = DetectorInformation.getCornerAngles(chordNum, giveTotalAngle = False) 
    xMults = [1,-1,-1,1]; yMults = [-1,-1,1,1]
    
    pinholeCorners = getPinholeCorners_cql3d()#pinhole corners in cql coordinates
    crossSecCorners = np.zeros((4,2))

    torTheta = thet2[chordNum-1]; polarTheta = thet1[chordNum-1]
    losCenter = np.array(np.array([np.cos(torTheta)*np.sin(polarTheta),
                        np.sin(torTheta)*np.sin(polarTheta), 
                        np.cos(polarTheta)]))
    for j in range(len(cornerAnglePairs)):
        anglePair = cornerAnglePairs[j]
        xMult = xMults[j]; yMult = yMults[j]
        
        losCorner = np.array(np.array([np.cos(torTheta + xMult*anglePair[0])*np.sin(polarTheta+yMult*anglePair[1]),
                        np.sin(torTheta+xMult*anglePair[0])*np.sin(polarTheta + yMult*anglePair[1]), 
                        np.cos(polarTheta + yMult*anglePair[1])]))
        projectionLength = np.dot(losCorner, losCenter)
        scalingFactor = L_tg/projectionLength
        scaledVector = losCorner * scalingFactor

        cornerVec = scaledVector + pinholeCorners[j]
        cornerR = np.sqrt(cornerVec[0]**2 + cornerVec[1]**2)
        cornerZ = cornerVec[2]

        disps = L_tg*np.tan(anglePair)
        crossSecCorners[j][0] = cornerR
        crossSecCorners[j][1] = cornerZ

    color = cm.viridis(tangencyAngle/leastTangent)
    ax.add_patch(Polygon(crossSecCorners, facecolor = color, edgecolor = color, zorder = 3))
    textColor = 'w'
    if tangencyAngle > .85*leastTangent:
        #for when the colormap to show degree of tangency makes it hard to read the sightline number
        textColor = 'k'
    ax.text(R_tg, Z_tg, str(chordNum),
        color = textColor, horizontalalignment='center', verticalalignment='center', fontsize = 14)


#returns the corners of the pinhole in cql3d coordinates
"""Needs to be review - etendues and hence angles are wrong"""
def getPinholeCorners_cql3d():
    polarTheta = thet1[62-1]; torTheta = thet2[62-1]
    pinholeNormal = losDirCenter = np.array([np.cos(torTheta)*np.sin(polarTheta), np.sin(torTheta)*np.sin(polarTheta), np.cos(polarTheta)])
    pinholeRad = constants.pinholeRadius_m
    pinholeCorners = np.array([(pinholeRad,pinholeRad,0),
        (pinholeRad,-pinholeRad,0),
        (-pinholeRad,-pinholeRad,0),
        (-pinholeRad,pinholeRad,0)])
    
    rotMat = getRotationMatrix(-torTheta,polarTheta,0)

    pinholeCornersCQLCoordinates = np.array([np.matmul(rotMat, pinholeCorner) for pinholeCorner in pinholeCorners]) + np.array([x_sxr, y_sxr, z_sxr])
    return pinholeCornersCQLCoordinates

#Rotation matrix for the improper Euler angles
#See https://en.wikipedia.org/wiki/Rotation_matrix
def getRotationMatrix(alpha, beta, gamma):
    from numpy import sin, cos
    mat = np.zeros((3,3))
    mat[0,:] = [cos(beta)*cos(gamma), sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma), cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)]
    mat[1,:] = [cos(beta)*sin(gamma), sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)]
    mat[2,:] = [-sin(beta), sin(alpha)*cos(beta), cos(alpha)*cos(beta)]

    return mat


#draws poloidal flux surfaces according to the levels passed to ax.contour
def drawFluxSurfaces(ax):
    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 200)
    zInterp = np.linspace(np.min(z), np.max(z), 200)
    psirzNormInterp = interp2d(r,z, psirzNorm, kind = 'cubic')(rInterp, zInterp)
    
    ax.contour(rInterp, zInterp, psirzNormInterp, np.square(np.arange(0,1,.1)), colors= 'k', linewidths= 1.5)

#returns the index of the array whose element is closest to value
def findNearestIndex(value, array):
    idx = (np.abs(array - value)).argmin()

    return idx


plotRaysAndCounts()
