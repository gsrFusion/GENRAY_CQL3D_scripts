import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

import DetectorInformation
import netCDF4

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import getTargetInfo
targetDir = getTargetInfo.getTargetDir()
print(targetDir)
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

import getInputFileDictionary
cqlinput = getInputFileDictionary.getInputFileDictionary('cql3d')
import getGfileDict
gfileDict = getGfileDict.getGfileDict()


#polar thetas as measured from the z axis
thet1 = cqlinput['setup']['thet1']*np.pi/180.
#toroidal thetas as measured from the x axis
thet2 = cqlinput['setup']['thet2']*np.pi/180.
#location of XR detector
x_sxr = cqlinput['setup']['x_sxr']/100.  # [m]
y_sxr = 0#by convetion of CQL3D
z_sxr = cqlinput['setup']['z_sxr']/100.  # [m]
#number of sight lines
nv = cqlinput['setup']['nv']

_all = np.arange(1,nv+1,1)
_imageBand = np.array([2,12,22,32,42,52,62, 3,13,23,33,43,53,63,  4,14,24,34,44,54,64,  5,15,25,35,45,55,65])
chordsToPlot = _imageBand

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 17)
plt.rc('figure', titlesize = 16)
plt.rc('legend',fontsize=20)

#plots the toroidal and poloidal lines of sight of the detectors
def plotLos():
    fig, ax = plt.subplots(dpi = 100)
    #fig.set_size_inches(6, 12)
    
    for chordNum in chordsToPlot:
        leftAngle, rightAngle, _, __ = DetectorInformation.getDivergenceAngles(chordNum)
        addLos(ax, x_sxr,y_sxr,thet2[chordNum-1], thet1[chordNum-1], leftAngle, rightAngle)

    setupToroidalLos(fig, ax)

    fig.tight_layout(rect=[0, 0.0, 1, 0.97])
    plt.show()

def setupToroidalLos(fig, ax):
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")

    ax.set_ylim([-.5, 2.6])
    ax.set_xlim([-2.45,1.5])
    ax.set_xticks([-2,-1,0,1])
    ax.set_yticks([0,1,2])
    ax.set_aspect('equal', adjustable='box')

    ax.yaxis.set_tick_params(which='major', width=2, length=5)
    ax.xaxis.set_tick_params(which='major', width=2, length=5)

    ax.add_patch(plt.Circle((0,0), 1, color = 'k', fill = False, lw = 2))
    ax.add_patch(plt.Circle((0,0), 1.67+.67, color = 'k', fill = False, lw = 2))

#adds the toroidal and poloidal lines of sight of a particular detector
#the toroidal line of sight is simply a line, but projecting the line of sight onto the poloidal plane is more complicated
#It is done by starting at the detector and following the derivative of the line while recording the associated (r,z) point
def addLos(ax, xDec, yDec, torTheta, polarTheta, leftAngle, rightAngle):
    #direction of the line of sight
    losDirCenter = np.array([np.cos(torTheta)*np.sin(polarTheta), np.sin(torTheta)*np.sin(polarTheta), np.cos(polarTheta)])
    losDirLeft = np.array([np.cos(torTheta + leftAngle)*np.sin(polarTheta), np.sin(torTheta+ leftAngle)*np.sin(polarTheta), np.cos(polarTheta)])
    losDirRight = np.array([np.cos(torTheta - rightAngle)*np.sin(polarTheta), np.sin(torTheta - rightAngle)*np.sin(polarTheta), np.cos(polarTheta)])

    currentX = np.array([x_sxr, x_sxr, x_sxr])
    currentY = np.array([y_sxr, y_sxr, y_sxr])
    currentZ = z_sxr
    
    centerXs = [currentX[0]]; leftXs = [currentX[1]]; rightXs = [currentX[2]];
    centerYs = [currentY[0]]; leftYs = [currentY[1]]; rightYs = [currentY[2]];
    centerZs = [currentZ]

    minR = 1
    maxR = np.sqrt(x_sxr**2 + y_sxr**2)
    boundingZ = 1.12
    majorRad = 1.67+.67

    dx = (maxR-minR)/100

    currentR = np.sqrt(currentX**2 + currentY**2)
    #while the current point is still inside a reasonable volume
    while(np.abs(currentZ) < boundingZ and minR < currentR[0] <= maxR):
        newX = -dx+currentX
        newY = -dx*np.array([losDirCenter[1]/losDirCenter[0], losDirLeft[1]/losDirLeft[0], losDirRight[1]/losDirRight[0]]) + currentY
        newZ = losDirCenter[2]*(-dx/losDirCenter[0]) + currentZ

        currentX = newX
        currentY = newY
        currentZ = newZ

        centerXs.append(newX[0]); leftXs.append(newX[1]); rightXs.append(newX[2])
        centerYs.append(newY[0]); leftYs.append(newY[1]); rightYs.append(newY[2])
        centerZs.append(newZ)

        currentR = np.sqrt(currentX**2 + currentY**2)
    centerXs = np.array(centerXs);centerYs = np.array(centerYs)
    leftXs = np.array(leftXs); leftYs = np.array(leftYs)
    rightXs = np.array(rightXs); rightYs = np.array(rightYs)
    #listR = np.sqrt(np.array(listx)**2 + np.array(listy)**2)  

    leftXs_rot = leftXs*np.cos(1.061) - leftYs*np.sin(1.061)
    leftYs_rot = leftXs*np.sin(1.061) + leftYs*np.cos(1.061)
    rightXs_rot = rightXs*np.cos(1.061) - rightYs*np.sin(1.061)
    rightYs_rot = rightXs*np.sin(1.061) + rightYs*np.cos(1.061)

    ax.fill_between(leftXs_rot, leftYs_rot, rightYs_rot, alpha=.5, facecolor='b', edgecolor='None')


    #axs[1].plot(leftXs, leftYs, color = 'g'); axs[1].plot(rightXs, rightYs, color = 'b')
    
 
plotLos()
