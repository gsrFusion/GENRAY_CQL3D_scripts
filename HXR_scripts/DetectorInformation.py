import constants
import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt

import constants
pinholeRad = constants.pinholeRadius_m
detectRad = constants.detectorRadius_m
pinToCollBack = constants.pinToCollBackplate_m

#returns the etendue in units of str cm^2
#I've gone through a lot of estimations on etendue, as can be seen
#So to figure which one is correct, I passed my detector geometry to TOFU
#The two TOFU_etendues come from two different numbers for the distance between the back of the collimator block to the pinhole
#It appears TOFU_etendues_22 is the most correct
def getEtendues():
    #using distance between pinhole and collimator from the CAD
    #assumed a 2.5mm distance between detectors and collimator blackplate
    TOFU_etendues = np.array([0.0003201676, 0.0003416119, 0.0003584421, 0.0003692164, 0.0003729285, 0.0003692164, 0.0003584421, 0.0003416119, 0.0003201676, 0.0003321364, 0.0003584421, 0.0003805424, 0.0003965721, 0.0004050146, 0.0004050146, 0.0003965721, 0.0003805424, 0.0003584421, 0.0003321364, 0.0003692164, 0.0003965721, 0.0004182538, 0.0004322315, 0.0004370648, 0.0004322315, 0.0004182538, 0.0003965721, 0.0003692164, 0.0003729285, 0.0004050146, 0.0004322315, 0.00045212, 0.0004626446, 0.0004626446, 0.00045212, 0.0004322315, 0.0004050146, 0.0003729285, 0.0004050146, 0.0004370648, 0.0004626446, 0.0004792177, 0.0004849634, 0.0004792177, 0.0004626446, 0.0004370648, 0.0004050146, 0.0003965721, 0.0004322315, 0.0004626446, 0.0004849634, 0.000496806, 0.000496806, 0.0004849634, 0.0004626446, 0.0004322315, 0.0003965721, 0.0004182538, 0.00045212, 0.0004792177, 0.000496806, 0.0005029096, 0.000496806, 0.0004792177, 0.00045212, 0.0004182538])
    """
    sideToSideDecEtendues = np.zeros(66)
    sideToSideDecEtendues_alt = np.zeros(66)
    projectedPinholeAreaEtendues = np.zeros(66)
    BryanEtendues = np.zeros(66)

    for i in range(1,67):
        sideToSideDecEtendues[i-1] = np.pi*.5**2 * np.sin(np.sum(getLeftRightUpDownAngles(i))/4)**2
        sideToSideDecEtendues_alt[i-1] = np.pi*.5**2 * np.sin(np.sum(getLeftRightUpDownAngles_alt(i))/4)**2

    avgAvgAngle = 0  
    for i in range(1,67):
        phi = getAngleFromPinholeNormal(i)
        pinholeProjectedArea = (pinholeRad*100)**2*np.cos(phi) #in cm^2

        angles = getCornerAngles(i, giveTotalAngle = True)

        avgAngle = np.sum(angles)/len(angles)
        avgAvgAngle += avgAngle
        projectedPinholeAreaEtendues[i-1] = 2*np.pi*np.sin(avgAngle)**2*pinholeProjectedArea
        BryanEtendues[i-1] = 2*np.pi*np.sin(avgAngle)**2*(detectRad*100)**2
    projectedPinholeAreaEtendues *= 2/np.pi #approxmation to account for decrease in solid angle due to square cross section instead of circular
    """
    """
    fig,ax = plt.subplots()
    #ax.plot(projectedPinholeAreaEtendues, label = "corner to corner with projected pinhole")
    #plt.plot(sideToSideDecEtendues, label = "side to side angles, detector area", color = 'r',linestyle = 'dotted')#
    #plt.plot(sideToSideDecEtendues_alt, label = "side to side angles, detector area. Alt origin", color = 'g',linestyle = 'dotted')#
    #plt.plot(BryanEtendues, label = "Old approximation", color = 'b',linestyle = 'dotted')
    plt.plot(TOFU_etendues_22, label = "TOFU", lw = 2, color ='k')
    plt.plot(TOFU_etendues_22_20mminserted, label = "TOFU 20mm inserted", lw = 2, color ='r')
    plt.plot(TOFU_etendues_22_40mminserted, label = "TOFU 40mm inserted", lw = 2, color ='g')
    plt.plot(TOFU_etendues_22_60mminserted, label = "TOFU 60mm inserted", lw = 2, color ='b')
    ax.set_ylabel("Etendue (cm*str)");ax.set_xlabel("Chord number")
    fig.tight_layout()
    ax.legend()
    plt.show()

    """
    return TOFU_etendues
   
#returns the x,y position in cm of the sightline's hole in the plastic holder (ie back of comminator block)
#converted from Carlos' matlab script to python
def getHoleLocation_cm(chordNum):
 
    xspace = np.sqrt(3)/2*2   # in cm for a 2 cm hexagon
    yspace = 2              #3 in cm for a 2 cm hexagon

    nx=13;                                                 # number of columns
    xposarr=np.linspace(-6,6,nx)*xspace                        # position array from -6 to 6 x sqrt(3)
    xheaderarr=[1,10,20,29,39,48,58,67,77,86,96,105,115]   # corresponding headers

    xheaderpos=-1
    
    for k in range(0, len(xheaderarr)):
        if chordNum>=xheaderarr[k] and chordNum<xheaderarr[k+1]:
            xheaderpos = k
    if xheaderpos==-1:
        xheaderpos=len(xheaderarr)-1# Catch last column

    xheader=xheaderarr[xheaderpos]                        # corresponding header hole number 
    xpos=xposarr[xheaderpos]                              # OUTPUT: x position of requested hole, -6 to 6 holes
    
    yhalf=0;                                                # half added in below for even columns
    if xheaderpos % 2 == 0:
        yhalf=1

    ny=19;                                                  # number of rows, including half rows
    yposarr=np.linspace(4.5,-4.5,ny)*yspace           # position array from 4.5 to -4.5 x 2 cm
    ypos=yposarr[2*(chordNum-xheader)+yhalf]           # OUTPUT: y position of requested hole, -4.5 to 4.5 holes
    
    return np.array([xpos, ypos])

#returns the starting position of a chord in DIIID coordinates
#copied from Carlos' matlab script
def getChordOrigin_DIIID_cm(chordNum):
    So=np.array([constants.pinholeX_DIIID_m, constants.pinholeY_DIIID_m, 0])*1e2  # in cm, location of GRI pinhole in cartesian coordinates
    Lcent=pinToCollBack*1e2;                                               # in cm, length from back of block to pinhole on central chord

    xholepos,yholepos=getHoleLocation_cm(chordNum)   # in cm, retrieve the hole position of chord
    Scent=So - Lcent*constants.LOS_62_DIIID        # in cm, position of center of hole on back of block, along -pinholeX_DIIID_m
    Schord=Scent + xholepos*np.array([constants.LOS_62_DIIID[1], -constants.LOS_62_DIIID[0],0])  + np.array([0,0,yholepos]) # in cm, position of chord hole on back of block
    return Schord

#returns the starting position of a chord in CQL3D coordinates
def getChordOrigin_CQL_cm(chordNum):
    Schord=getChordOrigin_DIIID_cm(chordNum)
    return np.matmul(constants.DIIIDtoCQLRotMat, Schord)

#returns [theta, phi] of the line of sight in CQL3D coordinates
#theta and phi are the typical spherical
def getChordAngles_CQL(chordNum):
    So=np.array([constants.pinholeX_DIIID_m, constants.pinholeY_DIIID_m, 0])*1e2  # in cm, location of GRI pinhole in cartesian coordinates
    Schord=getChordOrigin_DIIID_cm(chordNum)

    Lhat_DIID=(So-Schord)/np.linalg.norm(So-Schord)           # unit vector of chord
    Lhat_CQL = np.matmul(constants.DIIIDtoCQLRotMat, Lhat_DIID)

    theta = np.arctan2(np.sqrt(Lhat_CQL[0]**2 + Lhat_CQL[1]**2), Lhat_CQL[2]) % (2*np.pi)
    phi = np.arctan2(Lhat_CQL[1], Lhat_CQL[0]) % (2*np.pi)
    return theta, phi

#plots the positions of the detectors in the collimator block as according to getChordOrigin_DIIID_cm
#plotting in the YZ plane
def plotDetectorPositionsInCollimator_CQL(chords):
    fig ,ax = plt.subplots()
    for chord in chords:
        chordLoc = getChordOrigin_CQL_cm(chord)
        ax.scatter([chordLoc[1]], [chordLoc[2]], s = 15, color = 'r')
        
        holeLoc = getHoleLocation_cm(chord)

        ax.scatter([holeLoc[0]], [holeLoc[1]], s = 15, color = 'g')

    plt.show()


#returns the angle between a chord's line of sight and the pinhole normal
def getAngleFromPinholeNormal(chordNum):
    holeLoc = getHoleLocation_cm(chordNum) *1e-2 
    detectorToPinhole = np.array([-holeLoc[0], -holeLoc[1], pinToCollBack])
    pinholeNormal = np.array([0,0,1])
    detectorToPinholeUnit = makeUnit(detectorToPinhole)
    return np.arccos(np.dot(pinholeNormal, detectorToPinholeUnit))
    
#normalizes the input vector
def makeUnit(vec):
    return vec/norm(vec)

#returns the pair of angles (i.e. vertical and horizontal) from each corner of the detector to the opposite corner of the pinhole
#considers finite detector
def getCornerAngles(chordNum, giveTotalAngle = True):
    #get xy location of detector on collimator backplate
    holeLoc = getHoleLocation_cm(chordNum) *1e-2
    #angles relative to the z direction (i.e. normal to collimator plate) so that detector is looking through pinhole center
    xAngle =np.arctan(holeLoc[0]/pinToCollBack)
    yAngle =np.arctan(holeLoc[1]/pinToCollBack)
    
    #define all the detector corners (currently not rotated so normal is pointed in +z)
    detectorCorners = np.array([(detectRad, detectRad,0),
        (detectRad, -detectRad,0),
        (-detectRad, -detectRad,0),
        (-detectRad, detectRad,0)])
        
    #rotate detector normal to face pinhole and give correct distance behind pinhole
    r = R.from_euler('zyx', [[0, -xAngle, yAngle]], degrees=False)
    for i in range(len(detectorCorners)):
        detectorCorners[i] = np.matmul(r.as_matrix(), detectorCorners[i])[0] + np.array([holeLoc[0],holeLoc[1],-pinToCollBack])
    
    #define pinhole corners   
    pinHoleCorners = np.array([(pinholeRad,pinholeRad,0),
        (pinholeRad,-pinholeRad,0),
        (-pinholeRad,-pinholeRad,0),
        (-pinholeRad,pinholeRad,0)])
    
    #distance from pinhole center to detector center
    detectCenterToPinCenter = np.array([-holeLoc[0], -holeLoc[1], pinToCollBack])
    
    #vectors pointing from each corner to the opposite diagonal corner
    #these are the vector that define the field of view
    cornerCornerVecs = np.copy(pinHoleCorners)
    for i in range(len(cornerCornerVecs)):
        cornerCornerVecs[i] = pinHoleCorners[i]-detectorCorners[(i+2)%4]

    #(xAngle, yAngle) between the detector normal and the edge of the pinhole
    angles = np.array([None]*4)
    
    for i in range(len(angles)):
        if giveTotalAngle:
            totalAngle = np.arccos(np.dot(detectCenterToPinCenter,cornerCornerVecs[i])/(norm(detectCenterToPinCenter)*norm(cornerCornerVecs[i])))
            angles[i] = totalAngle        
        else:
            xMask = [True, False, True]
            yMask = [False, True, True]
            xAngle = np.arccos(np.dot(detectCenterToPinCenter[xMask],cornerCornerVecs[i][xMask])/(norm(detectCenterToPinCenter[xMask])*norm(cornerCornerVecs[i][xMask])))
            yAngle = np.arccos(np.dot(detectCenterToPinCenter[yMask],cornerCornerVecs[i][yMask])/(norm(detectCenterToPinCenter[yMask])*norm(cornerCornerVecs[i][yMask])))
            angles[i] = (xAngle, yAngle)
            
    if giveTotalAngle:
        pass#print(angles)#print(f"chord: {chordNum}, {[np.degrees(angle) for angle in angles]}")
    else:
        pass
        #print(f"chord: {chordNum}, {[(np.degrees(angle[0]), np.degrees(angle[1])) for angle in angles]}")
            
    return angles 
    

#returns the angles between the detector and the left, right, top, and bottom of the pinhole
#assumes point-like detector
def getLeftRightUpDownAngles(chordNum):
    holeLoc = getHoleLocation_cm(chordNum) *1e-2
    holeLocAbs = np.abs(holeLoc)

    losCenter = np.array([-holeLoc[0], -holeLoc[1], pinToCollBack])
    losCenterX = np.copy(losCenter); losCenterX[1] = 0
    losCenterY = np.copy(losCenter); losCenterY[0] = 0
    losCenterXUnit = makeUnit(losCenterX); losCenterYUnit = makeUnit(losCenterY)

    detectorLoc = np.array([holeLoc[0], holeLoc[1], -pinToCollBack])

    leftTargetVec = np.array([-pinholeRad,0,0]) - detectorLoc
    leftTargetXZVec = np.copy(leftTargetVec); leftTargetXZVec[1] = 0;
    leftTargetXZUnitVec = makeUnit(leftTargetXZVec)
    
    rightTargetVec = np.array([pinholeRad,0,0]) - detectorLoc
    rightTargetXZVec = np.copy(rightTargetVec); rightTargetXZVec[1] = 0;
    rightTargetXZUnitVec = makeUnit(rightTargetXZVec)

    upTargetVec = np.array([0,pinholeRad,0]) - detectorLoc
    upTargetYZVec = np.copy(upTargetVec); upTargetYZVec[0] = 0;
    upTargetYZUnitVec = makeUnit(upTargetYZVec)

    downTargetVec = np.array([0,-pinholeRad,0]) - detectorLoc
    downTargetYZVec = np.copy(downTargetVec); downTargetYZVec[0] = 0;
    downTargetYZUnitVec = makeUnit(downTargetYZVec)

    leftAngle = np.arccos(np.dot(leftTargetXZUnitVec, losCenterXUnit))
    rightAngle = np.arccos(np.dot(rightTargetXZUnitVec, losCenterXUnit))
    upAngle = np.arccos(np.dot(upTargetYZUnitVec, losCenterYUnit))
    downAngle = np.arccos(np.dot(downTargetYZUnitVec, losCenterYUnit))

    return [leftAngle,rightAngle,upAngle,downAngle]


#returns the angles between the detector and the left, right, top, and bottom of the pinhole
#assumes point-like detector
#coordinates used here are relative to pinhole center point, not CQL3D or DIIID coordinates
def getLeftRightUpDownAngles_alt(chordNum):
    holeLoc = getHoleLocation_cm(chordNum)/1e2
    
    #z is in the direction of the pinhole normal
    detectorLoc = np.array([holeLoc[0], holeLoc[1], -pinToCollBack])
    #from detector center to pinhole center unit vector
    centToCentUnitVec = makeUnit(-detectorLoc)

    leftPinholeWall = np.array([-pinholeRad,0,0]); leftTargetVec = leftPinholeWall - detectorLoc
    rightPinholeWall = np.array([pinholeRad,0,0]); rightTargetVec = rightPinholeWall - detectorLoc
    topPinholeWall = np.array([0,pinholeRad,0]); topTargetVec = topPinholeWall - detectorLoc
    bottomPinholeWall = np.array([0,-pinholeRad,0]); bottomTargetVec = bottomPinholeWall - detectorLoc

    leftTargetUnit = makeUnit(leftTargetVec)
    rightTargetUnit = makeUnit(rightTargetVec)
    topTargetUnit = makeUnit(topTargetVec)
    bottomTargetUnit = makeUnit(bottomTargetVec)

    leftAngle = np.arccos(np.dot(centToCentUnitVec,leftTargetUnit))
    rightAngle = np.arccos(np.dot(centToCentUnitVec,rightTargetUnit))
    topAngle = np.arccos(np.dot(centToCentUnitVec,topTargetUnit))
    bottomAngle = np.arccos(np.dot(centToCentUnitVec,bottomTargetUnit))

    return [leftAngle, rightAngle, topAngle, bottomAngle]

#getEtendues()
#print(f"old: {getLeftRightUpDownAngles(32)}, new: {getLeftRightUpDownAngles_alt(32)}")
#plotDetectorPositionsInCollimator_CQL(np.arange(1,115,1))
