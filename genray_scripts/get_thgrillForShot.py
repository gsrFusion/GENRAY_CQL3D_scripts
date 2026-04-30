"""
#This script determines what poloidal location to set the grill at
# thgrill in GENRAY is relative to the magnetic axis
# but obviously the grill is in the same physical position for each shot
# so this code basically converts between physical location to angle relative to the magnetic axis
"""

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import getGfileDict
gfileDict = getGfileDict.getGfileDict()
import numpy as np
import matplotlib.pyplot as plt

R_mag = gfileDict['rmaxis']
Z_mag = gfileDict['zmaxis']

grillR = 1.04241092
grillZ = -.11238484#-(1.67-grillR)*np.tan(np.radians(9))

dR = grillR - R_mag
dZ = grillZ - Z_mag

print(f'dR: {dR}, dZ: {dZ}')

thgrill = np.arctan2(dZ,dR)

if thgrill < 0:
    thgrill += 2*np.pi

#thgrill = np.pi/2
#if dZ > 0:
#    thgrill += np.pi/2 + np.tan(dZ/dR)
#else:
#    thgrill += np.tan(np.abs(dZ)/dR)

print(f'thgrill should be {np.degrees(thgrill)} deg')

fig,ax = plt.subplots()

ax.scatter([1.67, grillR, R_mag], [0, grillZ, Z_mag])
xlim = gfileDict["xlim"] #R points of the wall
ylim = gfileDict["ylim"] #Z points of the wall
ax.plot(xlim, ylim, 'r', lw = 2)#plot wall
ax.axhline(0,lw = 2, color ='k')
ax.set_aspect('equal')
plt.show()