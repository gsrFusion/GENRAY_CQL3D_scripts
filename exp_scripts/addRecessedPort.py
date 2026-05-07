###
# This script is what I used to add the recessed port for ECE
# it originally lived on the MFE workstation, but should work here too
# It's not overly elegant. If the DIII-D wall changes at all, this will probably not work since it relies on knowing where in the wall array the new points need to be inserted
###

import numpy as np
import matplotlib.pyplot as plt
import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from omfit_classes import omfit_eqdsk
import getTargetInfo
import getGfileDict

np.set_printoptions(linewidth=np.inf,suppress=True)

targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

topmostShotDir = getTargetInfo.getTopmostShotDir()

unmodifiedEqdskName = f'g{shotNum}'
modifiedEqdskName = f'g{shotNum}_port'

os.system(f'cp {topmostShotDir}/{unmodifiedEqdskName} {topmostShotDir}/{modifiedEqdskName}')

geqdsk = omfit_eqdsk.OMFITgeqdsk(f'{topmostShotDir}/{modifiedEqdskName}')
#gfileDict = getGfileDict.getGfileDict(f'{topmostShotDir}/{modifiedEqdskName}')

xlim = geqdsk['RLIM']
ylim = geqdsk['ZLIM']

print(np.round(xlim,5))
print(np.round(ylim,5))
print(len(xlim))

#this will need to be modified once DIII-D changes its wall again
newxlim = np.concatenate([xlim[:42],[xlim[43],2.525,2.525,xlim[44]], xlim[45:]])
newylim = np.concatenate([ylim[:42],[.2,.2,-.2,-.2], ylim[45:]])

fig,ax = plt.subplots()
ax.plot(xlim,ylim)
ax.scatter(xlim,ylim)
ax.plot(newxlim,newylim)
ax.scatter(newxlim,newylim)
ax.set_aspect('equal')
plt.show()

geqdsk['RLIM'] = newxlim
geqdsk['ZLIM'] = newylim
geqdsk['LIMITR']=len(geqdsk['RLIM'])

geqdsk.save()
