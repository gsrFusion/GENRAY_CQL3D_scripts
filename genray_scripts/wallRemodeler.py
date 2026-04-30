import numpy as np
import matplotlib.pyplot as plt
import getGfileDict
from omfit_classes import omfit_eqdsk
np.set_printoptions(linewidth=np.inf,suppress=True)

shotNum = '206629.01980'
if shotNum == '203619.04130':
    gfileDict = getGfileDict.getGfileDict('/home/grantr/Desktop/ECEplots/g203619.04129_deepPort')

if shotNum == '203912.02700':
    gfileDict = getGfileDict.getGfileDict('/home/grantr/Desktop/ECEplots/g203912.02700_deepPort')

if shotNum == '203912.02780':
    gfileDict = getGfileDict.getGfileDict('/home/grantr/Desktop/ECEplots/g203912.02780')

if shotNum == '206636.01960':
    gfileDict = getGfileDict.getGfileDict('/home/grantr/Desktop/ECEplots/g206636.01960')

if shotNum == '206629.01980':
    gfileDict = getGfileDict.getGfileDict('/home/grantr/Desktop/ECEplots/g206629.01980')

xlim = gfileDict['xlim']
ylim = gfileDict['ylim']

print(np.round(xlim,5))
print(np.round(ylim,5))
print(len(xlim))

newxlim = np.concatenate([xlim[:42],[xlim[43],2.525,2.525,xlim[44]], xlim[45:]])
newylim = np.concatenate([ylim[:42],[.2,.2,-.2,-.2], ylim[45:]])

fig,ax = plt.subplots()
ax.plot(xlim,ylim)
ax.scatter(xlim,ylim)
ax.plot(newxlim,newylim)
ax.scatter(newxlim,newylim)
ax.set_aspect('equal')
plt.show()

if shotNum == '203619.04130':
    geqdsk = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g203619.04129')
if shotNum == '203912.02700':
    geqdsk = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g203912.02700')
if shotNum == '203912.02780':
    geqdsk = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g203912.02780')
if shotNum == '206636.01960':
    geqdsk = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g206636.01960')
if shotNum == '206629.01980':
    geqdsk = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g206629.01980')

geqdsk['RLIM'] = newxlim
geqdsk['ZLIM'] = newylim
geqdsk['LIMITR']=len(geqdsk['RLIM'])

geqdsk.save()
if shotNum == '203619.04130':
    geqdsk_port = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g203619.04129')
if shotNum == '203912.02700':
    geqdsk_port = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g203912.02700')
if shotNum == '203912.02780':
    geqdsk_port = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g203912.02780')
if shotNum == '206636.01960':
    geqdsk_port = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g206636.01960')
if shotNum == '206629.01980':
    geqdsk_port = omfit_eqdsk.OMFITgeqdsk('/home/grantr/Desktop/ECEplots/g206629.01980')


fig,ax = plt.subplots()
geqdsk_port.plot()
#ax.scatter(geqdsk_port['RLIM'], geqdsk_port['ZLIM'])
#ax.set_aspect('equal')
plt.show()