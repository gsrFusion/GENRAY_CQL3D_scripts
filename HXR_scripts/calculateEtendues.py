# make sure to have an interactive backend
import sys
import os
import numpy as np
from numpy.linalg import norm
from numpy import sqrt
import matplotlib.pyplot as plt
import constants
import tofu as tf
import DetectorInformation

print("past imports")

# simpler config => faster computation and plotting
conf = tf.load_config('AUG-V1')
# instanciate a Diagnostic
coll = tf.data.Collection()

############################
#Set up pinhole
############################
# Names aperature object
key_ap = 'ap'
# Defines aperature object
# Here we do a square slit
dap = {
    # Defines aperature geometry
    'dgeom': {
        'cent': np.array([constants.pinholeX_cql_m, 0, 0]), # [m], center in global axis
        'nin': np.array([-1,0,0]), # inward normal vec
        'e0': np.array([0,-1,0]), # horizontal vec
        'e1': np.array([0, 0, 1]), # vertical vec
        'outline_x0': np.array([-constants.pinholeRadius_m,-constants.pinholeRadius_m,constants.pinholeRadius_m,constants.pinholeRadius_m]), # [m], hor. disp. of corners wrt cent
        'outline_x1': np.array([-constants.pinholeRadius_m,constants.pinholeRadius_m,constants.pinholeRadius_m,-constants.pinholeRadius_m]), # [m], ver. disp. of corners wrt cent
    }
}
coll.add_aperture(
    key=key_ap,
    **dap['dgeom'],
    )

etendues = np.zeros(66)    


############################
#Set up pixels
############################

#we get the lines of sight by taking the vector between the collimator block opening and the pinhole center
#this does not exactly match how the Carlos' GRI code gets the LOSs
#the differences should be relatively small
for i in np.arange(1,66+1,1):
    theta, phi = DetectorInformation.getChordAngles_CQL(i)

    holeLoc = DetectorInformation.getHoleLocation_cm(i)/1e2
    #chord origin where the location is taken from the collimator block
    chordOrigin = np.array([constants.pinholeX_cql_m + constants.pinToCollBackplate_m, holeLoc[0],holeLoc[1]])
    los = np.array([constants.pinholeX_cql_m,0,0]) - chordOrigin
    unitLos = los/norm(los)

    insertedChordOrigin = constants.distanceInserted*unitLos + chordOrigin


    theta = np.arctan2(np.sqrt(los[0]**2 + los[1]**2), los[2]) % (2*np.pi)
    phi = np.arctan2(los[1], los[0]) % (2*np.pi)
    
    R_y = np.array([[np.cos(theta), 0, np.sin(theta)],[0 ,1, 0],[-np.sin(theta), 0 , np.cos(theta)]])
    R_z = np.array([[np.cos(phi), -np.sin(phi),0], [np.sin(phi), np.cos(phi), 0],[0,0,1]])
    e_0 = np.matmul(R_z, np.matmul(R_y,np.array([1,0,0])))
    

    e_1 = np.cross(los, e_0)

    """
    cqlChordOrigin = DetectorInformation.getChordOrigin_CQL_cm(i)/1e2
    los = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])

    R_y = np.array([[np.cos(theta), 0, np.sin(theta)],[0 ,1, 0],[-np.sin(theta), 0 , np.cos(theta)]])
    R_z = np.array([[np.cos(phi), -np.sin(phi),0], [np.sin(phi), np.cos(phi), 0],[0,0,1]])
    e_0 = np.matmul(R_z, np.matmul(R_y,np.array([1/sqrt(2),1/sqrt(2),0])))

    e_1 = np.cross(los, e_0)
    """

    key_pixel = str(i)
    dcam1d = {
        # Defines misc
        'dmics': {
            'name': 'SPEAR',
            'manufacturer': 'Kromek',
        },

        # Defines camera geometry
        'dgeom': {
            'cents_x': np.array(insertedChordOrigin[0]),
            'cents_y': np.array(insertedChordOrigin[1]),
            'cents_z': np.array(insertedChordOrigin[2]),
            'nin_x': np.array(los[0]), 'nin_y': np.array(los[1]), 'nin_z': np.array(los[2]),
            'e0_x': np.array(e_0[0]), 'e0_y': np.array(e_0[1]), 'e0_z': np.array(e_0[2]),
            'e1_x': np.array(e_1[0]), 'e1_y': np.array(e_1[1]), 'e1_z': np.array(e_1[2]),
            'outline_x0': np.array([-constants.detectorRadius_m,-constants.detectorRadius_m,constants.detectorRadius_m,constants.detectorRadius_m]), # [m], hor. disp. of corners wrt cent
            'outline_x1': np.array([-constants.detectorRadius_m,constants.detectorRadius_m,constants.detectorRadius_m,-constants.detectorRadius_m]), # [m], ver. disp. of corners wrt cent, 
        },

        # Defines material
        'dmat': {
            'name': 'CZT', # detector material name
            'symbol': 'CZT', # detector material symbol
            'thickness': 5e-3, # [m], material thickness
            #'energy': None,
            'qeff': None, # quantum efficiency
        },
    }


    coll.add_camera_1d(
        key=key_pixel,
        dgeom=dcam1d['dgeom'],
        dmat=dcam1d['dmat'],
        )

    #Add pixel to diagnostic collection
    optics_i = [key_pixel,key_ap]
    coll.add_diagnostic(
    key='detector'+str(i),
    doptics=optics_i, # list of optics configuration
    compute=True, # compute etendues
    config=conf, # tokamak
    )
    etendues[i-1] = coll.ddata[f'detector{i}_{key_pixel}_etend']['data']*1e4

print(f"etendues: {np.round(etendues,10).tolist()}")
"""
from tofu.data import _class8_etendue_los as _etendue_los
zeroth = np.zeros(len(etendues))
first = np.zeros(len(etendues))
second = np.zeros(len(etendues))

for i in range(len(etendues)):
    dcompute,store = _etendue_los.compute_etendue_los(
        coll = coll,
        key = f'detector{i+1}',
        analytical = False,
        numerical = True,
        res = None,
        check = False,
        margin_par = None,
        margin_perp = None,
        add_points = None,
        convex=None,
        verb = None,
        plot = False,
        store = False,)

    zeroth[i] = dcompute[str(i+1)]['analytical'][0]
    first[i] = dcompute[str(i+1)]['analytical'][1]
    second[i] = dcompute[str(i+1)]['analytical'][2]

fig,ax = plt.subplots()
ax.plot(zeroth, label = 'zeroth')
ax.plot(first, label = 'first')
ax.plot(second, label = 'second')
ax.legend()
plt.show()
"""
#dax = coll.plot_diagnostic("diag0", plot_config=conf, data="etendue")
#plt.show()
