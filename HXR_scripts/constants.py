import numpy as np

pinholeRadius_m = 5e-3#m Half sidelength of the pinhole at its smallest
detectorRadius_m = 2.5e-3#m Half sidelength of the detector pixel

pinToCollBackplate_m = (220.345)*1e-3 #m  Distance from the pinhole at its smallest to the back of the collimator plate # from GRI CAD
distanceInserted = -2.6e-3 #m distance the detector is inserted past the collimator backplate

pinholeX_cql_m = 2.93 #pinhole location in the cql coordinate system
pinholeX_DIIID_m = 143.3/1e2 #pinhole location in the DIII-D coordinate system
pinholeY_DIIID_m = 255.8/1e2 #pinhole location in the DIII-D coordinate system

#clockwise toroidal angle between the pinhole location in DIII-D coordinates and CQL coordinates in radians
clockwiseAngleShift_DIIIDtoCQL_rad = np.arctan(pinholeY_DIIID_m/pinholeX_DIIID_m)
#rotation matrix to convert a DIII-D coordinate to a CQL coordinate
DIIIDtoCQLRotMat = np.array([[np.cos(-clockwiseAngleShift_DIIIDtoCQL_rad),-np.sin(-clockwiseAngleShift_DIIIDtoCQL_rad), 0],
                                [np.sin(-clockwiseAngleShift_DIIIDtoCQL_rad), np.cos(-clockwiseAngleShift_DIIIDtoCQL_rad), 0],
                                [0,0,1]])

LOS_62_DIIID=np.array([-0.802, -0.597, 0.009]) #line of sight vector of chord 62, the central chord, in DIIID coordinates
