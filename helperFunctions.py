###
# Series of helper functions that are used is a variety of other scripts
###


import numpy as np
import getGfileDict
from scipy.interpolate import interp1d, interp2d
import matplotlib.pyplot as plt
import getTargetInfo
import getInputFileDictionary
import netCDF4
import shotToEqdsk

"""
#This script determines what poloidal location to set the grill at
# thgrill in GENRAY is relative to the magnetic axis
# but obviously the grill is in the same physical position for each shot
# so this code basically converts between physical location to angle relative to the magnetic axis
"""
def getThgrill(targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)

    R_mag = gfileDict['rmaxis']
    Z_mag = gfileDict['zmaxis']

    grillR = 1.04241092
    grillZ = -.11238484#-(1.67-grillR)*np.tan(np.radians(9))

    dR = grillR - R_mag
    dZ = grillZ - Z_mag

    thgrill = np.arctan2(dZ,dR)

    if thgrill < 0:
        thgrill += 2*np.pi

    print(f'thgrill should be {np.degrees(thgrill)} deg')

    fig,ax = plt.subplots()

    ax.scatter([1.67, grillR, R_mag], [0, grillZ, Z_mag])
    xlim = gfileDict["xlim"] #R points of the wall
    ylim = gfileDict["ylim"] #Z points of the wall
    ax.plot(xlim, ylim, 'r', lw = 2)#plot wall
    ax.axhline(0,lw = 2, color ='k')
    ax.set_aspect('equal')
    plt.show()

    return np.round(np.degrees(thgrill),4)

#returns the radial location of where the most power is deposited in the first pass
def getPeakFirstPassDepostion(cqlrf_nc, genray_in):
    radialBinEdges = np.linspace(0,1,51)
    radialBinCenters = (radialBinEdges[1:]+radialBinEdges[:-1])/2
    powerDep = np.zeros(len(radialBinCenters))
    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))
    nparas = cqlrf_nc.variables['wnpar'][:]
    for ray in range(len(nparas)):
        #if np.ma.getdata(nparas[ray][0]) == 0:
        #    return np.NAN
        if genray_in['grill']['anmax(1)'] >= nparas[ray][0] >= genray_in['grill']['anmin(1)']:
            mostPowerDep = findNearestIndex(1 - 0.99, delpwr[ray]/delpwr[ray][0]) #find the index of the last ray point we want to plot
            firstBounceIndex = findBounceIndex(radialVariable[ray][:mostPowerDep],bounceToFind = 1)
            radialVariableCenters = (radialVariable[ray][:firstBounceIndex][1:] + radialVariable[ray][:firstBounceIndex][:-1])/2
            indices = np.digitize(radialVariableCenters, radialBinEdges, right = False)
            indices[indices>49] = 49
            powerDep[indices-1] += np.abs(np.diff( delpwr[ray][:firstBounceIndex]))

    return radialBinCenters[np.argmax(powerDep)]

#returns the width in meters of the largest current peak
def getJcdWidth(targetDir = None, fracOfPeak = 0.5, mainPeak = True,
                R_lfs = None):

    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()

    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')

    curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2

    if R_lfs is None:
        rya = np.ma.getdata(cql_nc.variables["rya"][:])
        R_lfs = convertRhopolToRmidplane(rya, targetDir, side = 'LFS')

    R_lfs_interp = np.linspace(R_lfs[0], R_lfs[-1],500)

    curr_interp = interp1d(R_lfs, curr)(R_lfs_interp)
    normed = curr_interp / np.max(curr_interp)

    if mainPeak:
        maxPeakIndex = np.argmax(normed)
        """
        indicesBelowThresh = np.where(normed <= fracOfPeak)[0]

        leftIndex = indicesBelowThresh[indicesBelowThresh < maxPeakIndex][-1]
        rightIndex = indicesBelowThresh[indicesBelowThresh > maxPeakIndex][0]
        """
        maxPeakHeight = curr_interp[maxPeakIndex]
        leftIndex = maxPeakIndex
        rightIndex = maxPeakIndex

        while normed[leftIndex] > fracOfPeak or normed[rightIndex] > fracOfPeak:
            if normed[leftIndex] > fracOfPeak:
                leftIndex -= 1
            if normed[rightIndex] > fracOfPeak:
                rightIndex += 1
        #"""
        width = R_lfs_interp[rightIndex] - R_lfs_interp[leftIndex]

    else:
        width = R_lfs_interp[normed >= fracOfPeak][-1] - R_lfs_interp[normed >= fracOfPeak][0]
    return width

#returns int(rho *J * dA)/int(J * dA)
def getAvgCurrentLocAndTotal(targetDir):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    
    curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
    rya = np.ma.getdata(cql_nc.variables["rya"][:])
    darea = np.ma.getdata(cql_nc.variables["darea"][:])/1e4#convert to m^2

    totalCD = np.sum(curr*darea)

    print(f'loc of max J: {rya[np.argmax(curr)]}, weighted max: {np.sum(curr*rya*darea)/totalCD}')

    """
    fig,ax = plt.subplots()
    ax.plot(rya, curr)
    ax.axvline(rya[np.argmax(curr)])
    ax.axvline(np.sum(curr*rya*darea)/np.sum(curr*darea))
    plt.show()
    """

    return np.sum(curr*rya*darea)/totalCD, totalCD

#returns the SPA of the forward lobe
def getSPA(targetDir):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()


    genray_in = None
    genray_in = getInputFileDictionary.getInputFileDictionary('genray_LH', targetDir=targetDir)
    cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')

    ngrill = genray_in['grill']['ngrill']
    lobes = np.arange(1,ngrill+1).astype(int)

    onePassDelpwrs= np.zeros(len(lobes))
    raysInLobes = np.zeros(len(lobes))
    startingDelpwrs = np.zeros(len(lobes))
    try:
        delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
        radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))
        nparas = np.ma.getdata(cqlrf_nc.variables['wnpar'][:])
        for i in range(len(lobes)):
            lobe = lobes[i]

            if lobe > genray_in['grill']['ngrill']:
                continue
                
            for ray in range(len(delpwr)):
            
                if genray_in['grill'][f'anmax({lobe})'] >= nparas[ray][0] >= genray_in['grill'][f'anmin({lobe})']:
                    firstBounceIndex = findBounceIndex(radialVariable[ray],bounceToFind = 1)
                    
                    startingDelpwrs[i] += delpwr[ray][0]
                    onePassDelpwrs[i] += delpwr[ray][firstBounceIndex]
                    raysInLobes[i] +=1
        avgSPA = 1-(onePassDelpwrs/startingDelpwrs)
        avgSPA[avgSPA < 0] = np.nan
        avgSPA[avgSPA > 1] = np.nan

        return avgSPA, onePassDelpwrs, startingDelpwrs
    
    except:
        import traceback
        traceback.print_exc()
        print(f'something went wrong')
        return -1

def getAverageDampingNpara(targetDir = None, lobe = 1, effic = True):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()

    if not effic:
        print('need to implement')
        uh = oh

    cqlrf_nc = netCDF4.Dataset(f'{targetDir}/cql3d_krf001.nc','r')
    genray_nc = netCDF4.Dataset(f'{targetDir}/genray.nc','r')
    genray_in = getInputFileDictionary.getInputFileDictionary('genray_LH', targetDir=targetDir)

    delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
    
    Nparas = np.copy(genray_nc.variables["wnpar"]) #n_|| of the ray at each point along the ray trace
    
    #Npara_binEdges = np.linspace(0,1,100)
    #Npara_binCenters = (Npara_binEdges[:-1] + Npara_binEdges[1:])/2
    #Npara_depositions = np.zeros(len(Npara_binCenters))
    totalDeposition = 0
    weighted_effic = 0
    for ray in range(len(delpwr)):
        
        if genray_in['grill'][f'anmax({lobe})'] >= Nparas[ray][0] >= genray_in['grill'][f'anmin({lobe})']:
            ray_delpwr = delpwr[ray]
            ray_delpwr_normed = ray_delpwr/ray_delpwr[0]
            ray_npara = Nparas[ray]
            
            ray_npara = ray_npara[ray_delpwr_normed>0.05]
            ray_npara[np.isnan(ray_npara)] = 0
            ray_delpwr_normed = ray_delpwr_normed[ray_delpwr_normed>0.05]

            ray_npara_centers = np.abs(ray_npara[:-1] + ray_npara[1:])/2
            ray_pwrchange = ray_delpwr_normed[:-1] - ray_delpwr_normed[1:]
            
            effic = 1/ray_npara_centers**2
            
            weighted_effic += np.sum((effic) * ray_pwrchange)
            totalDeposition += np.sum(ray_pwrchange)


            #for j in range(len(powerChangeCenters)):
            #    if powerChangeCenters[j] > 0:
            #        index = findNearestIndex(Nparas_ray_center[j], Npara_binCenters)
            #        Npara_depositions[index] += powerChangeCenters[j]
       
    #Npara_binCenters = Npara_binCenters[Npara_depositions > 1e8]
    #Npara_depositions = Npara_depositions[Npara_depositions > 1e8]

    weightedAvg_Npara = weighted_effic/totalDeposition#np.sum(Npara_depositions*Npara_binCenters) / np.sum(Npara_depositions)

    #fig,ax = plt.subplots()
    #one = Npara_depositions
    #two = Npara_depositions/Npara_binCenters**2

    #weight_one = np.sum(one*Npara_binCenters) / np.sum(one)
    #weight_two = np.sum(two*Npara_binCenters) / np.sum(two)
    """
    ax.plot(Npara_binCenters, one/np.max(one))
    ax.plot(Npara_binCenters, two/np.max(two))
    ax.axvline(weight_one, lw = 2, color = 'k')
    ax.axvline(weight_two, lw = 2, color = 'k')
    ax.axvline(weightedAvg_Npara, lw = 2, color = 'k')
    plt.show()
    """

    return weightedAvg_Npara

#returns the SPA of the forward lobe
def getNPA(cqlrf_nc, genray_in,numBounces, lobes = [1]):
    avgNPA = 0
    startingDelpwr = -1
    numForwardLobeRays = -1
    try:
        delpwr= cqlrf_nc.variables["delpwr"][:] #power in the ray at each point
        radialVariable = (np.copy(cqlrf_nc.variables["spsi"]))
        nparas = np.ma.getdata(cqlrf_nc.variables['wnpar'][:])
        for ray in range(len(delpwr)):
            for lobe in lobes:
                if lobe > genray_in['grill']['ngrill']:
                    continue
                if genray_in['grill'][f'anmax({lobe})'] >= nparas[ray][0] >= genray_in['grill'][f'anmin({lobe})']:
                    firstBounceIndex = findBounceIndex(radialVariable[ray],bounceToFind = numBounces)
                    
                    if numForwardLobeRays == -1:
                        numForwardLobeRays = 0
                        startingDelpwr = 0
                    startingDelpwr +=delpwr[ray][0]
                    avgNPA += delpwr[ray][firstBounceIndex]
                    numForwardLobeRays +=1
        if startingDelpwr == -1:
            return np.nan
        avgNPA /= startingDelpwr
        avgNPA = 1-avgNPA

        return avgNPA
    except:
        import traceback
        traceback.print_exc()
        print(f'something went wrong')
        return -1

#returns the index of the array whose element is closest to value
def findNearestIndex(value, array):
    idx = (np.abs(array - value)).argmin()

    return idx

#returns the index at which the input ray bounces in the SOL
def findBounceIndex(radius, bounceToFind = 1):
    checkpointValue = 1
    pastCheckpoint = False

    numBounces = 0
    for i in range(len(radius)):
        if radius[i] <= checkpointValue:
            pastCheckpoint = True
        if radius[i] > checkpointValue and pastCheckpoint:
            numBounces += 1
            pastCheckpoint = False
            if numBounces == bounceToFind:
                return i
    return i

#returns the driven parallel current / RF power
def getParaCurrentDriveEffiencyApW(cql_nc):
    powerDepProf = cql_nc.variables["powrft"][-1]#W/cm^3
    fluxSurfaceVolume = cql_nc.variables["dvol"][:]#cm^3
    totalPower = np.sum(powerDepProf * fluxSurfaceVolume)
    
    ccurtor = cql_nc.variables["ccurtor"][-1,:]
    rya = cql_nc.variables["rya"][:]
    integTorCurrent = np.trapz(ccurtor, x = rya)

    curr = cql_nc.variables["curr"][-1,:]
    dArea = cql_nc.variables["darea"][:]
    paraCurrent = np.sum(curr * dArea)

    print(f'total parallel current: {paraCurrent:.4e}')

    return paraCurrent/totalPower

#returns the location of where the driven current is maximum
#this may or may not line up with the where the current density is highest
def getLocationOfCurrentPeak(cql_nc):
    import matplotlib.pyplot as plt
    curr = cql_nc.variables["curr"][-1,:]
    dArea = cql_nc.variables["darea"][:]
    paraCurrentProf = (curr * dArea)
    rya = cql_nc.variables["rya"][:]
    return rya[np.argmax(np.abs(paraCurrentProf))]

#returns the location of where the driven current density is maximum
#this may or may not line up with the where the current is highest
def getLocationOfCurrentDensityPeak(cql_nc):
    import matplotlib.pyplot as plt
    curr = cql_nc.variables["curr"][-1,:]
    dArea = cql_nc.variables["darea"][:]
    paraCurrentProf = (curr * dArea)
    rya = cql_nc.variables["rya"][:]
    return rya[np.argmax(np.abs(curr))]

#draw poloidal flux surfaces
def drawFluxSurfaces(ax, gfileDict = None, rhosToPlot = [.2,.4,.6,.8,1], 
                     colors = 'k', zBounds = None, limPath = None):
    if gfileDict == None:
        gfileDict = getGfileDict.getGfileDict()
    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    if zBounds == None:
        zBounds = [np.min(z), np.max(z)]
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 200)
    zInterp = np.linspace(zBounds[0],zBounds[1], 200)
    psirzNormInterp = interp2d(r,z, psirzNorm, kind = 'cubic')(rInterp, zInterp)

    if limPath is not None:
        Rs, Zs = np.meshgrid(rInterp, zInterp)
        points = np.vstack((Rs.ravel(), Zs.ravel())).T
        mask = limPath.contains_points(points).reshape(Rs.shape)
        psirzNormInterp[np.logical_not(mask)] = np.nan     

    ax.contour(rInterp, zInterp, psirzNormInterp, np.square(rhosToPlot), colors= colors, linewidths=2.25)

#returns the total current
#this I think also includes the inductive current when E|| isnt set to 0
def getCurrent(cql_nc):
    curr = cql_nc.variables["curr"][-1,:]*1e4/1e6#convert to MA/m^2
    darea = cql_nc.variables["darea"][:]/1e4#convert to m^2
    totalCD = np.sum(curr*darea)
    
    return totalCD#in MA

#returns rya and the absolute value of the safety factor profile
def getCQLq(targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cql_nc = netCDF4.Dataset(f'{targetDir}/cql3d.nc','r')
    rya = cql_nc.variables["rya"][:]
    q_prof = cql_nc.variables["qsafety"][:]

    return rya, np.abs(q_prof)

#returns the electron density profile in cql3d
def getCQLne(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    enescal = 1
    try:
        enescal = cqlInputDict['setup']['enescal']
    except:
        print(f'enescale not in input dict')
    n_e = cqlInputDict['setup']['enein(1,1)']*1e6*enescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        n_e = interp1d(rhos, n_e)(rho_pol)
        rhos = rho_pol

    return rhos, n_e

#returns the CQL3D Te profile
def getCQLTe(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    tescal = 1
    try:
        tescal = cqlInputDict['setup']['tescal']
    except:
        print(f'tescal not in input dict')
    T_e = cqlInputDict['setup']['tein']*tescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        T_e = interp1d(rhos, T_e)(rho_pol)
        rhos = rho_pol

    return rhos, T_e

#returns the deuterium density profile in cql3d
def getCQLnD(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    enescal = 1
    try:
        enescal = cqlInputDict['setup']['enescal']
    except:
        print(f'enescale not in input dict')
    n_D = cqlInputDict['setup']['enein(1,2)']*1e6*enescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        n_D = interp1d(rhos, n_D)(rho_pol)
        rhos = rho_pol

    return rhos, n_D

#returns the electron density profile in cql3d
def getCQLnC(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    enescal = 1
    try:
        enescal = cqlInputDict['setup']['enescal']
    except:
        print(f'enescale not in input dict')
    n_C = cqlInputDict['setup']['enein(1,3)']*1e6*enescal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        n_C = interp1d(rhos, n_C)(rho_pol)
        rhos = rho_pol

    return rhos, n_C

#returns the CQL3D Ti profile
def getCQLTD(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    tiscal = 1
    try:
        tiscal = cqlInputDict['setup']['tiscal']
    except:
        print(f'tiscal not in input dict')
    T_i = cqlInputDict['setup']['tiin']*tiscal
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        T_i = interp1d(rhos, T_i)(rho_pol)
        rhos = rho_pol

    return rhos, T_i

#returns the CQL3D Zeff profile
def getCQLZeff(targetDir = None, rho_pol = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    cqlInputDict = getInputFileDictionary.getInputFileDictionary('cql3d',targetDir = targetDir)
    Zeff = cqlInputDict['setup']['zeffin']
    rhos = cqlInputDict['setup']['ryain']

    if rho_pol is not None:
        Zeff = interp1d(rhos, Zeff)(rho_pol)
        rhos = rho_pol

    return rhos, Zeff

def getRhopolOfq(q, gfileDict):
    qpsi = np.abs(gfileDict["qpsi"])

    tryPsi = np.linspace(0,1,len(qpsi))
    closest_qIndex = findNearestIndex(q, qpsi)
    return np.sqrt(tryPsi)[closest_qIndex]

def convertRhopolToRhotor(rho_pol, targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()

    scratch = targetDir.split('/')[6]
    machine = scratch.split('_')[0]
    shotNum = scratch.split('_')[1]

    from omfit_classes import omfit_eqdsk
    eqdskName = shotToEqdsk.getEqdskName(shotNum, machine = machine)
    gfile = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')

    psi_n = gfile['fluxSurfaces']['geo']['psin']
    rho_pol = np.sqrt(psi_n)
    rho_n = gfile['fluxSurfaces']['geo']['rhon']

    rhopolToRhotor = interp1d(rho_pol,rho_n)

    return rhopolToRhotor(rho_pol)

def convertRhotorToRhopol(rho_tor, targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    
    scratch = targetDir.split('/')[6]
    machine = scratch.split('_')[0]
    shotNum = scratch.split('_')[1]

    from omfit_classes import omfit_eqdsk
    eqdskName = shotToEqdsk.getEqdskName(shotNum, machine = machine)
    gfile = omfit_eqdsk.OMFITgeqdsk(f'{targetDir}/{eqdskName}')

    psi_n = gfile['fluxSurfaces']['geo']['psin']
    rho_pol = np.sqrt(psi_n)
    rho_n = gfile['fluxSurfaces']['geo']['rhon']

    rhopolToRhotor = interp1d(rho_n,rho_pol)

    return rhopolToRhotor(rho_tor)


#convert the input rho_pol coordinate to major radius at the midplane
def convertRhopolToRmidplane(rhos, targetDir = None, side = 'LFS'):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)

    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 200)
    psirzNormFunc = interp2d(r,z, psirzNorm, kind = 'cubic')

    zmaxis = gfileDict["zmaxis"]
    Rmaxis = gfileDict['rmaxis']

    if side == 'LFS':
        relevantR_interp = rInterp[rInterp >= Rmaxis]
    elif side == 'HFS':
        relevantR_interp = rInterp[rInterp <= Rmaxis]
    else:
        relevantR_interp = rInterp
        print(f'returning both HFS and LFS values')

    midplanePsi = psirzNormFunc(relevantR_interp, zmaxis)
    
    Rs_fromRho_psi = interp1d(midplanePsi, relevantR_interp, fill_value = 'extrapolate', bounds_error = False)(rhos**2)
    return Rs_fromRho_psi

def getBmidplaneAtR(Rs,targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)

    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]

    btot = np.sqrt(gfileDict['btrz']**2 + gfileDict['brrz']**2 + gfileDict['bzrz']**2)

    Bfunc = interp2d(r,z, btot, kind = 'cubic')

    zmaxis = gfileDict["zmaxis"]
    midplaneB = Bfunc(Rs, zmaxis)

    return midplaneB

#convert the input major radius coordinate to rho_pol
def convertRtoRhopol(Rs, targetDir = None):
    if targetDir == None:
        targetDir = getTargetInfo.getTargetDir()
    gfileDict = getGfileDict.getGfileDict(targetDir=targetDir)
    r = gfileDict["rgrid"]
    z = gfileDict["zgrid"]
    psirz = gfileDict["psirz"]
    
    psi_mag_axis = gfileDict["ssimag"]
    psi_boundary = gfileDict["ssibdry"]
    
    psirzNorm = (psirz - psi_mag_axis)/(psi_boundary-psi_mag_axis)

    rInterp = np.linspace(np.min(r), np.max(r), 200)
    psirzNormFunc = interp2d(r,z, psirzNorm, kind = 'cubic')

    zmaxis = gfileDict["zmaxis"]
    Rmaxis = gfileDict['rmaxis']

    midplanePsi = psirzNormFunc(rInterp, zmaxis)
    
    rhopols_from_Rs = np.sqrt(interp1d(rInterp, midplanePsi, fill_value = 'extrapolate', bounds_error = False)(Rs))
    return rhopols_from_Rs
    
