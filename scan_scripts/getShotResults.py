import numpy as np
import os, sys
import netCDF4
import matplotlib.pyplot as plt
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import helperFunctions as helper
import getInputFileDictionary


shotDict = {}

def garbageCollection():
    rootDir1 = f'/home/grantr/symlinks/c7_scratch/genray_batch/DIIID_shots/'
    folderList1 = fast_scandir(rootDir1)
    rootDir2 = f'/home/grantr/orcd/c7/pool/'
    folderList2 = fast_scandir(rootDir2)
    
    totalList = np.concatenate([folderList1, folderList2])

    fig, ax = plt.subplots()
    illegalWords = ['profiles','thgrill','mod','tribot','skewed','scal', 'gap','enorm','rev','edge','error','field','helicon','twocolor']
    SPA_crit = .5

    for folder in totalList:
        
        try:
            if not ('2.7' in folder):
                continue
            if any(word in folder.lower() for word in illegalWords):
                continue
            print(folder)
            cqlrf_nc = netCDF4.Dataset(f'{folder}/cql3d_krf001.nc','r')
            genray_in = getInputFileDictionary.getInputFileDictionary('genray', targetDir = folder)
            SPA = helper.getSPA(cqlrf_nc, genray_in, lobes = [1])
            
            shotAndTime = folder.split('/')[-1].split('_')[1]
            if len(shotAndTime.split('.')[0])!=6 or len(shotAndTime.split('.')[1])!=5:
                continue
            isGood = False
            if SPA > SPA_crit:
                isGood = True
            shotDict[shotAndTime] = isGood
        except:
            pass
    print(shotDict)

def fast_scandir(dirname):
    subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
    for dirname in list(subfolders):
        subfolders.extend(fast_scandir(dirname))
    return subfolders  


garbageCollection()
