machine = 'NTPT'

if machine == 'FENIX':
    shot = 'PTXPT'
    timeString = ''

    intermediaryDir = 'npara_thgrill_scan/'

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_n2.7Npara_170thgrill_1MW_newEqdsk'

if machine == 'DIIID':
    shot ='180403'
    timeString = '.04400'

    intermediaryDir = ''

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_n2.9Npara_1MW_id2'
    
    #"""
    shot ='203912'
    timeString = '.02700'

    intermediaryDir = ''

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_expSpectrum_first'
    #"""
    """
    shot ='180403'
    timeString = '.04400'

    intermediaryDir = ''

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_n2.9Npara_1MW_id2'
    """
    """
    shot ='206629'
    timeString = '.01980'

    intermediaryDir = ''

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_ECCD'
    """
    """
    shot ='147634'
    timeString = '.04565'

    intermediaryDir = 'Npara_height_scan/'

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_n2.9Npara_0.0grillHeight_1MW_LFS'
    """

if machine == 'KSTAR':
    shot ='39608'
    timeString = '.009250'

    intermediaryDir = ''#f'ECE_testing/'#3modules/'

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{timeString}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{timeString}_n2.9Npara_220thgrill_1MW'

if machine == 'MANTA':
    shot = 'posCS'

    intermediaryDir = f'npara_thgrill_scan/'

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}_n2.25Npara_170thgrill_1MW'

if machine == 'NTPT':
    shot = 'ARC'
    tri = '.V3APT'

    intermediaryDir = 'LFSVersion/'

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{tri}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{tri}_n1.9Npara_0.25grillHeight_10MW'
    #"""
    shot = 'DIIID'
    tri = '.147634PT'

    intermediaryDir = ''

    topmostShotDir = f'/home/grantr/symlinks/genray_batch/{machine}_shots/{machine}_{shot}{tri}'
    targetDir = f'{topmostShotDir}/{intermediaryDir}{machine}_{shot}{tri}_n2.8Npara_0.0grillHeight_1MW'
    #"""



#returns the directory where the shot of interest is stored
def getTargetDir():
    return targetDir.strip()

#returns the machine where the shot is from (DIIID, WEST, etc)
def getMachine():
    return topmostShotDir.split('/')[-1].split('_')[-2]

def getShotNum():
    return topmostShotDir.split('_')[-1]

def getTopmostShotDir():
    return topmostShotDir