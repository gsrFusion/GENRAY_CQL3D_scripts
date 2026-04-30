import numpy as np
import matplotlib.pyplot as plt
import getTargetInfo

targetDir = getTargetInfo.getTargetDir()
shotNum = getTargetInfo.getShotNum()
machine = getTargetInfo.getMachine()

import getGfileDict
gfileDict = getGfileDict.getGfileDict(targetDir = targetDir)

current = gfileDict['cpasma']
B0 = gfileDict['bcentr']

Jsign = np.sign(current)
Bsign = np.sign(B0)

print(f'current: {current}, B0: {B0}')

if Bsign == Jsign:
    print(f'To drive co-current, Npara should be negative')
else:
    print(f'To drive co-current, Npara should be positive')
