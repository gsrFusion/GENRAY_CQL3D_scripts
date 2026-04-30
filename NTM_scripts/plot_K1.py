import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import RectBivariateSpline as RBS

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('figure', titlesize = 22)
plt.rc('legend', fontsize = 14)

fig,ax = plt.subplots()

with open('/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/K_1Dict.pkl', 'rb') as f:
    K_1Dict = pickle.load(f)

widthPoints = K_1Dict['widthPoints']
alignmentPoints = K_1Dict['alignmentPoints']
K_1_width_alignment = K_1Dict['K_1_width_alignment']
#"""
p2 = ax.contourf(alignmentPoints,widthPoints, K_1_width_alignment, levels = [-.2,-.1,0,.1,.2,.3,.4,.5,.6,.7])
#p2 = ax.contourf(alignmentPoints, widthPoints, K_1_width_alignment,levels = 10)
cbar = fig.colorbar(p2, ax = ax, pad = .01)
cbar.set_label(r"$K_1$")



print(f'{K_1_width_alignment.shape, alignmentPoints.shape, widthPoints.shape}')
widthMask = np.where((widthPoints > 1)*(widthPoints < 2.5))[0]
alignmentMask = np.where(alignmentPoints<.1)[0]
print(alignmentPoints[alignmentMask])
toTest = K_1_width_alignment[widthMask,:]
toTest = toTest[:,alignmentMask]
print(f'{toTest.shape}')
print(f'max in test, min in test: {np.max(toTest), np.min(toTest)}')

ax.set_ylabel(r'Relative width $w/\delta_{CD}$')
ax.set_xlabel(r'Relative alignment $\Delta R/\delta_{CD}$')
fig.tight_layout()
plt.show()
#"""
fig,ax = plt.subplots()
interp_spline = RBS(widthPoints, alignmentPoints, K_1_width_alignment)
print(f'after RBS')
K_1_interp = interp_spline(widthPoints, alignmentPoints)

p2 = ax.contourf(alignmentPoints,widthPoints, K_1_width_alignment - K_1_interp, levels = 10)
#p2 = ax.contourf(alignmentPoints, widthPoints, K_1_width_alignment,levels = 10)
cbar = fig.colorbar(p2, ax = ax, pad = .01)
cbar.set_label(r"$K_1$")


ax.set_ylabel(r'Relative width $w/\delta_{CD}$')
ax.set_xlabel(r'Relative alignment $\Delta R/\delta_{CD}$')
fig.tight_layout()
plt.show()

fig,ax = plt.subplots()
ax.plot(alignmentPoints, interp_spline([1.5], alignmentPoints)[0])
plt.show()