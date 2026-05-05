###
# Recreates Figure 5.4 from Andrea Schmidt's thesis "Measurements and Modeling of Lower Hybrid Driven Fast Electrons on Alcator C-Mod"
# However, the original had a /5 in the label. This plots things on the same scale.
# i got this data by digitizer her plot, so not all details may be perfectly recreated
###

import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 14)
plt.rc('legend', fontsize = 13)

fig,ax = plt.subplots()

measurements = np.loadtxt('/home/grantr/codes/GENRAY_CQL3D_scripts/thesisImages/andrea_meas.csv', delimiter=',')
gen = np.loadtxt('/home/grantr/codes/GENRAY_CQL3D_scripts/thesisImages/andrea_GEN.csv', delimiter=',')

ax.plot(gen[:,0], gen[:,1]*5, label = 'Predicted ECE Spectrum', lw = 3, color = 'tab:blue')
ax.plot(measurements[:,0], measurements[:,1], label = 'Measured ECE Spectrum', lw = 3, color = 'k')

ax.text(338,68, 'Nonthermal\nfeature', rotation = 0, fontsize = 14, ha='center', va='center')
ax.annotate(
            "",
            xytext=(338,62),
            xy=(244,56),
            arrowprops=dict( linewidth = 2, 
                            linestyle = '-', 
                            color = 'k',
                            arrowstyle = 'simple',
                            joinstyle='miter',   # sharp corners
                            capstyle='butt'  
                            )
            )

ax.legend(loc = 'best')
ax.set_ylabel(r'Radiation temperature T$_{rad}$ (keV)')
ax.set_xlabel(r'Frequency (GHz)')
ax.set_xlim([100,700])
ax.set_ylim(bottom=0)
ax.set_title(f'Alcator C-Mod shot 1060728011', loc = 'right')

fig.tight_layout()

plt.show()