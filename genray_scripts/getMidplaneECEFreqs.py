import numpy as np
import matplotlib.pyplot as plt
import getGfileDict

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 14)
plt.rc('legend', fontsize = 12)

gfile = getGfileDict.getGfileDict()

rgrid = gfile['rgrid']
Btmid = gfile['btmid'][rgrid >= 1.045]
Bzmid = gfile['bzmid'][rgrid >= 1.045]
rgrid = rgrid[rgrid >= 1.045]

#Brmid = gfile['brmid']

Btot = np.sqrt(Btmid**2 + Bzmid**2)

e = 1.602e-19
m_e = 9.109e-31

w = e*Btot/m_e

fig,ax = plt.subplots(figsize=(7,5))
#ax.plot(rgrid, w/(2*np.pi*1e9), lw = 3, label = r'$\Omega_e(E=0$ keV)', color = 'tab:blue')
ax.plot(rgrid[rgrid < 2.25], 2*w[rgrid < 2.25]/(2*np.pi*1e9), lw = 3, label = r'$2\Omega_e(\gamma=1$)', color = 'tab:green')
#ax.plot(rgrid, 3*w/(2*np.pi*1e9), lw = 3, label = r'$3\Omega_e(E=0$ keV)', color = 'tab:purple')

#ax.plot(rgrid, w/(2*np.pi*1e9*1.2), lw = 3, label = r'$\Omega_e(E=100$ keV)', linestyle = 'dashed', color = 'tab:blue')
ax.plot(rgrid[rgrid < 2.25], 2*w[rgrid < 2.25]/(2*np.pi*1e9*1.2), lw = 3, label = r'$2\Omega_e(\gamma=1.2$)', linestyle = 'dashed', color = 'tab:green')
#ax.plot(rgrid, 3*w/(2*np.pi*1e9*1.2), lw = 3, label = r'$3\Omega_e(E=100$ keV)', linestyle = 'dashed', color = 'tab:purple')
ax.set_ylim(bottom = 0)
ax.set_ylabel('Electron cyclotron frequency (GHz)')
ax.set_xlabel('Major radius (m)')
#ax.legend(ncol = 2)
ax.legend(ncol = 1)
fig.tight_layout()
plt.savefig('DIIID_147634.04565_ECEfreqs.jpeg',dpi=300)

plt.show()