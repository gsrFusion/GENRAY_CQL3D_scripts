import matplotlib.pyplot as plt
import numpy as np

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 22)
plt.rc('legend', fontsize = 14)

def square(x, start, w):
    result = np.ones(len(x))
    result[x < start] = 0
    result[x > start + w] = 0
    return result

countRate = 500
allowedVariation = .5
minTimeBin = 1e3/(allowedVariation**2 * countRate)#in ms

time = np.linspace(0, 250, 1000)

MSEonTime = 15#ms
neutronFalloffTime = 20#ms
print(f'minTimeBin: {minTimeBin}')
fig,ax = plt.subplots()

t = 0
while t < max(time):
    beam = square(time, t, MSEonTime)
    HXR_on = square(time, t + MSEonTime + neutronFalloffTime, minTimeBin)
    HXR_off = square(time, t + MSEonTime + neutronFalloffTime+minTimeBin, minTimeBin)
    neutrons = square(time, t + MSEonTime,neutronFalloffTime)
    t = t + MSEonTime + neutronFalloffTime + minTimeBin*2
    
    ax.fill_between(time, 0,beam, color = 'tab:blue')#, label = 'MSE')
    ax.fill_between(time, 0,HXR_on, color = 'tab:red')#, label = 'HXR')
    ax.fill_between(time, 0,HXR_off, color = 'darkred')#, label = 'HXR')
    ax.fill_between(time, 0,neutrons, color = 'tab:purple')#, label = 'HXR')

ax.plot(time, beam+1e3, color = 'tab:blue', lw = 4, label = 'MSE')
ax.plot(time, beam+1e3, color = 'tab:purple', lw = 4, label = 'neutron falloff')
ax.plot(time, beam+1e3, color = 'tab:red', lw = 4, label = 'HXR, LH off')
ax.plot(time, beam+1e3, color = 'darkred', lw = 4, label = 'HXR, LH on')

ax.set_ylim([0,1.2])
ax.set_ylabel('')
ax.set_xlabel('time (ms)')
ax.legend(ncol = 2, loc = 'best')
fig.tight_layout()

plt.show()