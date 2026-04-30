import numpy as np
import matplotlib.pyplot as plt
from DetectorEfficiency import getEfficiencies


peakEnergies = [121,344,1408,964,1112,778,1085,244,867,443,411,1089,1299,1212]
peakIntensities = [28.58,26.5,21,14.61,13.644,12.942,10.207,7.583,4.245,2.821,2.234,1.727,1.623,1.422]

fig,ax = plt.subplots()

ax.scatter(peakEnergies, peakIntensities/np.max(peakIntensities), label = 'relative emitted intensities',s=8)


measuredIntensities = getEfficiencies(peakEnergies)*peakIntensities
relativeMeasuredIntensities = measuredIntensities/np.max(measuredIntensities)
ax.scatter(peakEnergies, relativeMeasuredIntensities, label = 'relative measured intensities',s=8)

print(np.sum(relativeMeasuredIntensities[3:]))

ax.legend()

plt.show()
