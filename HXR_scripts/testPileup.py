import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

pulseFWHM = .4e-6#.475e-6 are John's current cards

digitizingTime = 1#s
averagingError = .01
criterionFactor = 1.1


plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 18)
plt.rc('legend', fontsize = 14)


#returns the times at which counts are detected over the course of 1 second
def getDetectionTimes(countRate):
    pulseSigma = pulseFWHM/2.355

    std = np.sqrt(countRate)
    instantaneousCountRates = np.random.normal(loc = countRate, scale = std, size = int(1.5*countRate*digitizingTime))

    scales = 1/instantaneousCountRates
    #print(instantaneousCountRates)
    timesBetweenPulses = np.random.exponential(scale = scales)

    detectionTimes = np.cumsum(timesBetweenPulses)
    detectionTimes = detectionTimes[detectionTimes < digitizingTime]

    return detectionTimes

def getNumPileupFromDetectionTimes(detectionTimes):
    diffTimes = np.diff(detectionTimes)
    numPileup = len(diffTimes[diffTimes < pileupCriterion])
    return numPileup

def plotDetections(detectionTimes):
    pulseSigma = pulseFWHM/2.355
    fig,ax = plt.subplots()
    t = np.linspace(0,1,int(25e6))
    trace = np.zeros(len(t))
    for time in detectionTimes:
        print('here')
        #print(f'time: {time}, traceAddition: {np.exp(-.5*(t - time)**2/pulseSigma**2)}')
        trace += np.exp(-.5*(t - time)**2/pulseSigma**2)

    ax.plot(t, trace)
    plt.show()

def getAveragePileup(countRate):
    minSamples = 10

    totalSamples = minSamples
    totalPileup = 0

    previousAvgPileup = 1e-20

    for i in range(minSamples):
        totalPileup += getNumPileupFromDetectionTimes(getDetectionTimes(countRate))
    
    while abs(previousAvgPileup - totalPileup/totalSamples)/previousAvgPileup > averagingError:
        previousAvgPileup = totalPileup/totalSamples
        for i in range(minSamples):
            totalPileup += getNumPileupFromDetectionTimes(getDetectionTimes(countRate))
            totalSamples += 1

    print(f'totalSamples: {totalSamples} required')

    return totalPileup/totalSamples

countRates = np.linspace(5e4,7e5,20)
fwhms = np.array([.12,.25,.325,.475])*1e-6#[
fig,ax = plt.subplots()
for l in range(len(fwhms)):
    print(f'fwhm: {fwhms[l]}')
    pileupCriterion = criterionFactor*fwhms[l]
    portionPileup = np.zeros(len(countRates))
    pulseFWHM = fwhms[l]
    for k in range(len(countRates)):
        avgTotalCounts = countRates[k]*digitizingTime
        #plotDetections(getDetectionTimes(countRate))
        portionPileup[k] = getAveragePileup(countRates[k])/avgTotalCounts


    ax.plot(countRates, portionPileup, linewidth = 2, label = f'FWHM = {fwhms[l]*1e6: .2f} us')
ax.set_ylabel('Portion of pulses piled up')
ax.set_xlabel('Count rate (Hz)')
ax.axhline(.1, color ='k', linestyle = 'dashed')
ax.set_xscale('log')
ax.set_title(f'criterionFactor: {criterionFactor}')
ax.legend()
fig.tight_layout()
plt.show()