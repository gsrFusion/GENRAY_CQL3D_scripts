import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d

isNTM = True

fluxSurfaces = [[],[],[],[],[],[],[],[],[],[],[],[]]#[LCFS, rho_8, rho_6_normal, rho_4, rho_2, rho_2_normal]
csvLabels = ['LCFS', '.8','.6_normal','.6_NTM1','.6_NTM2','.6_NTM3','.6_NTM4','.6_NTM5','.6_NTM6','.4','.2_normal','.2',]

with open('/home/grantr/codes/GENRAY_CQL3D_scripts/NTM_scripts/dataStorage/NTM_figInfo.csv', mode ='r') as file:
  csvFile = csv.DictReader(file)
  for line in csvFile:
        for i in range(len(fluxSurfaces)):
            try: 
                fluxSurfaces[i].append((float(line[csvLabels[i]+'_X']), float(line[csvLabels[i]+'_Y'])))
            except:
                pass

def trickySort(points):
    # Start from the leftmost point
    start_idx = np.argmin(points[:, 0])
    order = [start_idx]
    unused = set(range(len(points))) - {start_idx}

    while unused:
        last = points[order[-1]]
        # Find nearest unused point
        next_idx = min(unused, key=lambda i: np.linalg.norm(points[i] - last))
        order.append(next_idx)
        unused.remove(next_idx)

    # Apply ordering
    points_sorted = points[order]
    return points_sorted

for i in range(len(fluxSurfaces)):
    fluxSurfaces[i] = np.array(fluxSurfaces[i])

    if not ('NTM' in csvLabels[i]):
        fluxSurfaces[i] = fluxSurfaces[i][np.argsort(np.arctan2(fluxSurfaces[i][:,1]-np.average(fluxSurfaces[i][:,1]),fluxSurfaces[i][:,0]-np.average(fluxSurfaces[i][:,0])))]
        fluxSurfaces[i] = np.append(fluxSurfaces[i], [fluxSurfaces[i][0]], axis= 0)

    else:
        fluxSurfaces[i] = trickySort(fluxSurfaces[i])
        fluxSurfaces[i] = np.append(fluxSurfaces[i], [fluxSurfaces[i][0]], axis= 0)

    indices = np.arange(len(fluxSurfaces[i]))
    newIndices = np.arange(1000)*np.max(indices)/1000
    Xs = fluxSurfaces[i][:,0]
    Ys = fluxSurfaces[i][:,1]

    newXs = interp1d(indices, Xs)(newIndices)
    newYs = interp1d(indices, Ys)(newIndices)

    newFluxSurface = np.array([newXs,newYs]).T
    fluxSurfaces[i] = newFluxSurface

    from scipy.signal import savgol_filter

    # window_length must be odd and <= len(points)
    numBuffer = 50

    X_toSmooth = np.concatenate([fluxSurfaces[i][:,0][-(numBuffer+1):-1], fluxSurfaces[i][:,0], fluxSurfaces[i][:,0][1:(numBuffer + 1)]])
    Y_toSmooth = np.concatenate([fluxSurfaces[i][:,1][-(numBuffer+1):-1], fluxSurfaces[i][:,1], fluxSurfaces[i][:,1][1:(numBuffer + 1)]])

    smoothed_X = savgol_filter(X_toSmooth, window_length=201, polyorder=3)
    smoothed_Y = savgol_filter(Y_toSmooth, window_length=201, polyorder=3)

    fluxSurfaces[i][:,0] = smoothed_X[numBuffer:-numBuffer]
    fluxSurfaces[i][:,1] = smoothed_Y[numBuffer:-numBuffer]

    fluxSurfaces[i][-1] = np.copy(fluxSurfaces[i][0])

lw = 2.5

fig,ax = plt.subplots()

ax.plot(fluxSurfaces[0][:,0], fluxSurfaces[0][:,1], lw = lw, color = 'k',solid_capstyle='butt')
ax.plot(fluxSurfaces[1][:,0], fluxSurfaces[1][:,1], lw = lw, color = 'k',solid_capstyle='butt')

ax.plot(fluxSurfaces[9][:,0], fluxSurfaces[9][:,1], lw = lw, color = 'k',solid_capstyle='butt')
ax.plot(fluxSurfaces[11][:,0], fluxSurfaces[11][:,1], lw = lw, color = 'k',solid_capstyle='butt')

if isNTM:
    ax.plot(fluxSurfaces[3][:,0], fluxSurfaces[3][:,1], lw = lw, color = 'k',solid_capstyle='butt')
    ax.plot(fluxSurfaces[4][:,0], fluxSurfaces[4][:,1], lw = lw, color = 'k',solid_capstyle='butt')
    ax.plot(fluxSurfaces[5][:,0], fluxSurfaces[5][:,1], lw = lw, color = 'k',solid_capstyle='butt')
    ax.plot(fluxSurfaces[6][:,0], fluxSurfaces[6][:,1], lw = lw, color = 'k',solid_capstyle='butt')
    ax.plot(fluxSurfaces[7][:,0], fluxSurfaces[7][:,1], lw = lw, color = 'k',solid_capstyle='butt')
    ax.plot(fluxSurfaces[8][:,0], fluxSurfaces[8][:,1], lw = lw, color = 'k',solid_capstyle='butt')
else:
    avg2_normalX = np.average(fluxSurfaces[-2][:,0])
    avg2_X = np.average(fluxSurfaces[-1][:,0])
    diffX = avg2_X - avg2_normalX
    ax.plot(fluxSurfaces[2][:,0] + diffX, fluxSurfaces[2][:,1], lw = lw, color = 'k',solid_capstyle='butt')

ax.set_aspect('equal')

plt.show()