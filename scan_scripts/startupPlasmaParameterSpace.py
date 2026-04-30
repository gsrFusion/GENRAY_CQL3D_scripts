import numpy as np
import matplotlib.pyplot as plt

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('legend', fontsize = 14)

plasmaParamsDict = {
    '187017': {3000: {'pinj': 5822.9504008016029, 'wplasma': 1023636.75, 'density': 44695617310930.133}, 4500: {'pinj': 5844.1242484969944, 'wplasma': 949066.58333333337, 'density': 62351612970240.32}, 1500: {'pinj': 5792.3752505010016, 'wplasma': 970512.83333333337, 'density': 37830057186616.367}},
    '186464': {2000: {'pinj': 5736.7730460921848, 'wplasma': 905766.66666666663, 'density': 65760803798728.172}, 3000: {'pinj': 5750.5430861723444, 'wplasma': 1001215.3333333334, 'density': 81908339701900.391}, 4500: {'pinj': 5751.4043086172342, 'wplasma': 890245.08333333337, 'density': 97957306978526.016}},
    '195188': {3300: {'pinj': 0.0, 'wplasma': 156196.20833333334, 'density': 21626198515121.871}, 1500: {'pinj': 1889.4814629258517, 'wplasma': 280474.625, 'density': 22494156588753.809}, 2750: {'pinj': 1894.872745490982, 'wplasma': 294420.25, 'density': 23208210833120.875}},
    '193592': {3000: {'pinj': 9356.5180360721442, 'wplasma': 818331.58333333337, 'density': 53151748312243.828}, 4000: {'pinj': 7516.6332665330665, 'wplasma': 768063.33333333337, 'density': 47139897666802.594}, 5000: {'pinj': 5638.3031062124246, 'wplasma': 661369.04166666663, 'density': 46178812464020.344}, 1500: {'pinj': 5417.8316633266531, 'wplasma': 609037.54166666663, 'density': 49846492369439.883}},
    '193593': {3000: {'pinj': 4539.137274549098, 'wplasma': 380740.0625, 'density': 29184615881973.953}, 4000: {'pinj': 3573.0187875751503, 'wplasma': 368953.0625, 'density': 30640800020681.254}, 5000: {'pinj': 1673.4203406813626, 'wplasma': 191776.64583333334, 'density': 27439183327310.387}, 1500: {'pinj': 5420.8512024048096, 'wplasma': 429350.66666666669, 'density': 38075120772619.453}},
    '194547': {1400: {'pinj': 2186.1154809619238, 'wplasma': 214204.79166666666, 'density': 17956524044653.625}, 4500: {'pinj': 4265.3451903807618, 'wplasma': 698484.33333333337, 'density': 45013680672094.648}, 2500: {'pinj': 4228.7179358717431, 'wplasma': 724817.41666666663, 'density': 47771054983607.641}, 5750: {'pinj': 2146.1548096192387, 'wplasma': 592251.0, 'density': 53793236436820.75}},
    '194550': {1400: {'pinj': 2186.1154809619238, 'wplasma': 214204.79166666666, 'density': 17956524044653.625}, 4500: {'pinj': 4265.3451903807618, 'wplasma': 698484.33333333337, 'density': 45013680672094.648}, 2500: {'pinj': 4228.7179358717431, 'wplasma': 724817.41666666663, 'density': 47771054983607.641}, 5750: {'pinj': 2146.1548096192387, 'wplasma': 592251.0, 'density': 53793236436820.75}},
    '194749': {2500: {'pinj': 5899.419338677355, 'wplasma': 452015.95833333331, 'density': 25049689392113.664}, 3650: {'pinj': 4711.2479959919838, 'wplasma': 743803.8125, 'density': 48494862072438.961}, 1300: {'pinj': 3715.9286072144287, 'wplasma': 239761.54166666666, 'density': 21061060872133.164}, 4900: {'pinj': 5956.6087174348695, 'wplasma': 462772.79166666669, 'density': 24825550499956.129}},
    '198760': {1200: {'pinj': 6281.6720901126409, 'wplasma': 548999.20833333337, 'density': 24047790417807.953}, 3500: {'pinj': 6327.6013767209015, 'wplasma': 597108.375, 'density': 39831260104805.641}, 4500: {'pinj': 3847.9834167709637, 'wplasma': 329329.125, 'density': 28957001373771.125}, 2500: {'pinj': 8256.0763454317894, 'wplasma': 788531.83333333337, 'density': 39679703586524.688}, 1750: {'pinj': 7662.0275344180227, 'wplasma': 595368.25, 'density': 28449422920929.227}},
    '198762': {2500: {'pinj': 6328.5068836045057, 'wplasma': 648396.25, 'density': 35588785537483.469}, 4500: {'pinj': 3861.6974342928661, 'wplasma': 522694.91666666669, 'density': 39633464432723.922}, 3500: {'pinj': 6338.489361702128, 'wplasma': 697173.25, 'density': 37530706980902.016}, 1500: {'pinj': 7812.7878598247808, 'wplasma': 875702.08333333337, 'density': 44519347197957.469}},
    '198868': {4000: {'pinj': 5391.8692115143931, 'wplasma': 627115.29166666663, 'density': 68055141528206.164}, 5500: {'pinj': 7422.442428035044, 'wplasma': 568175.5, 'density': 53018288086060.484}, 1500: {'pinj': 5357.984355444305, 'wplasma': 463021.33333333331, 'density': 30587454731696.664}, 2750: {'pinj': 5391.2008760951185, 'wplasma': 507350.96875, 'density': 49063986309648.766}}

}

dampingDict = {
     '187017': {3000: True, 4500: True, 1500: True},
    '186464': {2000: True, 3000: True, 4500: True},
    '195188': {3300: False, 1500:True, 2750: True},
    '193592': {3000:True, 5000: True, 1500: True},
    '193593': {3000: True, 4000:True, 5000:False , 1500:False},
    '194547': {1400: True, 4500: True,  2500: True, 5750: False},
    '194550': {1400: True, 4500: True, 2500: True, 5750: True},
    '194749': {2500:True , 3650: True, 1300: True, 4900:True },
    '198760': {1200: False, 3500: True, 4500:  False, 2500: True, 1750: True},
    '198762': {2500:True , 4500: True, 3500:True, 1500:True},
    '198868': {4000: True, 5500: True, 1500: True, 2750: True}

}

def densityPinjPlot():
    fig,ax = plt.subplots()
    goodParams = np.empty((0, 3))
    badParams = np.empty((0, 3))
    for shot in plasmaParamsDict.keys():
        damping = dampingDict[shot]
        plasmaParams = plasmaParamsDict[shot]

        print(f'{shot} ')
        print(f'{damping} \n')

        times = damping.keys()
        for time in times:
                paramArray = np.array([plasmaParams[time]['density'],plasmaParams[time]['pinj'],plasmaParams[time]['wplasma']])
                if damping[time]:
                    goodParams = np.append(goodParams, [paramArray], axis = 0)
                else:
                    badParams = np.append(badParams, [paramArray], axis = 0)

    ax.scatter(goodParams[:,0]/1e13,  goodParams[:,1]/1e3, label = 'good targets')
    ax.scatter(badParams[:,0]/1e13,  badParams[:,1]/1e3, label = 'worse targets')
    ax.legend(loc = 'best')
    ax.set_ylabel(r'P$_{beam}$ (MW)')
    ax.set_xlabel(r'$\bar{n_e}$ ($10^{19}$ $m^{-3}$)')
    fig.tight_layout()
    plt.show()

def wplasmaPlot():
    fig,ax = plt.subplots()
    goodWplasma = np.empty(1)
    badWplasma = np.empty(1)

    for shot in plasmaParamsDict.keys():
        damping = dampingDict[shot]
        plasmaParams = plasmaParamsDict[shot]

        print(f'{shot} ')
        print(f'{damping} \n')

        times = damping.keys()
        for time in times:
                if damping[time]:
                    goodWplasma = np.append(goodWplasma, [plasmaParams[time]['wplasma']], axis = 0)
                else:
                    badWplasma = np.append(badWplasma, [plasmaParams[time]['wplasma']], axis = 0)

    ax.scatter([1]*len(goodWplasma), goodWplasma, label = 'good targets')
    ax.scatter([-1]*len(badWplasma), badWplasma, label = 'worse targets')
    ax.legend(loc = 'best')
    ax.set_ylabel(r'W$_{plasma}$ (J)')
    ax.set_xlabel(r'')
    fig.tight_layout()
    plt.show()

wplasmaPlot()
