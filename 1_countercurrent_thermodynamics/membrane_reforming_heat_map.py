import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os

def countercurrent_reforming(T, p, nH2O):
    """Determine the equilibrium conversion limits of CH4 and H2O with 1CH4:nH2O flow rates in a countercurrent
    membrane reactor. Reactions
    flow 1: H2O -> H2 + 0.5 O2
    flow 2: CH4 + k O2 -> x_1 CH4 + x_2 CO + x_3 CO2 + x_4 H2 + x_5 H2O
    Returns, H2O conversion, CH4 conversion and product yields in flow_1"""
    # define the two flows as canter solution objects using cantera's gri30 database
    flow1 = ct.Solution('gri30.cti')
    flow2 = ct.Solution('gri30.cti')
    zeros = np.zeros(len(flow1.X))
    # Make arrays tof the pO2(kappa) in each flow
    PO2_flow2 = []
    PO2_flow1 = []
    kappa_max = 0
    # complete conversion of CO2 to CO gives Kappa = 0.5, kappa in range 0-0.5
    kappa = np.arange(0, nH2O*0.5+0.01*nH2O*0.5, 0.01*nH2O*0.5)
    for k in kappa:
        flow1.X = zeros
        flow2.X = zeros
        flow1.X = {'H2': nH2O, 'O2': nH2O*0.5 - k}
        flow2.X = {'CH4': 1, 'O2': k}
        flow1.TP = T, p
        flow2.TP = T, p
        flow1.equilibrate('TP')
        flow2.equilibrate('TP')
        pO2_flow1 = flow1.X[flow1.species_index('O2')] * p
        pO2_flow2 = flow2.X[flow2.species_index('O2')] * p
        # add the new pO2 values to the arrays
        PO2_flow1.append(pO2_flow1)
        # for the second flow add pO2 value to the start of array to reverse kappa for countercurrent
        PO2_flow2.insert(0, pO2_flow2)
        # if the arrays meet at any point, then we have reached the maximum exchange extent kappa_max
        if np.any(np.asarray(PO2_flow1)-np.asarray(PO2_flow2) <= 0) or k > 0.999*nH2O*0.5 or k > 1.999:
            kappa_max = k
            xCH4 = flow2.X[flow2.species_index('CH4')]
            xH2O = flow2.X[flow2.species_index('H2O')]
            xCO2 = flow2.X[flow2.species_index('CO2')]
            xCO = flow2.X[flow2.species_index('CO')]
            xH2 = flow2.X[flow2.species_index('H2')]
            # calculate conversion and yields
            X_CH4 = 1 - xCH4 / (xCH4 + xCO2 + xCO)
            Y_CO2 = xCO2 / (xCH4 + xCO2 + xCO)
            Y_CO = xCO / (xCH4 + xCO2 + xCO)
            Y_H2O = xH2O/(xH2+xH2O+2*xCH4)
            Y_H2 = xH2/(xH2+xH2O+2*xCH4)
            X_H2O = kappa_max*2/nH2O
            break
    return X_H2O, X_CH4, Y_CO2, Y_CO, Y_H2O, Y_H2


p = 1500000 # Pa

# make a heat map
# Make a heat map
nH2O = np.arange(1.0,5.0,0.025)
Trange = np.arange(573, 1183, 5)
"""
# the calculation is slow so I saved the results for plotting and load it from a file

Map_X_H2O = np.zeros((len(nH2O), len(Trange)))
Map_X_CH4 = np.zeros((len(nH2O), len(Trange)))
Map_Y_CO2 = np.zeros((len(nH2O), len(Trange)))
Map_Y_H2O = np.zeros((len(nH2O), len(Trange)))

for i, n_H2O in enumerate(nH2O):
    print("Progress: ", round(i*100 / len(nH2O),3), " %")
    for j, T in enumerate(Trange):
        X_H2O, X_CH4, Y_CO2, Y_CO, Y_H2O, Y_H2 = countercurrent_reforming(T, p, n_H2O)
        Map_X_H2O[i][j] = X_H2O
        Map_X_CH4[i][j] = X_CH4
        Map_Y_CO2[i][j] = Y_CO2
        Map_Y_H2O[i][j] = Y_H2O


np.savetxt('CT_reforming_X_H2O_map.csv', Map_X_H2O, delimiter=",")
np.savetxt('CT_reforming_X_CH4_map.csv', Map_X_CH4, delimiter=",")
np.savetxt('CT_reforming_Y_CO2_map.csv', Map_Y_CO2, delimiter=",")
np.savetxt('CT_reforming_Y_H2O_map.csv', Map_Y_H2O, delimiter=",")
"""

Map_X_H2O = np.loadtxt("CT_reforming_X_H2O_map.csv", delimiter=',')
Map_X_CH4 = np.loadtxt("CT_reforming_X_CH4_map.csv", delimiter=',')
Map_Y_CO2 = np.loadtxt("CT_reforming_Y_CO2_map.csv", delimiter=',')
Map_Y_H2O = np.loadtxt("CT_reforming_Y_H2O_map.csv", delimiter=',')


Trange, nH2O  = np.meshgrid(Trange, nH2O)



fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_X_CH4, rasterized=True, vmin =0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97, 0.99]
cont = ax0.contour(Trange, nH2O, Map_X_CH4, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$X_\mathrm{CH_4}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","CH4_conversion.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","CH4_conversion.pdf"), bbox_inches= 'tight')
plt.show()

fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_X_H2O, rasterized=True, vmin =0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97]
cont = ax0.contour(Trange, nH2O, Map_X_H2O, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$X_\mathrm{H2O}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","H2O_conversion.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","H2O_conversion.pdf"), bbox_inches= 'tight')
plt.show()


fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_Y_H2O, rasterized=True, vmin =0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
cont = ax0.contour(Trange, nH2O, Map_Y_H2O, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$Y_\mathrm{H2O}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","H2O_yield.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","H2O_yield.pdf"), bbox_inches= 'tight')
plt.show()

fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_Y_CO2, rasterized=True, vmin =0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
cont = ax0.contour(Trange, nH2O, Map_Y_CO2, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$Y_\mathrm{CO2}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","CO2_yield.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","CO2_yield.pdf"), bbox_inches= 'tight')
plt.show()

