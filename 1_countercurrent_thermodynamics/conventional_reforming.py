import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os

nCH4 = 1
nH2O = 1

feed_solution = ct.Solution('gri30.cti')
zeros = np.zeros(len(feed_solution.X))
feed_solution.X = zeros
feed_solution.X = {'CH4': 1, 'H2O': 1}
feedstock = ct.Mixture([(feed_solution, nH2O+nCH4)])



p = 1500000 # Pa

# make a heat map
# Make a heat map
nH2O = np.arange(1.0,5.0,0.05)
Trange = np.arange(573, 1183, 10)


Map_X_H2O = np.zeros((len(nH2O), len(Trange)))
Map_X_CH4 = np.zeros((len(nH2O), len(Trange)))
Map_Y_CO2 = np.zeros((len(nH2O), len(Trange)))
Map_Y_CO = np.zeros((len(nH2O), len(Trange)))

for i, n_H2O in enumerate(nH2O):
    for j, T in enumerate(Trange):
        feed_solution.X = {'CH4': 1, 'H2O': n_H2O}
        feed_solution.TP = T, p
        feedstock = ct.Mixture([(feed_solution, 1 + n_H2O)])
        feed_solution.equilibrate('TP')
        n_products = feedstock.phase_moles()[0]
        xCH4 = feed_solution.X[feed_solution.species_index('CH4')]
        xH2O = feed_solution.X[feed_solution.species_index('H2O')]
        xCO2 = feed_solution.X[feed_solution.species_index('CO2')]
        xCO = feed_solution.X[feed_solution.species_index('CO')]
        xH2 = feed_solution.X[feed_solution.species_index('H2')]
        # calculate conversion and yields
        X_CH4 = 1 - xCH4*n_products
        Y_CO2 = xCO2*n_products
        Y_CO = xCO*n_products
        X_H2O = 1 - xH2O * n_products / n_H2O
        Map_X_H2O[i][j] = X_H2O
        Map_X_CH4[i][j] = X_CH4
        Map_Y_CO[i][j] = Y_CO
        Map_Y_CO2[i][j] = Y_CO2



Trange, nH2O  = np.meshgrid(Trange, nH2O)



fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_X_CH4, rasterized=True, vmin = 0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97, 0.99]
cont = ax0.contour(Trange, nH2O, Map_X_CH4, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$X_\mathrm{CH_4}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","a_catalytic_CH4_conversion.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","a_catalytic_CH4_conversion.pdf"), bbox_inches= 'tight')
plt.show()

fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_X_H2O, rasterized=True, vmin =0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
cont = ax0.contour(Trange, nH2O, Map_X_H2O, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$X_\mathrm{H2O}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","a_catalytic_H2O_conversion.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","a_catalytic_H2O_conversion.pdf"), bbox_inches= 'tight')
plt.show()


fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_Y_CO, rasterized=True, vmin =0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5]
cont = ax0.contour(Trange, nH2O, Map_Y_CO, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$Y_\mathrm{CO}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","a_catalytic_CO_yield.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","a_catalytic_CO_yield.pdf"), bbox_inches= 'tight')
plt.show()

fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
ax0.set_xlabel('$T$ [K]')
ax0.set_ylabel('$\mathrm{H_2O/CH_4}$ [-]')

im = ax0.pcolormesh(Trange, nH2O, Map_Y_CO2, rasterized=True, vmin =0, vmax =1.0)
levels = [0.1, 0.2, 0.3, 0.4, 0.5]
cont = ax0.contour(Trange, nH2O, Map_Y_CO2, levels = levels, colors = 'black')
ax0.clabel(cont, fmt='%1.2f')
cb = fig.colorbar(im, ax=ax0)
cb.set_label(label ='$Y_\mathrm{CO2}$')
#plt.title(title)
#plt.plot(xCO2i, xCOi)
plt.savefig(os.path.join("plots","a_catalytic_CO2_yield.png"), dpi = 200, bbox_inches= 'tight')
plt.savefig(os.path.join("plots","a_catalytic_CO2_yield.pdf"), bbox_inches= 'tight')
plt.show()


fig = plt.figure(figsize=(4.5, 3.5), facecolor='white')
ax0 = plt.subplot()
