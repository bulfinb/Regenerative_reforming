import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import integrate
import os

"""In this file we analyse the experimental data from the comparison PBR experiment shown in figure 4 of the manuscript."""


# RELATIVE ERROR IN THE MASS BALANCE
# 2 % error in bottle concentration. 2 % error in gas analysis signals, and 1 % error in the mass flow.
# The errors in the gas-analysis are used twice in each mass balance.
# This neglects the numerical errors introduced during integration of the signals.
relative_error_MB = np.sqrt(2 * (0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))
print("relative error mass balance:", relative_error_MB)

# Bottle concentrations and molar volumes
x_CO2_Bottle = 0.0501           # rest Ar => Accuracy 2%
x_CH4_Bottle = 0.0496            # rest Ar => Accuracy 2%
V_0 = 22.414                    # Molar Volume [dm^3/mol] = [l/mol]


# Import Data
filename = "24052023_Demo_930C_dataset.csv"
collumn_names= ['time', 'T', 'p', 'v_CH4', 'v_CO2', 'v_Ar', 'x_H2', 'x_CO2', 'x_CO','x_O2', 'x_CH4', 'Ar MFC sp', 'CH4 MFC sp', 'CO2 MFC sp']
# Units
# time  s
# T     C
# p     bar (above atmospheric, i.e. p = p_amb - p_reactor)
# v_i   standard litres/min
# x_i   %
# MFC sp  - MFC set point standard litres/min

#  load data as pandas arrays and make an instance of the Tga_data object for each file
df = pd.read_csv( filename, skiprows=1, sep=',', names=collumn_names)

# convert measured % to mole fractions
df['x_H2'] = df['x_H2'] * 0.01
df['x_CO'] = df['x_CO'] * 0.01
df['x_CO2'] = df['x_CO2'] * 0.01
df['x_CH4'] = df['x_CH4'] * 0.01
df['x_O2'] = df['x_O2'] * 0.01

#Check CH4 and O2 signals are zero
print("CH4, mean x, max x =", df['x_CH4'].mean(),  df['x_CH4'].max())
print("O2, mean x, max x =", df['x_O2'].mean(),  df['x_O2'].max())

# Calculate molar inflow rates [mol/s]
df['F_CH4_0'] = df['v_CH4']*x_CH4_Bottle/V_0/60
df['F_CH4_Ar_0'] = df['v_CH4']*(1 - x_CH4_Bottle)/V_0/60
df['F_CO2_0'] = df['v_CO2']*x_CO2_Bottle/V_0/60
df['F_CO2_Ar_0'] = df['v_CO2']*(1 - x_CO2_Bottle)/V_0/60
df['F_Ar_0'] = df['v_Ar']/V_0/60
df['F_total_0'] = (df['v_CH4'] + df['v_CO2'] + df['v_Ar'])/V_0/60
df['F_Ar_total_0'] = df['F_CH4_Ar_0']+df['F_CO2_Ar_0']+df['F_Ar_0']

# calculate molar outflow rates.
# We see 0 CH4 and 0 H2 after the first cycle. We do have H2 cross sensitivity which gives a small signal.
# H2O has been condensed out which reduces the flow rate which needs to be accounted for.
df['F_H2_out'] = df['x_H2']*(df['F_Ar_total_0']+df['F_CO2_0']+df['F_CH4_0'])/(1-df['x_H2'])
df['F_H2O_out'] = 2*df['F_CH4_0'] - df['F_H2_out']
# total outflow minus H2O
df['F_total_out'] = df['F_CH4_0'] + df['F_CO2_0']  + df['F_Ar_total_0'] + df['F_H2_out']
df['F_CO2_out'] = df['x_CO2']*df['F_total_out']
df['F_CO_out'] = df['x_CO']*df['F_total_out']
df['F_CH4_out'] = df['x_CH4']*df['F_total_out']

# Next analyse the chemical looping cycles by breaking up the data at the switching points between the steps
# where the Ar flow rate set point was 0.1 slm
red1 = df[df['time'] < 1510 ]
ox1 = df[df['time'] > 1510 ]
ox1 = ox1[ox1['time'] < 3720 ]

red2 = df[df['time'] > 3720 ]
red2 = red2[red2['time'] < 4700 ]
ox2 = df[df['time'] > 4700 ]
ox2 = ox2[ox2['time'] < 6190 ]

red3 = df[df['time'] > 6190 ]
red3 = red3[red3['time'] < 7120 ]
ox3 = df[df['time'] > 7120 ]
ox3 = ox3[ox3['time'] < 8750 ]

red4 = df[df['time'] > 8750]
red4 = red4[red4['time'] < 9650 ]
ox4 = df[df['time'] > 9650 ]
ox4 = ox4[ox4['time'] < 11340 ]

red5 = df[df['time'] > 11340]
red5 = red5[red5['time'] < 12250 ]
ox5 = df[df['time'] > 12250 ]


# First cycle was test of CO breakthrough
frames_reduction = [red2, red3, red4, red5]
frames_oxidation = [ox2, ox3, ox4, ox5]


# Make lists to store the results
n_CO2_in = []       # [mol]
n_CO2_out = []       # [mol]
n_CO_out = []       # [mol]
C_balance = []      # [-]
X_CO2 = []          # [-]
X_CO2_peak = []     # [-]

# integrate over time each oxidation step
for frame in frames_oxidation:
    nCO2_0 = integrate.trapz(frame['F_CO2_0'], frame['time'])       # [mol]
    nCO2_f = integrate.trapz(frame['F_CO2_out'], frame['time'])     # [mol]
    nCO_f = integrate.trapz(frame['F_CO_out'], frame['time'])       # [mol]
    n_CO2_in.append(nCO2_0)                                         # [mol]
    n_CO2_out.append(nCO2_f)
    n_CO_out.append(nCO_f)                                          # [mol]
    X_CO2.append(round(nCO_f/(nCO2_f+nCO_f), 5))
    X_CO2_peak.append(frame['x_CO'].max()/0.0502)
    C_balance.append(round((nCO_f+nCO2_f)/(nCO2_0), 5))

print('Results:')
print('X_CO2 = ', X_CO2)
print('X_CO2_peak = ', X_CO2_peak)
print('C_balance = ', C_balance)

n_CH4_in = []        # [mol]
n_CH4_out = []       # [mol]
n_H2O_out = []      # [mol]
n_H2_out = []
red_n_CO2_out = []      # [mol]
red_n_CO_out = []       # [mol]
red_C_balance = []
X_CH4 = []           # [-]
Y_CO2 = []           # [-]
Y_CO = []            # [-]
Y_H2O = []           # [-]

# integrate over time for each reduction step
for frame in frames_reduction:
    nCH4_0 = integrate.trapz(frame['F_CH4_0'], frame['time'])
    nCH4_f = integrate.trapz(frame['F_CH4_out'], frame['time'])
    nCO2_f = integrate.trapz(frame['F_CO2_out'], frame['time'])
    nCO_f = integrate.trapz(frame['F_CO_out'], frame['time'])
    nH2_f = integrate.trapz(frame['F_H2_out'], frame['time'])
    nH2O_f = integrate.trapz(frame['F_H2O_out'], frame['time'])
    n_CH4_in.append(nCH4_0)
    n_CH4_out.append(nCH4_f)
    n_H2O_out.append(nH2O_f)
    red_n_CO2_out.append(nCO2_f) # is equal to extent of reaction
    red_n_CO_out.append(nCO_f)  # is equal to extent of reaction
    n_H2_out.append(nH2_f)
    red_C_balance.append((nCO2_f+nCO_f+nCH4_f)/nCH4_0)
    X_CH4.append(round(1-nCH4_f/nCH4_0, 5))
    Y_CO2.append(round(nCO2_f/nCH4_0, 5))
    Y_CO.append(round(nCO_f / nCH4_0, 5))
    Y_H2O.append(round(nH2O_f / (2 * nCH4_0), 5))


R_feed = np.asarray(n_CO2_in)/np.asarray(n_CH4_in)
print('R_feed (CO2/CH4) = ' , R_feed)
# oxygen mass balance
n_H2O_out = np.asarray(n_H2O_out)
red_n_CO2_out = np.asarray(red_n_CO2_out)
red_n_CO_out = np.asarray(red_n_CO_out)
n_CO_out = np.asarray(n_CO_out)

O_balance = (n_H2O_out + 2*red_n_CO2_out + red_n_CO_out) /n_CO_out

print('X_CH4 = ', X_CH4)
print('Y_CO2 = ', Y_CO2)
print('Y_CO = ', Y_CO)
print('Y_H2O = ', Y_H2O)
print('red_C_balance', red_C_balance)
print('O_balance = ', O_balance)




width = 0.25
#x_label = ['Cycle 1', 'Cycle 2', 'Cycle 3', 'Cycle 4', 'Cycle 5', 'Cycle 6', 'Cycle 7', 'Cycle 8', 'Cycle 9', 'Cycle 10']
x_label = ['CT1', 'CT2', 'CT3', 'Cocurrent']



fig = plt.figure(figsize=(4.3, 3.8))
#plt.plot([0.0+width/2, 1.0+width/2, 2.0 + width/2],[R_feed[1], R_feed[2], 0.666], marker = 's', lw = 0.0, color = 'black', label = 'CO$_2$/H$_2$')
plt.errorbar([1.0, 2.0, 3.0, 4.0], C_balance,  yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'C0', label = 'Oxidation $\mathrm{C}$-balance  ')
plt.errorbar([0.95, 1.95, 2.95, 3.95], red_C_balance,  yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'C1', label = 'Reduction $\mathrm{C}$-balance ')
plt.errorbar([ 1.05, 2.05, 3.05, 4.05],O_balance,yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'black', label = 'Cycle $\mathrm{O}$-balance ')

plt.plot([0.5, 4.5],[1.0,1.0], linewidth=1.0, ls = '--', color = 'grey', label = '__nolegend__')

#plt.xlabel("Cycle")
plt.ylabel("Mass balance [-]")
#plt.title("Change of concentration")

# plt.grid(linestyle='--')
plt.xticks([1.0, 2.0, 3.0, 4.0] , x_label)
plt.legend(loc='lower right', ncol=1)
plt.ylim([0.75, 1.25])
plt.tight_layout()
plt.savefig(os.path.join('plots', 'Demo_930_MB' + '.png'), dpi =400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'Demo_930_MB' + '.pdf'), bbox_inches='tight')
plt.show()


fig = plt.figure(figsize=(4.5, 3.9))

plt.bar(np.arange(len(X_CH4)), X_CH4, color='grey', width=width, edgecolor='black', label='$X_\\mathrm{CH4}$')
plt.bar(np.arange(len(Y_CO2)) + width, Y_CO2, color='C0', width=width, edgecolor='black', label='$Y_\\mathrm{CO_2}$')
plt.bar(np.arange(len(X_CO2)) + 2*width, X_CO2, color='C3', width=width, edgecolor='black', label='Oxidation $X_\\mathrm{CO_2}$')
#plt.plot([0.0+width/2, 1.0+width/2, 2.0 + width/2],[R_feed[1], R_feed[2], 0.666], marker = 's', lw = 0.0, color = 'black', label = 'CO$_2$/H$_2$')
#plt.plot([0.0, 1.0],[X_CO2_peak[1], X_CO2_peak[2]], marker = 's', markersize = 4,  lw = 0.0, color = 'dimgrey', label = '$X_\\mathrm{CO_2}$ peak')
#plt.plot([2.0, 2.0+width],[RWGS_X_CO2, RWGS_X_H2], marker = '+', markersize = 9,  lw = 0.0, color = 'dimgrey', label = 'cofeed TL')

#plt.xlabel("Cycle")
plt.ylabel("Conversion extent [-]")
#plt.title("Change of concentration")
plt.ylim([0,1.3])
# plt.grid(linestyle='--')
plt.xticks(np.arange(len(X_CO2)) + width, x_label)
plt.legend(loc='upper right',
               ncol=2)
plt.tight_layout()

plt.savefig(os.path.join('plots', 'Demo_930_conversion' + '.png'))
plt.savefig(os.path.join('plots', 'Demo_930_conversion' + '.pdf'), bbox_inches='tight')
plt.show()