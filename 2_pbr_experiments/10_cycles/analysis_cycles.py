import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cantera as ct
from scipy import integrate
import os

print('Counter Current, 10 cycles')

# RELATIVE ERROR IN THE MASS BALANCE
# 2 % error in bottle concentration. 2 % error in gas analysis signals, and 1 % error in the mass flow.
# The errors in the gas-analysis are used twice in each mass balance.
# This neglects the numerical errors introduced during integration of the signals.
relative_error_MB = np.sqrt(2 * (0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))
print("relative error mass balance:", relative_error_MB)
# relative error in n_f outflows is smaller as we only use on gas signal
relative_error_nf = np.sqrt((0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))

# Bottle concentrations and molar volumes
x_CO2_Bottle = 0.0501           # rest Ar => Accuracy 2%
x_CH4_Bottle = 0.0496            # rest Ar => Accuracy 2%
V_0 = 22.414                    # Molar Volume [dm^3/mol] = [l/mol]



# Import Data
filename = "30052023_10cycle_930C_dataset.csv"
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
df['x_H2'] = df['x_H2']*0.01
df['x_CO'] = df['x_CO']*0.01
df['x_CO2'] = df['x_CO2']*0.01
df['x_CH4'] = df['x_CH4']*0.01
df['x_O2'] = df['x_O2']*0.01

#Check CH4 and O2 are zero
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


# Analyse the chemical looping cycles by breaking up the data at the switching points between the steps
# where the Ar flow rate set point was 0.1 slm
timeswitches = []
timeswitches.append(0)
time_skip = 10*60  # skip ten minutes to save time
safety = 60        # Safety after switching

i = 0
while i < len(df['Ar MFC sp']):
    Ar_flow = df['Ar MFC sp'].iloc[i]
    if Ar_flow == 0.1:
        timeswitches.append(df['time'].iloc[i])
        i += time_skip  # skip 10min
    else:
        i += 5  # skip 5sec to next value

timeswitches.append(df['time'].iloc[-1])
print("timeswitches ", timeswitches)
# first cycle is reduction (H2) => first timestamp is end of cycle

# store the data for each cycle steps as pandas dataframes in lists
frames = []
frames_reduction = []
frames_oxidation = []


j = 0
while j < (len(timeswitches)-1):
    # loop to slice up data into sections of oxidation and reduction store in frames
    start_time = timeswitches[j]
    end_time = timeswitches[j+1]
    dfcut = df[df['time'] > start_time]
    dfcut = dfcut[dfcut['time'] < (end_time + safety)]
    frames.append(dfcut)
    j += 1

# use the set flow rates to separate oxidation from reduction
for frame in frames:
    if frame['CH4 MFC sp'].mean() > 0.1:
        frames_reduction.append(frame)
    if frame['CO2 MFC sp'].mean() > 0.1:
        frames_oxidation.append(frame)






# Make lists to store the results
cycle_number = [1,2,3,4,5,6,7,8,9,10]
n_CO2_in = []       # [mol]
n_CO2_out = []       # [mol]
n_CO_out = []       # [mol]
X_CO2 = []          # [-]
C_balance = []


# integrate over time each oxidation step
for frame in frames_oxidation:
    nCO2_0 = integrate.trapz(frame['F_CO2_0'], frame['time'])       # [mol]
    nCO2_f = integrate.trapz(frame['F_CO2_out'], frame['time'])     # [mol]
    nCO_f = integrate.trapz(frame['F_CO_out'], frame['time'])       # [mol]
    n_CO2_in.append(nCO2_0)                                         # [mol]
    n_CO2_out.append(nCO2_f)
    n_CO_out.append(nCO_f)                                          # [mol]
    X_CO2.append(round(nCO_f / (nCO2_f + nCO_f), 5))
    C_balance.append(round((nCO_f+nCO2_f)/(nCO2_0), 5))


n_CH4_in = []        # [mol]
n_CH4_out = []       # [mol]
n_H2O_out = []      # [mol]
n_H2_out = []
red_n_CO2_out = []      # [mol]
red_n_CO_out = []       # [mol]
X_CH4 = []           # [-]
Y_CO2 = []           # [-]
red_C_balance = []

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
    X_CH4.append(round(1 - nCH4_f / nCH4_0, 5))
    Y_CO2.append(round(nCO2_f / nCH4_0, 5))
    red_C_balance.append((nCO2_f + nCO_f + nCH4_f) / nCH4_0)

n_CH4_in = np.asarray(n_CH4_in)
n_CO2_in =  np.asarray(n_CO2_in)
n_CH4_out = np.asarray(n_CH4_out)
n_CO_out = np.asarray(n_CO_out)
n_CO2_out = np.asarray(n_CO2_out)
n_H2O_out = np.asarray(n_H2O_out)
red_n_CO2_out = np.asarray(red_n_CO2_out)
red_n_CO_out = np.asarray(red_n_CO_out)


O_balance = (n_H2O_out + 2*red_n_CO2_out + red_n_CO_out) /n_CO_out

# Make plots
width = 0.25
x_label = ['1','2', ' 3', ' 4', ' 5', ' 6', ' 7', ' 8', ' 9', '10']

fig = plt.figure(figsize=(4.3, 3.8))
#plt.plot([0.0+width/2, 1.0+width/2, 2.0 + width/2],[R_feed[1], R_feed[2], 0.666], marker = 's', lw = 0.0, color = 'black', label = 'CO$_2$/H$_2$')
plt.errorbar([1.0, 2.0, 3.0, 4.0 , 5.0, 6.0, 7.0 , 8.0, 9.0,10.0], C_balance,  yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'C0', label = 'Oxidation $\mathrm{C}$-balance ')
plt.errorbar([0.95, 1.95, 2.95, 3.95, 4.95, 5.95, 6.95, 7.95, 8.95, 9.95], red_C_balance,  yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'C1', label = 'Reduction $\mathrm{C}$-balance')
plt.errorbar([ 1.05, 2.05, 3.05, 4.05,5.05, 6.05, 7.05 , 8.05, 9.05, 10.05],O_balance,yerr=relative_error_MB, elinewidth=1.0, capsize=2,
             marker = 's', lw = 0.0, color = 'black', label = 'Cycle $\mathrm{O}$-balance')

plt.plot([0.5, 10.5],[1.0,1.0], linewidth=1.0, ls = '--', color = 'grey', label = '__nolegend__')

plt.xlabel("Cycle #")
plt.ylabel("Mass balance [-]")
#plt.title("Change of concentration")

# plt.grid(linestyle='--')
plt.xticks([1.0, 2.0, 3.0, 4.0 , 5.0, 6.0, 7.0 , 8.0, 9.0, 10.0] , x_label)
plt.legend(loc='lower right', ncol=1)
plt.ylim([0.75, 1.15])
plt.tight_layout()
plt.savefig(os.path.join('plots', '10cycle_930_MB' + '.png'), dpi =400, bbox_inches='tight')
plt.savefig(os.path.join('plots', '10cycle_930_MB' + '.pdf'), bbox_inches='tight')
plt.show()


fig = plt.figure(figsize=(4.5, 3.9))

plt.bar(np.arange(len(X_CH4)), X_CH4, color='grey', width=width, edgecolor='black', label='$X_\\mathrm{CH4}$')
plt.bar(np.arange(len(Y_CO2)) + width, Y_CO2, color='C0', width=width, edgecolor='black', label='$Y_\\mathrm{CO_2}$')
plt.bar(np.arange(len(X_CO2)) + 2*width, X_CO2, color='C3', width=width, edgecolor='black', label='Oxidation $X_\\mathrm{CO_2}$')
#plt.plot([0.0+width/2, 1.0+width/2, 2.0 + width/2],[R_feed[1], R_feed[2], 0.666], marker = 's', lw = 0.0, color = 'black', label = 'CO$_2$/H$_2$')
#plt.plot([0.0, 1.0],[X_CO2_peak[1], X_CO2_peak[2]], marker = 's', markersize = 4,  lw = 0.0, color = 'dimgrey', label = '$X_\\mathrm{CO_2}$ peak')
#plt.plot([2.0, 2.0+width],[RWGS_X_CO2, RWGS_X_H2], marker = '+', markersize = 9,  lw = 0.0, color = 'dimgrey', label = 'cofeed TL')

plt.xlabel("Cycle #")
plt.ylabel("Conversion extent [-]")
#plt.title("Change of concentration")
plt.ylim([0,1.3])
# plt.grid(linestyle='--')
plt.xticks(np.arange(len(X_CO2)) + width, x_label)
plt.legend(loc='upper right',
               ncol=2)
plt.tight_layout()

plt.savefig(os.path.join('plots', '10cycle_930_conversion' + '.png'))
plt.savefig(os.path.join('plots', '10cycle_930_conversion' + '.pdf'), bbox_inches='tight')
plt.show()


# CUMULATIVE MASS BALANCE CYLES 2-10
# remove cycle 1
frames_oxidation.pop(0)
frames_reduction.pop(0)

n_CO2_in = []       # [mol]
n_CO2_out = []       # [mol]
n_CO_out = []       # [mol]


# integrate over time each oxidation step
for frame in frames_oxidation:
    nCO2_0 = integrate.trapz(frame['F_CO2_0'], frame['time'])       # [mol]
    nCO2_f = integrate.trapz(frame['F_CO2_out'], frame['time'])     # [mol]
    nCO_f = integrate.trapz(frame['F_CO_out'], frame['time'])       # [mol]
    n_CO2_in.append(nCO2_0)                                         # [mol]
    n_CO2_out.append(nCO2_f)
    n_CO_out.append(nCO_f)                                          # [mol]


n_CH4_in = []        # [mol]
n_CH4_out = []       # [mol]
n_H2O_out = []      # [mol]
n_H2_out = []
red_n_CO2_out = []      # [mol]
red_n_CO_out = []       # [mol]
X_CH4 = []           # [-]
Y_CO2 = []           # [-]
red_C_balance = []

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

n_CH4_in = np.asarray(n_CH4_in)
n_CO2_in =  np.asarray(n_CO2_in)
n_CH4_out = np.asarray(n_CH4_out)
n_CO_out = np.asarray(n_CO_out)
n_CO2_out = np.asarray(n_CO2_out)
n_H2O_out = np.asarray(n_H2O_out)
red_n_CO2_out = np.asarray(red_n_CO2_out)
red_n_CO_out = np.asarray(red_n_CO_out)

# CUMULATIVE MASS BALANCE FOR TABLE 2
print("2-10 cycle cumulative mass balance results")
n_CH4_in = n_CH4_in.sum()
n_CO2_in =  n_CO2_in.sum()
n_CH4_out = n_CH4_out.sum()
n_CO_out = n_CO_out.sum()
n_CO2_out = n_CO2_out.sum()
n_H2O_out = n_H2O_out.sum()
red_n_CO2_out = red_n_CO2_out.sum()
red_n_CO_out = red_n_CO_out.sum()

print('R_feed (CO2/CH4) =    ' , n_CO2_in / n_CH4_in)
print('R_prod (CO/CH4) =     ', n_CO_out / n_CH4_in)
print('Reduction X_CH4 =     ', 1 - n_CH4_out/n_CH4_in)
print('Reduction Y_CO2 =     ', red_n_CO2_out/n_CH4_in)
print('Reduction Y_H2O =     ', n_H2O_out/(2 * n_CH4_in))
print('Reduction C_balance = ', (red_n_CO_out + red_n_CO2_out + n_CH4_out)/n_CH4_in)
print('Oxidation X_CO2 =     ', 1 - n_CO2_out/n_CO2_in)
print('Oxidation Y_CO =      ', n_CO_out/n_CO2_in)
print('Oxidation C_balance = ', (n_CO_out+n_CO2_out)/n_CO2_in)
print('Total O_balance =     ', (n_H2O_out + 2*red_n_CO2_out + red_n_CO_out) /n_CO_out)



