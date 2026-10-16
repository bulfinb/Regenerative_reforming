import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate

print('Counter Current, 3 biogas cycles')

# RELATIVE ERROR IN THE MASS BALANCE
# 2 % error in bottle concentration. 2 % error in gas analysis signals, and 1 % error in the mass flow.
# The errors in the gas-analysis are used twice in each mass balance.
# This neglects the numerical errors introduced during integration of the signals.
relative_error_MB = np.sqrt(2 * (0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))
print("relative error mass balance:", relative_error_MB)
# relative error in n_f outflows is smaller as we only use on gas signal
relative_error_nf = np.sqrt((0.02 ** 2 + 0.02 ** 2 + 0.01 ** 2))

# Bottle concentrations and molar volumes
x_CO2_Bottle = 0.0501  # rest Ar => Accuracy 2%
x_CH4_Bottle = 0.0496  # rest Ar => Accuracy 2%
V_0 = 22.414  # Molar Volume [dm^3/mol] = [l/mol]

# Import Data
filename = "31052023_biogas_930C_dataset.csv"
collumn_names = ['time', 'T', 'p', 'v_CH4', 'v_CO2', 'v_Ar', 'x_H2', 'x_CO2', 'x_CO', 'x_O2', 'x_CH4', 'Ar MFC sp',
                 'CH4 MFC sp', 'CO2 MFC sp']
# Units
# time  s
# T     C
# p     bar (above atmospheric, i.e. p = p_amb - p_reactor)
# v_i   standard litres/min
# x_i   %
# MFC sp  - MFC set point standard litres/min

#  load data as pandas arrays and make an instance of the Tga_data object for each file
df = pd.read_csv(filename, skiprows=1, sep=',', names=collumn_names)

# convert measured % to mole fractions
df['x_H2'] = df['x_H2'] * 0.01
df['x_CO'] = df['x_CO'] * 0.01
df['x_CO2'] = df['x_CO2'] * 0.01
df['x_CH4'] = df['x_CH4'] * 0.01
df['x_O2'] = df['x_O2'] * 0.01

# Check CH4 and O2 are zero
print("CH4, mean x, max x =", df['x_CH4'].mean(), df['x_CH4'].max())
print("O2, mean x, max x =", df['x_O2'].mean(), df['x_O2'].max())

# Calculate molar inflow rates [mol/s]
df['F_CH4_0'] = df['v_CH4'] * x_CH4_Bottle / V_0 / 60
df['F_CH4_Ar_0'] = df['v_CH4'] * (1 - x_CH4_Bottle) / V_0 / 60
df['F_CO2_0'] = df['v_CO2'] * x_CO2_Bottle / V_0 / 60
df['F_CO2_Ar_0'] = df['v_CO2'] * (1 - x_CO2_Bottle) / V_0 / 60
df['F_Ar_0'] = df['v_Ar'] / V_0 / 60
df['F_total_0'] = (df['v_CH4'] + df['v_CO2'] + df['v_Ar']) / V_0 / 60
df['F_Ar_total_0'] = df['F_CH4_Ar_0'] + df['F_CO2_Ar_0'] + df['F_Ar_0']

# calculate molar outflow rates.
# We see 0 CH4 and 0 H2 after the first cycle.
# H2O has been condensed out which reduces the flow rate which needs to be accounted for.

df['F_H2_out'] = df['x_H2'] * (df['F_Ar_total_0'] + df['F_CO2_0'] + df['F_CH4_0']) / (1 - df['x_H2'])
df['F_total_out'] = df['F_CH4_0'] + df['F_CO2_0'] + df['F_Ar_total_0'] + df['F_H2_out']
df['F_CH4_out'] = df['x_CH4'] * df['F_total_out']
df['F_H2O_out'] = 2 * (df['F_CH4_0'] - df['F_CH4_out']) - df['F_H2_out']
# total outflow minus H2O

df['F_CO2_out'] = df['x_CO2'] * df['F_total_out']
df['F_CO_out'] = df['x_CO'] * df['F_total_out']

# Analyse the chemical looping cycles by breaking up the data at the switching points between the steps
# where the Ar flow rate set point was 0.1 slm
timeswitches = []
timeswitches.append(0)
time_skip = 10 * 60  # skip ten minutes to save time
safety = 60  # Safety after switching

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

# Plot the data
# Plot the data
fig = plt.figure(figsize=(10, 4.3))
ax = plt.subplot()
ax.set_xlabel('Time [s]')
ax.set_ylabel('$x_\mathrm{i}$ [%]')

axr = ax.twinx()
axr.yaxis.set_ticks_position('left')
# axr.yaxis.set_ticks_position('left')
axrr = ax.twinx()
axrr.spines['right'].set_color('C1')
axrr.yaxis.label.set_color('C1')
axrr.tick_params(axis='y', colors='C1')
# axrr.spines.right.set_position(("axes", 1.1))

axr.set_ylabel('$v_\mathrm{i}$ [l/min]')
axrr.set_ylabel('$T$ [°C]', color='C1')

axrr.plot(df['time'], df['T'], lw=1.0, color='C1', label='packed bed')

# plot mole fractions
ax.plot(df['time'], 100 * df['x_H2'], lw=1.0, color='C4', label='$x_\mathrm{H_2}$')
ax.plot(df['time'], 100 * df['x_CO'], lw=1.0, color='C3', label='$x_\mathrm{CO}$')
ax.plot(df['time'], 100 * df['x_CO2'], lw=1.0, color='C0', label='$x_\mathrm{CO_2}$')
ax.plot(df['time'], 100 * df['x_CH4'], lw=1.0, color='C2', label='$x_\mathrm{CH_4}$')

axr.plot(df['time'], df['v_CH4'], lw=1.0, color='black', ls='-', label="$v_\mathrm{CH_4}$")
axr.plot(df['time'], df['v_CO2'], lw=1.0, color='grey', ls='-', label="$v_\mathrm{CO_2}$")
axr.plot(df['time'], df['v_Ar'], lw=1.1, color='grey', ls='--', label="$v_\mathrm{Ar}$")

# axr.plot(data['time'], data['v_out'], lw = 1.0, color='grey', ls ='-.',  label = "$v_\mathrm{total}$")

ax.plot([0, 7700], [0, 0], lw=1.0, color='black', label='__nolegend__')
for time in timeswitches:
    ax.axvline(time, -3, 6, ls='--', color='grey')

ax.set_xlim(0, 7700)
ax.set_ylim(-2.5, 5.3)
axr.set_ylim(0, 4.9)

ax.set_yticks([0, 1, 2, 3, 4, 5])
axr.set_yticks([0, 0.4, 0.8, 1.2])
axrr.spines.right.set_bounds((910.0, 945))
axrr.set_ylim(910, 1010)
axrr.set_yticks([910, 920, 930, 940])

ax.yaxis.set_label_coords(-0.05, 0.65)
axr.yaxis.set_label_coords(-0.07, 0.13)
axrr.yaxis.set_label_coords(1.06, 0.13)

# ax.set_ylim(0,max(H)+max(H)*0.2)
# axr.set_ylim(0,max(S)+max(S)*0.4)
# ax.set_xlim(min(delta)-0.01,max(delta)+0.01)
ax.legend(loc='upper right',
          ncol=1)
axr.legend(loc='lower right',
           ncol=1)
# axrr.legend(loc='upper center', bbox_to_anchor=(1.1, 1.2), ncol=3, fancybox=True, shadow=True, title='Temperature')


plt.tight_layout()
plt.savefig(os.path.join('plots', '' + '31052023_Biogas_930C_dataset.png'), dpi=400, bbox_inches='tight')
plt.savefig(os.path.join('plots', '31052023_Biogas_930C_dataset' + '.pdf'), bbox_inches='tight')
plt.show()

# first cycle is reduction (H2) => first timestamp is end of cycle
# store the data for each cycle steps as pandas dataframes in lists
frames = []
frames_reduction = []
frames_oxidation = []

j = 0
while j < (len(timeswitches) - 1):
    # loop to slice up data into sections of oxidation and reduction store in frames
    start_time = timeswitches[j]
    end_time = timeswitches[j + 1]
    dfcut = df[df['time'] > start_time]
    dfcut = dfcut[dfcut['time'] < (end_time + safety)]
    frames.append(dfcut)
    j += 1

# use the set flow rates to separate oxidation from reduction
for frame in frames:
    if frame['CH4 MFC sp'].mean() > 0.1:
        frames_reduction.append(frame)
    if frame['CO2 MFC sp'].mean() > 0.4:
        frames_oxidation.append(frame)

print(len(frames_reduction))
print(len(frames_oxidation))

# Make lists to store the results
cycle_number = [1, 2, 3]
n_CO2_in = []  # [mol]
n_CO2_out = []  # [mol]
n_CO_out = []  # [mol]

# integrate over time each oxidation step
for frame in frames_oxidation:
    nCO2_0 = integrate.trapezoid(frame['F_CO2_0'], frame['time'])  # [mol]
    nCO2_f = integrate.trapezoid(frame['F_CO2_out'], frame['time'])  # [mol]
    nCO_f = integrate.trapezoid(frame['F_CO_out'], frame['time'])  # [mol]
    print(nCO2_0, nCO2_f + nCO_f)
    n_CO2_in.append(nCO2_0)  # [mol]
    n_CO2_out.append(nCO2_f)
    n_CO_out.append(nCO_f)  # [mol]

# frames_oxidation.pop(0)
n_CH4_in = []  # [mol]
n_CH4_out = []  # [mol]
n_H2O_out = []  # [mol]
n_H2_out = []  # [mol]
red_n_CO2_out = []  # [mol]
red_n_CO_out = []  # [mol]
red_n_CO2_in = []  # [mol]
n_total_in_red = []  # [mol]

# integrate over time for each reduction step
for frame in frames_reduction:
    nCH4_0 = integrate.trapezoid(frame['F_CH4_0'], frame['time'])
    nCO2_0 = integrate.trapezoid(frame['F_CO2_0'], frame['time'])
    nCH4_f = integrate.trapezoid(frame['F_CH4_out'], frame['time'])
    nCO2_f = integrate.trapezoid(frame['F_CO2_out'], frame['time'])
    nCO_f = integrate.trapezoid(frame['F_CO_out'], frame['time'])
    nH2_f = integrate.trapezoid(frame['F_H2_out'], frame['time'])
    nH2O_f = integrate.trapezoid(frame['F_H2O_out'], frame['time'])
    n_total_in = integrate.trapezoid(frame['F_total_0'], frame['time'])
    n_total_out = integrate.trapezoid(frame['F_total_out'], frame['time'])
    print(nCH4_0, nCO2_0, nCO2_f, nCO_f, n_total_in, n_total_out)
    n_CH4_in.append(nCH4_0)
    n_CH4_out.append(nCH4_f)
    n_H2O_out.append(nH2O_f)
    red_n_CO2_in.append(nCO2_0)
    red_n_CO2_out.append(nCO2_f)  # is equal to extent of reaction
    red_n_CO_out.append(nCO_f)  # is equal to extent of reaction
    n_H2_out.append(nH2_f)
    n_total_in_red.append(n_total_in)

# frames_reduction.pop(0)

# Per-cycle carbon retention check, using the still-per-cycle lists above
# (before they get collapsed to 3-cycle cumulative totals below), to see
# whether net coke retention is growing cycle-to-cycle or holding steady.
print("--- Per-cycle carbon retention check (coke accumulation) ---")
per_cycle_net_C_red = []
per_cycle_net_C_ox = []
per_cycle_net_C_total = []
for i in range(3):
    net_C_red = (n_CH4_in[i] + red_n_CO2_in[i]) - (n_CH4_out[i] + red_n_CO2_out[i] + red_n_CO_out[i])
    net_C_ox = n_CO2_in[i] - (n_CO2_out[i] + n_CO_out[i])
    net_C_total = net_C_red + net_C_ox
    per_cycle_net_C_red.append(net_C_red)
    per_cycle_net_C_ox.append(net_C_ox)
    per_cycle_net_C_total.append(net_C_total)
    print(f'Cycle {cycle_number[i]}: net_C_retained_red = {net_C_red:.6f} mol, '
          f'net_C_retained_ox = {net_C_ox:.6f} mol, net_C_retained_total = {net_C_total:.6f} mol')
print('Cumulative net C retained after 3 cycles [mol] = ', np.sum(per_cycle_net_C_total))

# Mean per-cycle inlet amounts, for use as simulation inputs (matches the
# "mean integrated inlet amounts" convention used in study_bulfin.py).
n_Ar_in_red = [n_total_in_red[i] - n_CH4_in[i] - red_n_CO2_in[i] for i in range(3)]
print("--- Mean per-cycle inlet amounts (simulation inputs) ---")
print('n_CH4_in_red (mean/cycle) [mol] =      ', np.mean(n_CH4_in))
print('n_CO2_in_red (mean/cycle) [mol] =      ', np.mean(red_n_CO2_in))
print('n_Ar_in_red (mean/cycle) [mol] =       ', np.mean(n_Ar_in_red))
print('n_CO2_in_ox (mean/cycle) [mol] =       ', np.mean(n_CO2_in))

n_CH4_in = np.asarray(n_CH4_in).sum()
n_CO2_in = np.asarray(n_CO2_in).sum()
n_CH4_out = np.asarray(n_CH4_out).sum()
n_CO_out = np.asarray(n_CO_out).sum()
n_CO2_out = np.asarray(n_CO2_out).sum()
n_H2O_out = np.asarray(n_H2O_out).sum()
red_n_CO2_in = np.asarray(red_n_CO2_in).sum()
red_n_CO2_out = np.asarray(red_n_CO2_out).sum()
red_n_CO_out = np.asarray(red_n_CO_out).sum()

# bed loading, from reactor_bed
M_CeO2 = 140.116 + 2 * 15.999  # g/mol
m_CeO2 = 292.3  # g
n_CeO2 = m_CeO2 / M_CeO2  # mol

print("--- Input molar quantities (summed over 3 cycles) ---")
print('n_CH4_in_red [mol] =  ', n_CH4_in)
print('n_CO2_in_red [mol] =  ', red_n_CO2_in)
print('n_CO2_in_ox [mol] =   ', n_CO2_in)
print('n_CeO2 [mol] =        ', n_CeO2, ' (m_CeO2 =', m_CeO2, 'g, M_CeO2 =', M_CeO2, 'g/mol)')

# oxygen mass balance
print("3 cycle cumulative mass balance results")
print('R_feed (CO2/CH4) =    ', n_CO2_in / n_CH4_in)
print('R_prod (CO/CH4) =     ', n_CO_out / n_CH4_in)
print('Reduction X_CH4 =     ', 1 - n_CH4_out / n_CH4_in)
print('Reduction Y_CO2 =     ', (red_n_CO2_out - red_n_CO2_in) / n_CH4_in)
print('Reduction Y_H2O =     ', n_H2O_out / (2 * n_CH4_in))
print('Reduction C_balance = ', (red_n_CO_out + red_n_CO2_out + n_CH4_out) / (n_CH4_in + red_n_CO2_in))
print('Reduction C_loss = ', (n_CH4_in + red_n_CO2_in) - (red_n_CO_out + red_n_CO2_out + n_CH4_out))
print('m_C [mg] = ', ((n_CH4_in + red_n_CO2_in) - (red_n_CO_out + red_n_CO2_out + n_CH4_out)) * 1000 * 12)
print('Oxidation X_CO2 =     ', 1 - n_CO2_out / n_CO2_in)
print('Oxidation Y_CO =      ', n_CO_out / n_CO2_in)
print('Oxidation C_balance = ', (n_CO_out + n_CO2_out) / n_CO2_in)
print('Total O_balance =     ', (n_H2O_out + 2 * (red_n_CO2_out - red_n_CO2_in) + red_n_CO_out) / n_CO_out)
