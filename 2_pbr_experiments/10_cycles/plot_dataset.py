#!/usr/bin/env python
import pandas as pd
import os
import matplotlib.pyplot as plt

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


fig = plt.figure(figsize=(10,4.3))
ax = plt.subplot()
ax.set_xlabel('Time [s]')
ax.set_ylabel('$x_\mathrm{i}$ [%]')

axr = ax.twinx()
axr.yaxis.set_ticks_position('left')
#axr.yaxis.set_ticks_position('left')
axrr = ax.twinx()
axrr.spines['right'].set_color('C1')
axrr.yaxis.label.set_color('C1')
axrr.tick_params(axis='y', colors='C1')
#axrr.spines.right.set_position(("axes", 1.1))
axr.set_ylabel('$v_\mathrm{i}$ [l/min]')
axrr.set_ylabel('$T$ [°C]', color ='C1')

axrr.plot(df['time'], df['T'], lw = 1.0, color = 'C1', label = 'packed bed')

# plot mole fractions
ax.plot(df['time'], df['x_H2'], lw = 1.0, color = 'C4', label = '$x_\mathrm{H_2}$')
ax.plot(df['time'], df['x_CO'], lw = 1.0,color = 'C3', label = '$x_\mathrm{CO}$')
ax.plot(df['time'], df['x_CO2'], lw = 1.0, color = 'C0', label = '$x_\mathrm{CO_2}$')
ax.plot(df['time'], df['x_CH4'], lw = 1.0, color = 'C2', label = '$x_\mathrm{CH_4}$')

axr.plot(df['time'], df['v_CH4'], lw = 1.0, color='black', ls ='-', label = "$v_\mathrm{CH_4}$")
axr.plot(df['time'], df['v_CO2'], lw = 1.0, color='dimgrey', ls ='-', label = "$v_\mathrm{CO_2}$")
axr.plot(df['time'], df['v_Ar'], lw = 1.1, color='grey', ls ='--',  label = "$v_\mathrm{Ar}$")

#axr.plot(data['time'], data['v_out'], lw = 1.0, color='grey', ls ='-.',  label = "$v_\mathrm{total}$")

ax.plot([0,16000], [0,0], lw = 1.0,color = 'black', label = '__nolegend__')

ax.set_xlim(0,26000)
ax.set_ylim(-2.5,5.3)
axr.set_ylim(0,3.6)
axrr.spines.right.set_bounds((915, 945))
axrr.set_ylim(915,1000)

ax.set_yticks([0, 1, 2, 3, 4, 5])
axr.set_yticks([0, 0.3, 0.6, 0.9])
axrr.set_yticks([915, 925, 935])

ax.yaxis.set_label_coords(-0.05,0.65)
axr.yaxis.set_label_coords(-0.07,0.13)
axrr.yaxis.set_label_coords(1.06,0.13)



ax.hlines(y=0.0, xmin=0, xmax=29000, linewidth=1.0, color='grey')
ax.legend(loc='upper right',
          ncol=1)
axr.legend(bbox_to_anchor=(0.88, 0.25),
          ncol=1)


plt.tight_layout()
plt.savefig(os.path.join('plots', '' + '30052023_10cycle_dataset.png'), dpi = 400, bbox_inches='tight')
plt.savefig(os.path.join('plots', '30052023_10cycle_dataset.pdf'), bbox_inches='tight')
plt.show()



