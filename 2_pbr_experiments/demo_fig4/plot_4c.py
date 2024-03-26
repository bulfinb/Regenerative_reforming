#!/usr/bin/env python
import pandas as pd
import os
import matplotlib.pyplot as plt


filename = "24052023_Demo_930C_dataset.csv"
collumn_names= ['time', 'T', 'p', 'v_CH4', 'v_CO2', 'v_Ar', 'x_H2', 'x_CO2', 'x_CO','x_O2', 'x_CH4', 'Ar MFC sp', 'CH4 MFC sp', 'CO2 MFC sp']

#  load data as pandas arrays and make an instance of the Tga_data object for each file
df = pd.read_csv( filename, skiprows=1, sep=',', names=collumn_names)

# Plot the data
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

#plot temperature
axrr.plot(df['time']-3730, df['T'], lw = 1.0, color = 'C1', label = 'packed bed')

# plot mole fractions
ax.plot(df['time']-3730, df['x_H2'], lw = 1.0, color = 'C4', label = '$x_\mathrm{H_2}$')
ax.plot(df['time']-3730, df['x_CO'], lw = 1.0,color = 'C3', label = '$x_\mathrm{CO}$')
ax.plot(df['time']-3730, df['x_CO2'], lw = 1.0, color = 'C0', label = '$x_\mathrm{CO_2}$')
ax.plot(df['time']-3730, df['x_CH4'], lw = 1.0, color = 'C2', label = '$x_\mathrm{CH_4}$')

# plot flow rates
axr.plot(df['time']-3730, df['v_CH4'], lw = 1.0, color='black', ls ='-', label = "$v_\mathrm{CH_4}$")
axr.plot(df['time']-3730, df['v_CO2'], lw = 1.0, color='dimgrey', ls ='-', label = "$v_\mathrm{CO_2}$")
axr.plot(df['time']-3730, df['v_Ar'], lw = 1.1, color='grey', ls ='--',  label = "$v_\mathrm{Ar}$")

#axr.plot(data['time'], data['v_out'], lw = 1.0, color='grey', ls ='-.',  label = "$v_\mathrm{total}$")

ax.plot([0,18000], [0,0], lw = 1.0,color = 'black', label = '__nolegend__')

# label zones
ax.vlines(x=[960, 2430, 3400,5020, 5940, 7610, 8550], ymin=[0,0,0,0,0,0,0], ymax =[6,6,6,6,6,6,6], colors='grey', ls=':', lw=1, label='__nolegend__')
# label sections
ax.text(300, 4.95, r"red", color ='grey')
ax.text(600, 4.93, r"$\uparrow$", color ='C0', fontsize= 12)
ax.text(1100, 4.95, r"oxidation", color ='grey')
ax.text(2050, 4.93, r"$\downarrow$", color ='red', fontsize= 12)
ax.text(2700, 4.95, r"red", color ='grey')
ax.text(3000, 4.93, r"$\uparrow$", color ='C0', fontsize= 12)
ax.text(3600, 4.95, r"oxidation", color ='grey')
ax.text(4550, 4.93, r"$\downarrow$", color ='red', fontsize= 12)
ax.text(5300, 4.95, r"red", color ='grey')
ax.text(5600, 4.93, r"$\uparrow$", color ='C0', fontsize= 12)
ax.text(6200, 4.95, r"oxidation", color ='grey')
ax.text(7150, 4.93, r"$\downarrow$", color ='red', fontsize= 12)
ax.text(7900, 4.95, r"red", color ='grey')
ax.text(8200, 4.93, r"$\uparrow$", color ='C0', fontsize= 12)
ax.text(8800, 4.95, r"cocurrent oxidation", color ='grey')
ax.text(10800, 4.93, r"$\uparrow$", color ='C0', fontsize= 12)


ax.set_xlim(0,13100)
ax.set_ylim(-2.5,5.3)
axr.set_ylim(0,4.8)


ax.set_yticks([0, 1, 2, 3, 4, 5])
axr.set_yticks([0, 0.4, 0.8, 1.2])
axrr.spines.right.set_bounds((910.0, 945))
axrr.set_ylim(910,1010)
axrr.set_yticks([910, 920, 930, 940])

ax.yaxis.set_label_coords(-0.05,0.65)
axr.yaxis.set_label_coords(-0.07,0.13)
axrr.yaxis.set_label_coords(1.06,0.13)

#ax.set_ylim(0,max(H)+max(H)*0.2)
#axr.set_ylim(0,max(S)+max(S)*0.4)
#ax.set_xlim(min(delta)-0.01,max(delta)+0.01)
ax.legend(loc='upper right',
          ncol=1)
axr.legend(loc='lower right',
          ncol=1)
#axrr.legend(loc='upper center', bbox_to_anchor=(1.1, 1.2), ncol=3, fancybox=True, shadow=True, title='Temperature')



plt.tight_layout()
plt.savefig(os.path.join('plots', '' + 'figure_4c.png'), dpi = 400, bbox_inches='tight')
plt.savefig(os.path.join('plots', 'figure_4c.pdf'), bbox_inches='tight')
plt.show()



