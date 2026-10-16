#!/usr/bin/env python
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np


filename = "rx_test_027.lvm"

# Only select the collumns we need. Mario uses a lot of post processing that is not relavent for our data
collumn_names= ['time', 'CH4 MFC', 'v_CH4', 'CO2 MFC', 'v_CO2', 'Ar MFC', 'v_Ar', 'v_out', 'x_CH4', 'x_CO2', 'x_H2', 'x_CO','x_O2', 'p', 'T']
usecollumns = [0,4,5,7,8,10,11,13,18,20,22,24,26,32,34]
rename_dict = {}
for i, col in enumerate(usecollumns):
    rename_dict[col] = collumn_names[i]

print(rename_dict)
# v_i   litres/min
# p     bar
# T     C
# x_i   %
#  load data as pandas arrays and make an instance of the Tga_data object for each file
df = pd.read_csv(filename, sep='\t', decimal =',', usecols = [0,4,5,7,8,10,11,13,18,20,22,24,26,32,34], header= None)
df.rename(rename_dict, axis = 'columns', inplace = True)
print(df)

# subtract background H2
df['x_H2'] = df['x_H2'] - 0.15

print(df)

# subtract cross sensitivity for CO and CO2 from H2 signal with time delay
CO_H2_cross = 0.07
CO2_H2_cross = 0.02
time_delay = 35
CO_cross = np.interp(df['time'], df['time']-time_delay,df['x_CO'])*CO_H2_cross
CO2_cross = np.interp(df['time'], df['time']-time_delay,df['x_CO2'])*CO2_H2_cross
df['x_H2'] = df['x_H2'] - CO_cross - CO2_cross


df.to_csv("31052023_biogas_930C_dataset.csv", index = False,
            columns = ['time', 'T', 'p', 'v_CH4', 'v_CO2', 'v_Ar', 'x_H2', 'x_CO2', 'x_CO','x_O2', 'x_CH4', 'Ar MFC', 'CH4 MFC', 'CO2 MFC' ],
            header=['time [s] ', 'T [C] ', 'p-p_amb [bar] ', 'v_CH4 [slm] ', 'v_CO2 [slm] ', 'v_Ar [slm] ',
                    'x_H2 [%] ', 'x_CO2 [%] ', 'x_CO [%] ','x_O2 [%] ', 'x_CH4 [%] ', 'Ar MFC sp [slm]', 'CH4 MFC sp [slm]', 'CO2 MFC sp [slm]'] )

# Identify the reduction/oxidation step boundaries. Purges between steps show
# up as an Ar flow rate set point of 0.1 slm; the switch point is taken as the
# end of each purge (last sample with Ar MFC == 0.1), i.e. the moment the next
# step's reactant gas actually starts flowing.
ar_purge = df['Ar MFC'] == 0.1
purge_end = ar_purge & ~ar_purge.shift(-1, fill_value=False)
timeswitches = [df['time'].iloc[0]] + df.loc[purge_end, 'time'].tolist() + [df['time'].iloc[-1]]

# Classify each interval as reduction (CH4 feed) or oxidation (CO2 feed), by
# the mean MFC set point over that interval.
step_spans = []
for start_time, end_time in zip(timeswitches[:-1], timeswitches[1:]):
    mask = (df['time'] >= start_time) & (df['time'] < end_time)
    if df.loc[mask, 'CH4 MFC'].mean() > 0.1:
        step_spans.append((start_time, end_time, 'reduction'))
    elif df.loc[mask, 'CO2 MFC'].mean() > 0.4:
        step_spans.append((start_time, end_time, 'oxidation'))

# Plot the data
fig = plt.figure(figsize=(10,3.87))
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

# Shade the reduction/oxidation steps and label their flow direction
# (countercurrent: reduction and oxidation gas feeds enter from opposite
# ends of the bed).
step_colors = {'reduction': 'lightblue', 'oxidation': 'lightgrey'}
step_labels = {'reduction': r'reduction $\uparrow$', 'oxidation': r'oxidation $\downarrow$'}
label_transform = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
for start_time, end_time, step_type in step_spans:
    ax.axvspan(start_time, end_time, color=step_colors[step_type], alpha=0.4, zorder=0)
    ax.text((start_time + end_time) / 2, 0.96, step_labels[step_type],
            transform=label_transform, ha='center', va='top', fontsize=10)

# plot mole fractions
ax.plot(df['time'], df['x_H2'], lw = 1.0, color = 'C4', label = '$x_\mathrm{H_2}$')
ax.plot(df['time'], df['x_CO'], lw = 1.0,color = 'C3', label = '$x_\mathrm{CO}$')
ax.plot(df['time'], df['x_CO2'], lw = 1.0, color = 'C0', label = '$x_\mathrm{CO_2}$')
ax.plot(df['time'], df['x_CH4'], lw = 1.0, color = 'C2', label = '$x_\mathrm{CH_4}$')





axr.plot(df['time'], df['v_CH4'], lw = 1.0, color='black', ls ='-', label = "$v_\mathrm{CH_4}$")
axr.plot(df['time'], df['v_CO2'], lw = 1.0, color='dimgrey', ls ='-', label = "$v_\mathrm{CO_2}$")
axr.plot(df['time'], df['v_Ar'], lw = 1.1, color='grey', ls ='--',  label = "$v_\mathrm{Ar}$")

#axr.plot(data['time'], data['v_out'], lw = 1.0, color='grey', ls ='-.',  label = "$v_\mathrm{total}$")

ax.plot([0,7700], [0,0], lw = 1.0,color = 'black', label = '__nolegend__')

ax.set_xlim(0,7700)
ax.set_ylim(-2.5,6)
axr.set_ylim(0,3.6)


ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
axr.set_yticks([0, 0.3, 0.6, 0.9])
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
plt.savefig(os.path.join('plots', '' + '31052023_Biogas_930C_dataset.png'), dpi = 400, bbox_inches='tight')
plt.savefig(os.path.join('plots', '31052023_Biogas_930C_dataset' + '.pdf'), bbox_inches='tight')
plt.show()



