# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
constants = phd.thermo.load_constants()
constants['Oid'] = -17.3
colors = phd.viz.phd_style()

data = pd.read_csv('../data/Garcia2011_Brewster2014_data.csv')

FIT_STRAIN = 1
ALL_DATA = 0

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('repressors per cell')
ax.set_ylabel('fold-change')

R = np.logspace(0, 4, 500)

op_colors = {'Oid': colors['black'], 'O1':colors['purple'], 'O2':colors['orange'],
            'O3':colors['blue']}
face_colors = {'Oid':colors['light_grey'], 'O1':colors['light_purple'],
               'O2':colors['light_orange'], 'O3':colors['light_blue']}
glyphs = {'garcia':'o', 'brewster': 'd'}

# Theory curves

for o in ['O1', 'O2', 'O3', 'Oid']:
    arch = phd.thermo.SimpleRepression(R=R, ep_r=constants[o],
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'],
                                       effector_conc=0).fold_change()
    ax.plot(R, arch, '-', lw=1, label=o, color=op_colors[o])

ax.plot([], [], 'ko', markerfacecolor='w', label='Garcia & Phillips 2011', ms=2)
ax.plot([], [], 'kD', markerfacecolor='w', label='Brewster et al. 2014', ms=2)
for g, d in data.groupby(['operator', 'author']):
    ax.plot(d['repressor'], d['fold_change'], marker=glyphs[g[1]], markerfacecolor=face_colors[g[0]],
            color=op_colors[g[0]], markeredgewidth=0.75, label='__nolegend__',
            linestyle='none', ms=4)

ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig('../figs/old_gods_titration.pdf', bbox_inches='tight')
