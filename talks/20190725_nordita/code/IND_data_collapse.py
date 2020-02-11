# -*- coding: utf-8 -*-
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.legend
import phd.viz
import phd.thermo
import seaborn as sns
constants = phd.thermo.load_constants()
colors = phd.viz.phd_style()
_colors = sns.color_palette('magma', n_colors=3)

# Load the data
data = pd.read_csv('../data/summarized_data.csv')
data = data[data['class']=='IND'].copy()
stats = pd.read_csv('../data/KaKi_epAI_summary.csv')
stats = stats[stats['operator']=='O2'].copy()

# Define colors and glyphs
op_colors = {'O1':_colors[0], 'O2':_colors[1], 'O3':_colors[2]}
glyphs = {'F164T': 's', 'Q294V':'o', 'Q294K':'^', 'Q294R':'x'}

# Other plotting constants. 
bohr_range = np.linspace(-10, 10, 500)

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(1,1,figsize=(3, 3))
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')

# ##############################################################################
# MASTER CURVE
# ##############################################################################
ax.plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', lw=1, 
        label='__nolegend__')

# ##############################################################################
#  COLLAPSE DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'operator']):
    # Extract the correct statistics.
    _stats = stats[(stats['mutant']==g[0])]
    ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
    epAI = _stats[_stats['parameter']=='ep_AI']['median'].values[0]

    # Compute the bohr parameter
    bohr = phd.thermo.SimpleRepression(R=260, ep_r=constants[g[1]], ka=ka,
                                       ki=ki, 
                                       ep_ai=epAI,
                                       effector_conc=d['IPTGuM']).bohr_parameter()

    # Plot the bohr parameter and errors
    ax.errorbar(bohr, d['mean'], d['sem'], lw=0.75, color=op_colors[g[1]],
                fmt=glyphs[g[0]], linestyle='none', 
                capsize=1, label='__nolegend__', ms=4)

# ##############################################################################
# COMPLICATED LEGEND STUFF
# ##############################################################################
for r, c in op_colors.items():
    ax.plot([], [], lw=2, color=c, label=r)
leg = ax.legend(title='operator')
leg.get_title().set_fontsize(6)

# Add the second legend
glyph_leg = []
for m, g in glyphs.items():
    glyph_leg += ax.plot([], [], marker=g, markerfacecolor='w', color='k', markeredgewidth=0.75)

leg = matplotlib.legend.Legend(ax, glyph_leg, ['F164T', 'Q294V', 'Q294K', 'Q294R'], 
                               loc='center left')
ax.add_artist(leg)
plt.tight_layout()
plt.savefig('../figs/IND_collapse.pdf', bbox_inches='tight')
