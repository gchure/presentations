# -*- coding: utf-8 -*-
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.legend
import phd.viz
import phd.thermo
constants = phd.thermo.load_constants()
colors = phd.viz.phd_style()

# Load the data
data = pd.read_csv('../data/summarized_data.csv')
data = data[data['class']=='DNA'].copy()
stats = pd.read_csv('../data/DNA_binding_energy_summary.csv')
stats = stats[stats['repressors']==260].copy()

# Define colors and glyphs
rep_colors = {60:colors['red'], 124:colors['black'], 260:colors['orange'],
              1220:colors['blue']} 
face_colors = {60:colors['light_red'], 124:colors['light_grey'], 
              260:colors['light_orange'], 1220:colors['light_blue']} 
glyphs = {'Q21M': 's', 'Y20I':'o', 'Q21A':'^'}

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
for g, d in data.groupby(['mutant', 'repressors']):
    # Extract the correct statistics.
    epRA = stats[(stats['mutant']==g[0]) & (stats['parameter']=='ep_RA')]

    # Compute the bohr parameter
    c, ep = np.meshgrid(d['IPTGuM'], epRA[['median', 'hpd_min', 'hpd_max']].values)
    bohr = phd.thermo.SimpleRepression(R=g[1], ep_r=ep, ka=constants['Ka'],
                                       ki=constants['Ki'], 
                                       ep_ai=constants['ep_AI'],
                                       effector_conc=c).bohr_parameter()

    # Plot the bohr parameter and errors
    ax.errorbar(bohr[0, :], d['mean'], d['sem'], lw=0.75, color=rep_colors[g[1]],
                markerfacecolor=face_colors[g[1]], fmt=glyphs[g[0]], linestyle='none', 
                capsize=1, label='__nolegend__', ms=4)
    ax.hlines(d['mean'], bohr[1, :], bohr[2, :], lw=0.75, color=rep_colors[g[1]],
            label='__nolegend__')

# ##############################################################################
# COMPLICATED LEGEND STUFF
# ##############################################################################
for r, c in rep_colors.items():
    ax.plot([], [], lw=2, color=c, label=int(r))
leg = ax.legend(title='rep. per cell')
leg.get_title().set_fontsize(6)

# Add the second legend
glyph_leg = []
for m, g in glyphs.items():
    glyph_leg += ax.plot([], [], marker=g, markerfacecolor='w', color='k', markeredgewidth=0.75)

leg = matplotlib.legend.Legend(ax, glyph_leg, ['Q21M', 'Y20I', 'Q21A'], 
                               loc='center left')
ax.add_artist(leg)
plt.tight_layout()
plt.savefig('../figs/DNA_collapse.pdf', bbox_inches='tight')
