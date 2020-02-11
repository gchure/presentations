# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load data and correct names
data = pd.read_csv('../data/RazoMejia2018_data.csv')
data = data[data['repressors']> 0].copy()
data['repressors'] *= 2
data.rename(columns={'IPTG_uM':'IPTGuM'}, inplace=True)

# Compute the summary statistics
summary = data.groupby(['operator', 'repressors', 
            'IPTGuM']).agg(('mean','sem')).reset_index()

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(1, 1, figsize=(3, 3))

# Set limits
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')
ax.set_ylim(-0.1, 1.2)
ax.set_xlim(-10, 12)

# Define the colors
rep_colors = {22:colors['dark_purple'], 60:colors['red'], 
              124:colors['black'], 260:colors['orange'],
              1220:colors['blue'], 1740:colors['green'] }
face_colors = {22:colors['light_purple'], 60:colors['light_red'], 
              124:colors['light_grey'], 260:colors['light_orange'],
              1220:colors['light_blue'], 1740:colors['light_green'] }
op_glyphs = {'O1':'s', 'O2':'o', 'O3':'D'}
DATA = 1
FIT_ONLY = 1
FIT_STRAIN = [260, 'O2']
# ##############################################################################
# MASTER CURVE
# ##############################################################################
F = np.linspace(-10, 12, 200)
master = (1 + np.exp(-F))**-1
ax.plot(F, master, 'k-', lw=1)


# ##############################################################################
# DATA
# ##############################################################################
if DATA == 1:
    for g, d in summary.groupby(['repressors', 'operator']):
        bohr = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]],
                                          ka=constants['Ka'], ki=constants['Ki'],
                                          ep_ai=constants['ep_AI'],
                            effector_conc=d['IPTGuM']).bohr_parameter()
        if (g[0] == FIT_STRAIN[0]) & (g[1]== FIT_STRAIN[1]): 
            face = 'w'
            plot = 1
        else:
            face = face_colors[g[0]]
            plot = 0

        if FIT_ONLY == 1:
            if plot == 1:
                ax.errorbar(bohr, d['fold_change_A']['mean'], d['fold_change_A']['sem'],
                fmt=op_glyphs[g[1]], lw=0.75, capsize=1, color=rep_colors[g[0]],
                markerfacecolor=face, markeredgewidth=0.5)
        else:
            ax.errorbar(bohr, d['fold_change_A']['mean'], d['fold_change_A']['sem'],
                fmt=op_glyphs[g[1]], lw=0.75, capsize=1, color=rep_colors[g[0]],
                markerfacecolor=face, markeredgewidth=0.5)

# ##############################################################################
# LEGEND AND SAVING
# ##############################################################################
plt.tight_layout()

if (DATA == 0):
    plt.savefig('../figs/master_curve.pdf', bbox_inches='tight')
elif (FIT_ONLY == 1):
    plt.savefig('../figs/collapse_fit_strain.pdf', bbox_inches='tight')
else:
    plt.savefig('../figs/collapse_all_data.pdf', bbox_inches='tight')