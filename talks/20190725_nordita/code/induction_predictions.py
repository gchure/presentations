# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
import seaborn as sns
import pickle
import imp
imp.reload(phd.viz)
colors = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load data and correct names
data = pd.read_csv('../data/RazoMejia2018_data.csv')
data = data[data['repressors']> 0].copy()
data['repressors'] *= 2
data.rename(columns={'IPTG_uM':'IPTGuM'}, inplace=True)

# Load the MCMC chains to draw credible regions. 
with open('../data/main_text_KaKi.pkl', 'rb') as pkl:
    chain = pickle.load(pkl)
ka_chain = np.exp(-chain[:, 0])[::10]
ki_chain = np.exp(-chain[:, 1])[::10]


# Identifier as to plot all of the data or not. 
ALL_DATA = 1
FIT_STRAIN = ['O2', 260]
IND_ONLY = 0
COLLAPSE = 1


# Compute the summary data
summary = data.groupby(['operator', 'repressors', 
            'IPTGuM']).agg(('mean','sem')).reset_index()

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
# Set up the figure. 
fig, ax = plt.subplots(2, 2, figsize=(5.5, 4))
for i, a in enumerate(ax.ravel()[:-1]):
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_xlabel('IPTG [µM]', style='italic')
    a.set_ylim([-0.015, 1.25])
    a.set_xlim([0, 1E4])

if COLLAPSE == 0:
    ax[1, 1].axis('off')
else:
    bohr_range = np.linspace(-10, 10)
    ax[1, 1].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, '-', color=colors['black'])
    ax[1, 1].set_xlabel('free energy [$k_BT$]', style='italic')
    ax[1, 1].set_ylabel('fold-change', style='italic')
    ax[1, 1].set_ylim([-0.15, 1.25])

# Add labels
for a in ax.ravel()[:-1]:
    a.set_ylabel('fold-change', style='italic')


# Define the axes and colors
axes = {'O1':ax[0, 0], 'O2':ax[0, 1], 'O3':ax[1, 0]}
rep_colors = {22:colors['dark_purple'], 60:colors['red'], 
              124:colors['black'], 260:colors['orange'],
              1220:colors['blue'], 1740:colors['green'] }
face_colors = {22:colors['light_purple'], 60:colors['light_red'], 
              124:colors['light_grey'], 260:colors['light_orange'],
              1220:colors['light_blue'], 1740:colors['light_green'] }

glyphs = {'O1': 's', 'O2': 'o', 'O3': '^'}
# Titles
for t, a in axes.items():
    a.set_title(f'{t}  ' + r'$\Delta\varepsilon_{RA} = %s\, k_BT$' %constants[t], style='italic', loc='left', y=0.95)

# ##############################################################################
# THEORY CURVES
# ##############################################################################
c_range = np.logspace(-2, 4, 500)
c_range[0] = 0
for g, d in data.groupby(['operator', 'repressors']):
    # Isolate the axis
    _ax = axes[g[0]]

    # Compute the architecture.
    arch = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]],
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=c_range)
    
    cred_region = np.zeros((2, len(c_range)))
    for i, c in enumerate(c_range):
        _arch = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]], 
            ka=ka_chain[::10], ki=ki_chain[::10], ep_ai=constants['ep_AI'], 
            effector_conc=c).fold_change()
        cred_region[:, i] = phd.stats.compute_hpd(_arch, 0.95)

    _ax.plot(c_range, arch.fold_change(), color=rep_colors[g[1]], label=int(g[1]),
    lw=0.75)
    _ax.fill_between(c_range, cred_region[0, :], cred_region[1, :],
                    color=rep_colors[g[1]], alpha=0.35, label='__nolegend__')

# ##############################################################################
# DATA 
# ##############################################################################
for g, d in summary.groupby(['operator', 'repressors']):
    # Get the axis
    _ax = axes[g[0]]

    # Define the colors. 
    if (g[0] == FIT_STRAIN[0]) & (g[1] == FIT_STRAIN[1]):
        face = 'w'
        plot = 1
    else:
        plot = 0
        face = face_colors[g[1]]
    
    # Plot the data. 
    if ALL_DATA == 0: 
        if plot == 1:
            _ax.errorbar(d['IPTGuM'], d['fold_change_A']['mean'], 
                    d['fold_change_A']['sem'], lw=1, capsize=1,
            fmt='.', ms=6, markeredgewidth=0.5, markerfacecolor=face, 
            color=rep_colors[g[1]], label='__legend__')
    else:
         _ax.errorbar(d['IPTGuM'], d['fold_change_A']['mean'], 
                     d['fold_change_A']['sem'], lw=1, capsize=1,
            fmt='.', ms=6, markeredgewidth=0.5, markerfacecolor=face, 
            color=rep_colors[g[1]], label='__nolegend__')

# ##############################################################################
# COLLAPSE
# ##############################################################################
if IND_ONLY == 0:
    ax[1,1].errorbar([], [], [], fmt='s', ms=2, color=colors['light_grey'], label='O1')
    ax[1,1].errorbar([], [], [], fmt='o', ms=2, color=colors['light_grey'], label='O2')
    ax[1,1].errorbar([], [], [], fmt='^', ms=2, color=colors['light_grey'], label='O3')
    ax[1,1].legend()
    for g, d in summary.groupby(['operator', 'repressors', 'IPTGuM']):
        # Compute the bohr parameter (point estimate)
        bohr = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]], ka=constants['Ka'],
                                       ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                       effector_conc=g[-1]).bohr_parameter()
        if (g[0] == 'O2') & (g[1]==260):
            face = 'white'
        else:
            face = face_colors[g[1]]
            label = '__nolegend__'
        ax[1, 1].errorbar(bohr, d['fold_change_A']['mean'], d['fold_change_A']['sem'],
                         fmt=glyphs[g[0]], color=rep_colors[g[1]], lw=0.5, 
                         capsize=1, ms=2, markeredgewidth=0.5, 
                         markerfacecolor=face, label='__nolegend__')

leg = ax[0, 0].legend(title='rep. / cell')
leg.get_title().set_fontsize(6)
plt.tight_layout()

if (ALL_DATA == 0) & (COLLAPSE == 0) & (IND_ONLY==1):
    plt.savefig('../figs/induction_fit_strain_only.pdf', facecolor='white')
elif (ALL_DATA == 1) & (COLLAPSE == 0) & (IND_ONLY==1):
    plt.savefig('../figs/induction_all_data.pdf', facecolor='white')
elif (ALL_DATA == 1) & (COLLAPSE == 1) & (IND_ONLY==1):
    plt.savefig('../figs/induction_all_data_collapse.pdf', facecolor='white')
elif (ALL_DATA == 1) & (COLLAPSE == 1) & (IND_ONLY==0):
    plt.savefig('../figs/induction_all_data_collapse_all_data.pdf', facecolor='white')
#%%



