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
colors, palette = phd.viz.phd_style()
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
ALL_DATA = 0
FIT_STRAIN = ['O2', 260]

# Compute the summary data
summary = data.groupby(['operator', 'repressors', 
            'IPTGuM']).agg(('mean','sem')).reset_index()

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
# Set up the figure. 
fig, ax = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(6, 2))
for a in ax:
    a.set_xscale('symlog', linthreshx=1E-2)
    a.set_xlabel('IPTG [µM]')
    a.set_ylim([-0.1, 1.25])

# Add labels
ax[0].set_ylabel('fold-change')


# Define the axes and colors
axes = {'O1':ax[0], 'O2':ax[1], 'O3':ax[2]}
rep_colors = {22:colors['dark_purple'], 60:colors['red'], 
              124:colors['black'], 260:colors['orange'],
              1220:colors['blue'], 1740:colors['green'] }
face_colors = {22:colors['light_purple'], 60:colors['light_red'], 
              124:colors['light_grey'], 260:colors['light_orange'],
              1220:colors['light_blue'], 1740:colors['light_green'] }

# Titles
for t, a in axes.items():
    a.set_title(f'{t}  ' + r'$\Delta\varepsilon_{RA} = %s\, k_BT$' %constants[t], y=1.04, 
            backgroundcolor=colors['gray'])

# # ##############################################################################
# # THEORY CURVES
# # ##############################################################################
# c_range = np.logspace(-2, 4, 500)
# c_range[0] = 0
# for g, d in data.groupby(['operator', 'repressors']):
#     # Isolate the axis
#     _ax = axes[g[0]]

#     # Compute the architecture.
#     arch = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]],
#                                        ka=constants['Ka'], ki=constants['Ki'],
#                                        ep_ai=constants['ep_AI'], 
#                                        effector_conc=c_range)
    
#     cred_region = np.zeros((2, len(c_range)))
#     for i, c in enumerate(c_range):
#         _arch = phd.thermo.SimpleRepression(R=g[1], ep_r=constants[g[0]], 
#             ka=ka_chain, ki=ki_chain, ep_ai=constants['ep_AI'], 
#             effector_conc=c).fold_change()
#         cred_region[:, i] = phd.stats.compute_hpd(_arch, 0.95)

#     _ax.plot(c_range, arch.fold_change(), color=rep_colors[g[1]], label=int(g[1]),
#     lw=0.75)
#     _ax.fill_between(c_range, cred_region[0, :], cred_region[1, :],
#                     color=rep_colors[g[1]], alpha=0.35, label='__nolegend__')

# # ##############################################################################
# # DATA 
# # ##############################################################################
# for g, d in summary.groupby(['operator', 'repressors']):
#     # Get the axis
#     _ax = axes[g[0]]

#     # Define the colors. 
#     if (g[0] == FIT_STRAIN[0]) & (g[1] == FIT_STRAIN[1]):
#         face = 'w'
#         plot = 1
#     else:
#         plot = 0
#         face = face_colors[g[1]]
    
#     # Plot the data. 
#     if ALL_DATA == 0: 
#         if plot == 1:
#             _ax.errorbar(d['IPTGuM'], d['fold_change_A']['mean'], 
#                     d['fold_change_A']['sem'], lw=1, capsize=1,
#             fmt='.', ms=6, markeredgewidth=0.5, markerfacecolor=face, 
#             color=rep_colors[g[1]], label='__legend__')
#     else:
#          _ax.errorbar(d['IPTGuM'], d['fold_change_A']['mean'], 
#                      d['fold_change_A']['sem'], lw=1, capsize=1,
#             fmt='.', ms=6, markeredgewidth=0.5, markerfacecolor=face, 
#             color=rep_colors[g[1]], label='__nolegend__')

# leg = ax[1].legend(title='rep. / cell', fontsize=5)
# leg.get_title().set_fontsize(5)

# plt.tight_layout()
# if ALL_DATA == 0:
#     plt.savefig('../figs/induction_fit_strain_only.pdf')
# else:
#     plt.savefig('../figs/induction_all_data.pdf')


# ##############################################################################
#  PROPERTIES
# ##############################################################################
ALL_DATA = 1
# Load the MCMC chains to draw credible regions. 
with open('../data/main_text_KaKi.pkl', 'rb') as pkl:
    chain = pickle.load(pkl)
ka_chain = np.exp(-chain[:, 0])[::10]
ki_chain = np.exp(-chain[:, 1])[::10]

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
#%%
fig, ax = plt.subplots(2, 3, figsize=(7, 4.5), dpi=100)
ax[0,0].axis(False)
ax = [ax[0, 1], ax[0, 2], ax[1, 0], ax[1, 1], ax[1, 2]]

for a in ax:
    a.set_xscale('log')
    a.set_xlabel('repressors per cell', fontsize=10)
    a.tick_params(labelsize=9)
ax[0].set_yscale('log')
ax[3].set_yscale('log')

# Add specific y axis labels.
ax[0].set_ylabel('leakiness', fontsize=10)
ax[1].set_ylabel('saturation', fontsize=10)
ax[2].set_ylabel('dynamic range', fontsize=10)
ax[3].set_ylabel('$[EC_{50}]$ [µM]', fontsize=10)
ax[4].set_ylabel('effective Hill coefficient', fontsize=10)

# Fix the axis limits
ax[0].set_ylim([1E-4, 1.2])
ax[1].set_ylim([1E-4, 1.02])
ax[2].set_ylim([1E-4, 1.02])
ax[3].set_ylim([1E-2, 500])
ax[4].set_ylim([1.15, 1.95])

# Define teh operator colors.
op_colors = {'O1':colors['purple'], 'O2':colors['orange'], 'O3': colors['blue']}
op_fill_colors = {'O1':colors['pale_purple'], 'O2':colors['pale_orange'], 
                  'O3': colors['pale_blue']}

# Define the property axes 
axes = {'leakiness':ax[0], 'saturation':ax[1], 'dynamic_range':ax[2],
        'EC50':ax[3], 'effective_hill':ax[4]}
# # ##############################################################################
# #  THEORY CURVES
# # ##############################################################################
rep_range = np.logspace(0, 4, 200)

for i, o in enumerate(['O1', 'O2', 'O3']):
    # Compute the modes
    arch = phd.thermo.SimpleRepression(R=rep_range, ep_r=constants[o],
                               ka=constants['Ka'], ki=constants['Ki'],
                               effector_conc=0, 
                               ep_ai=constants['ep_AI']).compute_properties()
    for k, v in arch.items():
        axes[k].plot(rep_range, v, color=op_colors[o], label=o)

    # Credible regions
    _ = np.zeros((2, len(rep_range)))
    cred_regions = {'saturation': _.copy(),
                    'dynamic_range':_.copy(), 'EC50':_.copy(),
                    'effective_hill':_.copy()}
    for j, r in enumerate(rep_range):
        arch = phd.thermo.SimpleRepression(R=r,ep_r=constants[o], 
                                           ka=ka_chain, ki=ki_chain,
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=0).compute_properties()
        for k, v in arch.items():
            if k != 'leakiness':
                cred_regions[k][:, j] = phd.stats.compute_hpd(v, 0.95)
    for k, v in cred_regions.items():
        axes[k].fill_between(rep_range, v[0, :], v[1, :], color=op_colors[o],
                            alpha=0.3, label='__nolegend__')

# # ##############################################################################
# # DATA
# # ##############################################################################
if ALL_DATA == 1:
    for g, d in data.groupby(['repressors', 'operator']):
        leakiness = d[d['IPTGuM']==0]['fold_change_A'].values
        sat =d[d['IPTGuM']==d['IPTGuM'].max()]['fold_change_A'].values
        try:
            dyn_rng = sat - leakiness  
        except:
            dyn_rng = sat - leakiness[:-1]
        ax[0].errorbar(g[0], np.mean(leakiness), 
                    np.std(leakiness)/len(leakiness), lw=0.75, capsize=1, 
                    color=op_colors[g[-1]], label='__nolegend__', linestyle='none',
                    ms=9, fmt='.', markerfacecolor=op_fill_colors[g[1]]) 
        ax[1].errorbar(g[0], np.mean(sat), np.std(sat)/len(sat), lw=0.75, 
                    capsize=1, color=op_colors[g[-1]], label='__nolegend__',
                    linestyle='none', ms=9, fmt='.', markerfacecolor=op_fill_colors[g[1]]) 
        ax[2].errorbar(g[0], np.mean(dyn_rng), np.std(dyn_rng)/len(dyn_rng), 
                          lw=0.75, capsize=1, color=op_colors[g[1]], 
                          label='__nolegend__', linestyle='none', ms=9, fmt='.',
                          markerfacecolor=op_fill_colors[g[1]]) 

        # Inferred points
        with open(f'../data/SI_I_{g[-1]}_R{int(g[0])}.pkl', 'rb') as pkl:
            chain = pickle.load(pkl)
        ka_chain = np.exp(-chain[:,0])
        ki_chain = np.exp(-chain[:,1])
        ka_median = np.median(ka_chain)
        ki_median = np.median(ki_chain)

        # Compute the inferred EC50 and effective hill
        arch = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]],
                                           ka=ka_median, ki=ki_median, 
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=0).compute_properties()
        ax[3].plot(g[0], arch['EC50'], 's', color=op_colors[g[1]], ms=4,
        markerfacecolor=op_fill_colors[g[1]], markeredgewidth=1)
        ax[4].plot(g[0], arch['effective_hill'], 's', color=op_colors[g[1]], 
                   ms=4, markerfacecolor=op_fill_colors[g[1]], markeredgewidth=1)

        # Compute the credible regions:
        arch = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]],
                                           ka=ka_chain, ki=ki_chain, 
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=0).compute_properties()
        ec50_cred = phd.stats.compute_hpd(arch['EC50'], 0.95)
        hill_cred = phd.stats.compute_hpd(arch['effective_hill'], 0.95)
    
        # Plot the vlines. 
        ax[3].vlines(g[0], ec50_cred[0], ec50_cred[1], lw=0.75,
                    color=op_colors[g[1]])
        ax[4].vlines(g[0], hill_cred[0], hill_cred[1], lw=0.75,
                    color=op_colors[g[1]])

# # ##############################################################################
# # LEGEND AND SAVING
# # ##############################################################################
# leg = ax[0].legend(title='operator')
# leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../figs/properties_data.pdf', bbox_inches='tight')


# %%
