# -*- coding: utf-8 -*-
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import phd.thermo
import phd.viz
colors, color_list = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load in the data sets
data = pd.read_csv('../data/summarized_data.csv')
data = data[data['class']=='DNA']
stats = pd.read_csv('../data/DNA_binding_energy_summary.csv')
stats = stats[stats['repressors']==260]
bohr = pd.read_csv('../data/empirical_F_statistics.csv')
empirical_bohr = bohr[bohr['class']=='DNA']

# Define some plotting constants. 
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
bohr_range = np.linspace(-8, 8, 200)
F = (1 + np.exp(-bohr_range))**-1

# Define the colors
rep_colors = {60:colors['red'], 
              124:colors['black'], 260:colors['orange'],
              1220:colors['blue']}
face_colors = {22:colors['light_purple'], 60:colors['light_red'], 
              124:colors['light_grey'], 260:colors['light_orange'],
              1220:colors['light_blue'], 1740:colors['light_green'] }
ALL_DATA = 1
title_bbox = dict(facecolor='none', edgecolor=colors['light_grey'], lw=0.25)
# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 3, figsize=(6, 4), dpi=150)

for i in range(3):
    ax[0, i].set_xscale('symlog', linthreshx=0.006)
    ax[-1, i].set_xscale('symlog', linthreshx=0.006)
    ax[0, i].set_ylim([-0.2, 1.2])
    ax[0, i].set_xlim([-0.001, 1E4])
    ax[-1, i].set_xlim([-0.001, 1E4])
    ax[-1, i].set_ylim([-8, 8])
    ax[0, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    ax[-1, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    for j in range(2):
        if i!=2:
            ax[i, j+1].set_yticklabels([])

    # Add labels   
    ax[0, i].set_xlabel("IPTG [$\mu$M]")
    ax[-1, i].set_xlabel("IPTG [$\mu$M]")

# Add ylabels
ax[0, 0].set_ylabel('fold-change')
ax[1, 0].set_ylabel('$\Delta F$ [$k_BT$]')

# Define the axes
axes = {'Q21M':0, 'Y20I':1, 'Q21A':2}
corr_ax = {'Y17I':1, 'Q18M':0, 'Q18A':2}

# Add titles
for m, a in corr_ax.items():
    ax[0, a].set_title(m, y=1.03, bbox=title_bbox) 

# ##############################################################################
# GUIDE CURVES FOR ∆F
# ##############################################################################
for i in range(3):
    ax[-1, i].hlines(0, -0.01, 1E4, 'k', linestyle=':', lw=0.75)

##############################################################################
#FOLD-CHANGE CURVES
##############################################################################
for r, cor in rep_colors.items():
    for m, a in axes.items():
        _stats = stats[(stats['mutant']==m) & (stats['parameter']=='ep_RA')]
        _c, _ep= np.meshgrid(c_range, _stats[['hpd_min', 'hpd_max']].values)
        arch = phd.thermo.SimpleRepression(R=r, ep_r=_ep, ka=constants['Ka'],
                                           ki=constants['Ki'], 
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=_c).fold_change()
        ax[0, a].fill_between(c_range, arch[0, :], arch[1, :], color=cor,
                             alpha=0.4) 

# # ##############################################################################
# # COLLAPSE CURVES
# # ##############################################################################
# for i in range(3):
#     ax[1, i].plot(bohr_range, F, 'k-', lw=0.75)

##############################################################################
# FREE ENERGY PREDICTIONS
# ##############################################################################
for m, a in axes.items():
    _stats = stats[(stats['mutant']==m) & 
                 (stats['parameter']=='ep_RA')][['hpd_min', 'hpd_max']].values[0]
    ax[-1, a].fill_between(c_range, _stats[0] - constants['O2'], 
                         _stats[1] - constants['O2'], color=rep_colors[260], 
                    alpha=0.4)

# ##############################################################################
# FOLD-CHANGE DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'repressors']):
    if g[1] == 260:
        face = 'w'
        alpha = 1
        legend = int(g[1])
    else:
        face = face_colors[g[1]]
        if ALL_DATA == 0:
            alpha = 0 
            legend = '__nolegend__'
        else:
            alpha = 1
            legend = int(g[1]) 

    ax[0, axes[g[0]]].errorbar(d['IPTGuM'], d['mean'], d['sem'], fmt='.',
                               color=rep_colors[int(g[1])], markerfacecolor=face,
                               ms=6, markeredgewidth=0.5, capsize=1, lw=0.25, linestyle='none',
                               label=legend, alpha=alpha)

# # ##############################################################################
# # COLLAPSE DATA
# # ##############################################################################
# for g, d in data.groupby(['mutant', 'repressors']):
#     ep_r = stats[(stats['mutant']==g[0]) & 
#             (stats['parameter']=='ep_RA')]['median'].values[0]
#     bohr = phd.thermo.SimpleRepression(R=g[1], ep_r=ep_r, ka=constants['Ka'],
#                                        ki=constants['Ki'], 
#                                        ep_ai=constants['ep_AI'], 
#                                        effector_conc=d['IPTGuM']).bohr_parameter()
#     if g[1] == 260:
#         face = 'w'
#     else:
#         face = rep_colors[g[1]]
#     ax[1, axes[g[0]]].errorbar(bohr, d['mean'], d['sem'], fmt='.', 
#                                linestyle='none', lw=1, capsize=1, ms=5, markeredgewidth=0.5,
#                                color=rep_colors[g[1]], markerfacecolor=face)

# ##############################################################################
# INFERRED F 
# ##############################################################################
for g, d in empirical_bohr.groupby(['mutant', 'repressors', 'IPTGuM']):
    _param = d[d['parameter']=='delta_bohr']
    mu = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (mu < sig) | (1 - mu < sig):
        color = 'slategray'
        alpha = 0
        lw = 0
        fmt = 'x'
    else:
        color = rep_colors[g[1]]
        alpha = 1 
        lw = 0.75
        fmt = '.'
    if g[1] == 260:
        face = 'w'
        zorder=1000
    elif fmt == 'x':
        zorder = 1
        face=face_colors[g[1]]
    else:
        if ALL_DATA == 1:
            alpha = 1
        else: 
            alpha = 0
        face = face_colors[g[1]]
        zorder=100
    _ax = ax[-1, axes[g[0]]]
    _ax.plot(_param['IPTGuM'], _param['median'], marker=fmt, linestyle='none', 
        color=color, markerfacecolor=face, alpha=alpha, ms=8, zorder=zorder, 
        markeredgewidth=0.5)
    _ax.vlines(_param['IPTGuM'], _param['hpd_min'], _param['hpd_max'], 
            lw=lw, color=color,  alpha=alpha, zorder=zorder)

# ##############################################################################
#  LEGEND INFORMATION
# ##############################################################################
leg = ax[0, 0].legend(title='rep. / cell',handletextpad=0.1)
leg.get_title().set_fontsize(6)
plt.subplots_adjust(hspace=0.5, wspace=0.1)
if ALL_DATA == 1:
    plt.savefig('../figs/DNA_muts_all_data_fit.pdf', bbox_inches='tight', facecolor='white')
# else:
#     plt.savefig('../figs/DNA_muts_fit_strain.pdf', bbox_inches='tight')



#%%



# %%
