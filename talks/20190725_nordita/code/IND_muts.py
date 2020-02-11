# -*- coding; utf-8 -*- 
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.thermo
import phd.viz
import seaborn as sns
constants = phd.thermo.load_constants()
colors = phd.viz.phd_style()
_colors = sns.color_palette('magma', n_colors=3)
# Load and restrict the various data sets
data = pd.read_csv('../data/summarized_data.csv')
data = data[data['class']=='IND'].copy()
kaki_only_stats = pd.read_csv('../data/KaKi_only_summary.csv')
kaki_only_stats = kaki_only_stats[kaki_only_stats['operator']=='O2'].copy()
kaki_epAI_stats = pd.read_csv('../data/kaki_epAI_summary.csv')
kaki_epAI_stats = kaki_epAI_stats[kaki_epAI_stats['operator']=='O2']
kaki_epAI_samps = pd.read_csv('../data/KaKi_epAI_samples.csv')
kaki_epAI_samps = kaki_epAI_samps[kaki_epAI_samps['operator']=='O2'].copy()
bohr = pd.read_csv('../data/empirical_F_statistics.csv')
bohr = bohr[bohr['class']=='IND'].copy()

# Define constants for plotting 
c_range = np.logspace(-3, 4, 200)
c_range[0] = 0
bohr_range = np.linspace(-8, 8, 200)
F = (1 + np.exp(-bohr_range))**-1

# ##############################################################################
# PLOTTING IDENTIFIERS
# ##############################################################################
KAKI_ONLY = 0
KAKI_epAI = 1
ALL_DATA = 1

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 4, figsize=(7.5, 4))

# Add appropriate scaling
for i in range(4):
    ax[0, i].set_xscale('symlog', linthreshx=0.01)
    ax[-1, i].set_xscale('symlog', linthreshx=0.01)
    ax[-1, i].set_ylim([-6, 5])
    ax[0, i].set_ylim([-0.2, 1.2])

    ax[0, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4]) 
    ax[0, i].set_yticks([0, 0.5, 1])
    ax[-1, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4]) 
    ax[-1, i].set_xlim([-0.001, 1E4])


# Define the axes
axes = {'F164T':3, 'Q294V':2, 'Q294K':1, 'Q294R':0}
axes_corr = {'F161T':3, 'Q291V':2, 'Q291K':1, 'Q291R':0}
op_colors = {'O1':_colors[0], 'O2':_colors[1], 'O3':_colors[2]}

# Add labels
for m, a in axes_corr.items():
    ax[0, a].set_title(m, fontsize=8, loc='left', style='italic',
                       y=0.95)
ax[0, 0].set_ylabel('fold-change', fontsize=8, style='italic')
ax[1, 0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=8, style='italic')

for i in range(4):
    ax[0, i].set_xlabel('IPTG [µM]', fontsize=8, style='italic')                   
    ax[-1, i].set_xlabel('IPTG [µM]', fontsize=8, style='italic')                   


for i in range(2):
    for j in range(3):
        ax[i, j+1].set_yticklabels([])

# ##############################################################################
# FOLD CHANGE CURVES 
# ##############################################################################
for i, o in enumerate(('O1', 'O2', 'O3')):
    ep_r = constants[o]
    for m, a in axes.items():
       # Plot the kaki only fits. 
       _kaki = kaki_only_stats[kaki_only_stats['mutant']==m]
       ka = _kaki[_kaki['parameter']=='Ka']['median'].values[0]
       ki = _kaki[_kaki['parameter']=='Ki']['median'].values[0]
       arch = phd.thermo.SimpleRepression(R=260, ep_r=ep_r, ka=ka, ki=ki, 
                ep_ai=4.5, effector_conc=c_range).fold_change()

       if KAKI_ONLY == 1:
            ax[0, a].plot(c_range, arch, '--', color=op_colors[o], lw=1)  

       if KAKI_epAI == 1:
             # Plot the credible regions for the kaki and epAI fits.
            _kaki = kaki_epAI_samps[kaki_epAI_samps['mutant']==m]
            cred_region = np.zeros((2, len(c_range)))
            for j, c in enumerate(c_range):
                 arch = phd.thermo.SimpleRepression(R=260, ep_r=ep_r, ka=_kaki['Ka'], 
                         ki=_kaki['Ki'],  ep_ai=_kaki['ep_AI'], effector_conc=c).fold_change()
                 cred_region[:, j] = phd.stats.compute_hpd(arch, 0.95)
    
            ax[0, a].fill_between(c_range, cred_region[0, :], cred_region[1, :], 
                              color=op_colors[o], alpha=0.4)

# ##############################################################################
#  FREE ENERGY PREDICTIONS
# ##############################################################################
ref_pact  = phd.thermo.MWC(effector_conc=c_range, ka=constants['Ka'],
                           ki=constants['Ki'], ep_ai=constants['ep_AI']).pact()

if KAKI_ONLY == 1:
    for g, d in data.groupby(['mutant']):
        _stats = kaki_only_stats[kaki_only_stats['mutant']==g]  
        ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
        ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
        pact = phd.thermo.MWC(effector_conc=c_range, ka=ka, ki=ki, 
                             ep_ai=constants['ep_AI']).pact()
        delF = -np.log(pact / ref_pact)
        ax[-1, axes[g]].plot(c_range, delF, '--', color=op_colors[o], 
                            label='__nolegend__')

if KAKI_epAI == 1:
    for g, d in data.groupby(['mutant']):
        _kaki_samps = kaki_epAI_samps[kaki_epAI_samps['mutant']==g]
        cred_region = np.zeros((2, len(c_range)))
        for i, c in enumerate(c_range):
            pact = phd.thermo.MWC(effector_conc=c, ka=_kaki_samps['Ka'],
                                  ki=_kaki_samps['Ki'], ep_ai=_kaki_samps['ep_AI']).pact()
            delF = -np.log(pact / ref_pact[i])
            cred_region[:, i] = phd.stats.compute_hpd(delF, 0.95)

        # Plot!
        ax[-1, axes[g]].fill_between(c_range, cred_region[0, :], cred_region[1, :],
                                    color=op_colors['O2'], alpha=0.4)

# ##############################################################################
# FOLD-CHANGE DATA 
# ##############################################################################
for g, d in data.groupby(['mutant', 'operator']):
    _ax = ax[0, axes[g[0]]]
    if g[1] == 'O2':
        face = 'w'
    else:
        face = op_colors[g[1]]

    if ALL_DATA == 0:
        if g[1] == 'O2':
            _ax.errorbar(d['IPTGuM'], d['mean'], d['sem'], fmt='.', 
            markerfacecolor=face, linestyle='none', color=op_colors[g[1]], capsize=1,
            label=g[1], markeredgewidth=0.5, markersize=5, lw=0.5)
    else:
        _ax.errorbar(d['IPTGuM'], d['mean'], d['sem'], fmt='.', 
            markerfacecolor=face, linestyle='none', color=op_colors[g[1]], capsize=1,
            label=g[1], markeredgewidth=0.5, markersize=5, lw=0.5)


# ##############################################################################
# DELTA F DATA
# ##############################################################################
for g, d in bohr.groupby(['mutant', 'operator', 'IPTGuM']):
    _ax = ax[1, axes[g[0]]]
    _param = d[d['parameter']=='delta_bohr']
    mu = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (mu < sig) | (1 - mu < sig):
        color = 'slategray'
        alpha = 0
        lw = 0
        fmt = 'x'
    else:
        color = op_colors[g[1]]
        alpha = 1 
        lw = 0.5
        fmt = '.'
    if g[1] == 'O2':
        face = 'w'
        zorder=1000
    else:
        face = color
        zorder=100
    if g[-1] == 0:
        cap_min = -0.1
        cap_max = 0.001
    else:
        cap_min = g[-1] * 0.8
        cap_max = g[-1] * 1.2
    if ALL_DATA == 0:
        if g[1] == 'O2':
            _ax.plot(_param['IPTGuM'], _param['median'], linestyle='none', marker=fmt, color=color, 
                    markerfacecolor=face , alpha=alpha, ms=5, zorder=zorder,
                    markeredgewidth=0.5)
            _ax.vlines(_param['IPTGuM'], _param['hpd_min'], _param['hpd_max'], color=color,
                    lw=lw, zorder=zorder)
            _ax.hlines(_param['hpd_min'], cap_min, cap_max, lw=lw, zorder=zorder, 
                       color=color)
            _ax.hlines(_param['hpd_max'], cap_min, cap_max, lw=lw, zorder=zorder, 
               color=color)
    else:
        _ax.plot(_param['IPTGuM'], _param['median'], linestyle='none', marker=fmt, color=color, 
                    markerfacecolor=face , alpha=alpha, ms=5, zorder=zorder,
                    markeredgewidth=0.5)
        _ax.vlines(_param['IPTGuM'], _param['hpd_min'], _param['hpd_max'], color=color,
                    lw=lw, zorder=zorder)
        _ax.hlines(_param['hpd_min'], cap_min, cap_max, lw=lw, zorder=zorder, 
                       color=color)
        _ax.hlines(_param['hpd_max'], cap_min, cap_max, lw=lw, zorder=zorder, 
               color=color)

# ##############################################################################
# LEGEND INFORMATION
# ##############################################################################
leg = ax[0, 0].legend(fontsize=5, ncol=3, columnspacing=0.001, handletextpad=0.01) #, title='operator'b)
plt.subplots_adjust(wspace=0.15, hspace=0.5)


if (ALL_DATA==0) & (KAKI_ONLY==1) & (KAKI_epAI==0):
    plt.savefig('../figs/IND_fit_strain_kaki.pdf', bbox_inches='tight', facecolor='white')
if (ALL_DATA==0) & (KAKI_ONLY==0) & (KAKI_epAI==1):
    plt.savefig('../figs/IND_fit_strain_kaki_epAI.pdf', bbox_inches='tight', facecolor='white')
if (ALL_DATA==0) & (KAKI_ONLY==1) & (KAKI_epAI==1):
    plt.savefig('../figs/IND_fit_strain_both.pdf', bbox_inches='tight', facecolor='white')
if (ALL_DATA==1) & (KAKI_ONLY==0) & (KAKI_epAI==1):
    plt.savefig('../figs/IND_all_data_kaki_epAI.pdf', bbox_inches='tight', facecolor='white')
if (ALL_DATA==1) & (KAKI_ONLY==1) & (KAKI_epAI==0):
    plt.savefig('../figs/IND_all_data_kaki.pdf', bbox_inches='tight', facecolor='white')
if (ALL_DATA==1) & (KAKI_ONLY==1) & (KAKI_epAI==1):
    plt.savefig('../figs/IND_all_data.pdf', bbox_inches='tight', facecolor='white')
if (ALL_DATA==1) &  (KAKI_ONLY==0) & (KAKI_epAI==0):
    plt.savefig('../figs/IND_all_data_only.pdf', bbox_inches='tight', facecolor='white')
if (ALL_DATA==0) &  (KAKI_ONLY==0) & (KAKI_epAI==0):
    plt.savefig('../figs/IND_fit_only.pdf', bbox_inches='tight', facecolor='white')


#%%
