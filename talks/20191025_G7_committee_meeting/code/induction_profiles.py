#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.stats
import phd.viz
import phd.thermo
import pickle
colors, color_list = phd.viz.phd_style()
constants = phd.thermo.load_constants()
title_bbox = dict(facecolor='none', edgecolor=colors['light_grey'], lw=0.1)
#%%
# Define the figure generation variables. 
ALL_DATA = False
ALL_PRED = True
C_RANGE = np.logspace(-2, 4, 200)

# Define the repressor colors. 
edge_reps = {22:colors['dark_purple'],  60:colors['dark_orange'],  
            124: colors['black'], 260: colors['dark_red'], 
            1220: colors['dark_blue'], 1740: colors['dark_green']}
fill_reps = {22:colors['light_purple'],  60:colors['light_orange'],   
            124: colors['light_grey'], 260: colors['light_red'], 
            1220: colors['light_blue'], 1740: colors['light_green']}
# %%
# Load the summarzed fold-change data
data = pd.read_csv('../data/RazoMejia2018_data.csv', comment='#')

# Summarize the data for display.
data = data[data['repressors'] > 0]
data['repressors'] *= 2
data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM': 'IPTGuM'},
            inplace=True)

summarized = data.groupby(['operator', 'repressors', 
                           'IPTGuM']).agg(('mean', 'sem')).reset_index()


# Load the MCMC samples from the induction paper for O2 R260
with open('../data/SI_I_O2_R260.pkl', 'rb') as f:
    out = pickle.load(f)
    fit_ka = np.exp(-out[:, 0][::10])
    fit_ki = np.exp(-out[:, 1][::10])
# %%
# Instantiate the figure canvas. 
fig, ax = plt.subplots(1, 3, figsize=(8, 3), sharey=True, sharex=True)

# Define the operator axes
_axes = {'O1':ax[0], 'O2':ax[1], 'O3':ax[2]}

# Format the axes and add labels. 
for a in ax:
    a.set_xscale('log')
    a.set_xlabel('IPTG [$\mu$M]')
ax[0].set_ylabel('fold-change')
ax[0].set_title('operator O1 (high affinity)', fontsize=8, color=colors['light_grey'], y=1.03, 
                bbox=title_bbox)
ax[1].set_title('operator O2 (moderate affinity)', fontsize=8, color=colors['light_grey'], y=1.03, 
                bbox=title_bbox)
ax[2].set_title('operator O3 (low affinity)', fontsize=8, color=colors['light_grey'], y=1.03, 
                bbox=title_bbox)

# Plot the theoretical predictions
if ALL_PRED:
    for g, d in data.groupby(['repressors', 'operator']):
        cred_region = np.zeros((2, len(C_RANGE)))
        for i, c in enumerate(C_RANGE):
            fc = phd.thermo.SimpleRepression(R=g[0], ep_r=constants[g[1]], 
                                             ka=fit_ka, ki=fit_ki, 
                                             ep_ai=constants['ep_AI'],
                                             effector_conc=c).fold_change()
            cred_region[:, i] = phd.stats.compute_hpd(fc, 0.95) 
        _axes[g[1]].fill_between(C_RANGE, cred_region[0, :], cred_region[1, :],
                                color=fill_reps[g[0]], alpha=0.75, label=g[0])

# Add the legend now. 
leg = ax[0].legend(fontsize=6, title='repressors', loc='upper left')


# Plot the data. 
for g, d in summarized.groupby(['repressors', 'operator']):
    if (ALL_DATA==False):
        if (g[0] != 260) | (g[1] != 'O2'):
            pass
        else:
            _ax = _axes[g[1]]
            _ax.errorbar(d['IPTGuM'], d['fold_change']['mean'], 
                    d['fold_change']['sem'], capsize=1, lw=0.75, color=edge_reps[g[0]],
                    fmt='.', ms=6, markerfacecolor=fill_reps[g[0]], markeredgewidth=0.75,
                    label='__nolegend__')
    else:
       _ax = _axes[g[1]]
       _ax.errorbar(d['IPTGuM'], d['fold_change']['mean'], 
            d['fold_change']['sem'], capsize=1, lw=0.75, color=edge_reps[g[0]],
            fmt='.', ms=6, markerfacecolor=fill_reps[g[0]], markeredgewidth=0.75,
            label='__nolegend__')

leg.get_title().set_fontsize(6)


# Save it
for a in ax:
    a.set_ylim([-0.05, 1.05])

plt.tight_layout()
if ALL_DATA == False:
    name = 'fit_data'
else:
    name = 'all_data'


plt.savefig(f'../figs/{name}_induction_profiles.pdf', bbox_inches='tight',
            facecolor='white')
# %%


# %%


# %%
