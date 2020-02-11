#%% 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
import pickle
colors, palette = phd.viz.phd_style()
show_data = False
constants = phd.thermo.load_constants()

#%% Load the titration data and prune
data = pd.read_csv('../data/RazoMejia2018_data.csv', comment='#')
data = data[(data['repressors'] > 0)]
data['repressors'] *= 2
data.rename(columns={'fold_change_A':'fold_change'}, inplace=True)

# Comput the summary of the data. 
data = data.groupby(['repressors', 'operator', 'IPTG_uM']).agg(
                    ('mean', 'sem')).reset_index()

# %%
# Load the MCMC samples
with open('../data/SI_I_O2_R260.pkl', 'rb') as file:
    chains = pickle.load(file)
    chain_df = pd.DataFrame(chains, columns=['ep_a', 'ep_i', 'sigma'])
    chain_df['Ka'] = np.exp(-chain_df['ep_a'])
    chain_df['Ki'] = np.exp(-chain_df['ep_i'])

# Define the inducer concentration range
c_range = np.logspace(-2, 4, 200)

band_colors = ['red', 'brown', 'green', 'orange', 'purple', 'blue']
reps = [22, 60, 124, 260, 1220, 1740]
rep_edge_colors = {r:colors[f'dark_{c}'] for r, c in zip(reps, band_colors)}
rep_fill_colors = {r:colors[f'light_{c}'] for r, c in zip(reps, band_colors)}
# %%
# Instantiate the figure canvas. 
fig, ax = plt.subplots(1, 3, figsize=(10, 3), dpi=100, sharex=True, sharey=True)
for a in ax:
    a.set_xscale('log')
    a.set_xlabel('IPTG concentration [µM]', fontsize=10)
    a.set_xlim([c_range[0], c_range[-1]])
    a.set_ylim([-0.05, 1.15])
    a.tick_params(labelsize=9)
ax[0].set_ylabel('fold-change in gene expression', fontsize=10)

# Set titles
for i, o in enumerate(['O1', 'O2', 'O3']):
    phd.viz.titlebox(ax[i], f'operator {o} | ' + r'$\Delta\varepsilon_{RA}=%s\, k_BT$' % constants[o],
                     color=colors['black'], bgcolor=colors['grey'], size=10,
                     pad=0.07)

# Plot the cred regions
for i , o in enumerate(['O1', 'O2', 'O3']):
    for j, r in enumerate([22, 60, 124, 260, 1220, 1740]):
        cred_region = np.zeros((2, len(c_range)))
        for k, c in enumerate(c_range):
            fc = phd.thermo.SimpleRepression(R=r, ep_r=constants[o], 
                    effector_conc=c, ka=chain_df['Ka'].values[::10],
                    ki=chain_df['Ki'].values[::10], ep_ai=4.5).fold_change()
            cred_region[:, k] = phd.stats.compute_hpd(fc, 0.95)
        ax[i].fill_between(c_range, cred_region[0, :], cred_region[1, :],
                color=colors[band_colors[j]], alpha=0.5)

# Plot the data. 
op_ax = {'O1':ax[0], 'O2':ax[1], 'O3':ax[2]}
for g, d in data.groupby(['operator', 'repressors']):
    _ax = op_ax[g[0]]
    if (g[0] == 'O2') & (g[1]==260):
        face='w'
        _ax.errorbar(d['IPTG_uM'], d['fold_change']['mean'], 
            yerr=d['fold_change']['sem'], fmt='.', ms=10, 
            markerfacecolor=face,
            markeredgecolor=rep_edge_colors[g[1]],
            color=rep_edge_colors[g[1]], label=int(g[1]))
    else:
        face = rep_fill_colors[g[1]]
        if show_data == True:
            _ax.errorbar(d['IPTG_uM'], d['fold_change']['mean'], 
                yerr=d['fold_change']['sem'], fmt='.', ms=10, 
                markerfacecolor=face,
                markeredgecolor=rep_edge_colors[g[1]],
                color=rep_edge_colors[g[1]], label=int(g[1]))

leg = ax[1].legend(title='repressors', fontsize=9)        
leg.get_title().set_fontsize(9)
if show_data == True:
    data_name = 'data'
else:
    data_name = 'nodata'

plt.subplots_adjust(wspace=0.1)
plt.savefig(f'../figs/induction_{data_name}_predictions.pdf', bbox_inches='tight',
            dpi=100, facecolor=None)
# %%
