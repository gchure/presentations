#%% 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
import pickle
colors, palette = phd.viz.phd_style()

show_data = True
show_fit = True
show_predictions = True
# %%
#%% Load the titration data and prune
data = pd.read_csv('../data/RazoMejia2018_data.csv', comment='#')
data = data[(data['repressors']==130) & (data['operator']=='O2')]
data['repressors'] *= 2
data.rename(columns={'fold_change_A':'fold_change'}, inplace=True)

# Comput the summary of the data. 
data = data.groupby(['IPTG_uM']).agg(('mean', 'sem')).reset_index()

# %%
# Load the MCMC samples
with open('../data/SI_I_O2_R260.pkl', 'rb') as file:
    chains = pickle.load(file)
    chain_df = pd.DataFrame(chains, columns=['ep_a', 'ep_i', 'sigma'])
    chain_df['Ka'] = np.exp(-chain_df['ep_a'])
    chain_df['Ki'] = np.exp(-chain_df['ep_i'])

# Define the inducer concentration range
c_range = np.logspace(-2, 4, 200)
# %%
# Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=100)
ax.set_xscale('log')
ax.set_xlabel('IPTG concentration [µM]', fontsize=10)
ax.set_ylabel('fold-change in gene expression', fontsize=10)
ax.set_ylim([-0.05, 1.15])
ax.set_xlim([1E-2, 1E4])
ax.tick_params(labelsize=9)
phd.viz.titlebox(ax, r'operator O2   |    $\Delta\varepsilon_{RA} = -13.9\, k_BT$',
                 color=colors['black'], bgcolor=colors['grey'], pad=0.05,
                 size=10)
if show_data == True:
    ax.errorbar(data['IPTG_uM'], data['fold_change']['mean'], 
               yerr=data['fold_change']['sem'], fmt='.', linestyle='none',
               linewidth=1.5, markersize=10, markeredgecolor=colors['orange'],
               markerfacecolor='white', markeredgewidth=1, 
               color=colors['orange'], label=int(260))

if show_fit == True:
    cred_region  = np.zeros((2, len(c_range)))
    for i, c in enumerate(c_range):
        fc = phd.thermo.SimpleRepression(R=260, ep_r=-13.9, effector_conc=c,
                                    ka=chain_df['Ka'].values, 
                                    ki=chain_df['Ki'].values,
                                    ep_ai=4.5).fold_change()
        cred_region[:, i] = phd.stats.compute_hpd(fc, 0.95)

    ax.fill_between(c_range, cred_region[0, :], cred_region[1, :], 
                    color=colors['light_orange'], label='__nolegend__', alpha=0.5)

if show_predictions == True:
    fill_colors = ['red', 'brown', 'green', 'purple', 'blue']
    for i, r in enumerate([22, 60, 124, 1220, 1740]):
        cred_region  = np.zeros((2, len(c_range)))
        for j, c in enumerate(c_range):
            fc = phd.thermo.SimpleRepression(R=r, ep_r=-13.9, effector_conc=c,
                                            ka=chain_df['Ka'].values[::10], 
                                            ki=chain_df['Ki'].values[::10],
                                            ep_ai=4.5).fold_change()
            cred_region[:, j] = phd.stats.compute_hpd(fc, 0.95)       
        ax.fill_between(c_range, cred_region[0, :], cred_region[1, :], 
                        color=colors[fill_colors[i]], label=int(r), alpha=0.5)

leg = ax.legend(loc='upper left', title = 'repressors', fontsize=9) 
leg.get_title().set_fontsize(9)

if (show_data == True):
    data_name = 'data'
else:
    data_name = 'nodata'
if (show_fit == True):
    fit_name = 'fit'
else:
    fit_name = 'nofit'
if (show_predictions == True):
    pred_name = 'pred'
else:
    pred_name = 'nopred'
plt.savefig(f'../figs/O2_{data_name}_{fit_name}_{pred_name}.pdf', 
            bbox_inches='tight')
# %%


# %%
