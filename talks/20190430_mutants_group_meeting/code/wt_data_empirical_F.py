# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
import phd.stats
colors = phd.viz.phd_style()
constants = phd.thermo.load_constants()

# Load all of the data. 
data = pd.read_csv('../data/RazoMejia2018_2019_full_data.csv')
data = data[data['method']=='flow cytometry']
ops = [constants[o] for o in data['operator'].values]
data['ref_bohr'] = phd.thermo.SimpleRepression(R=data['repressors'], ep_r=ops,
                                            ka=constants['Ka'], 
                                            ki=constants['Ki'],
                                            effector_conc=data['IPTGuM'],
                                            ep_ai=constants['ep_AI']).bohr_parameter()
stats = pd.read_csv('../data/wt_empirical_F_statistics.csv')
stats = stats[stats['method']=='flow cytometry']


# Define the bohr range and plotting arguments. 
bohr_range = np.linspace(-10, 10, 200)

# Define colors as necessary
meth_colors = {'microscopy':colors['purple'], 'flow cytometry':colors['orange']}
meth_face = {'microscopy':colors['light_purple'], 
            'flow cytometry':colors['light_orange']}

INF = 1
# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(1, 2, figsize=(5, 2.5))
ax[0].set_xlabel('free energy [$k_BT$]')
ax[0].set_ylabel('fold-change')

ax[1].set_xlabel('predicted free energy [$k_BT$]')
ax[1].set_ylabel('inferred free energy [$k_BT$]')

# ##############################################################################
# COLLAPSE CURVE AND PERFECT AGREEMENT
# ##############################################################################
ax[0].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', lw=1,
           label='__nolegend__')
ax[1].plot(bohr_range, bohr_range, 'k-', lw=1,
           label='perfect agreement')


# ##############################################################################
#  DATA
# ##############################################################################
i = 0
for g, d in data.groupby(['method']):
    if i == 0:
        label = 'individual replicate'
        i += 1
    else:
        label = '__nolegend__'

    ax[0].plot(d['ref_bohr'], d['fold_change'], '.', alpha=0.5, ms=1, 
                color=meth_face[g], label=label)

# ##############################################################################
#  INFERRED DATA
# ##############################################################################
labels = []
for g, d in stats.groupby(['repressors', 'operator', 'IPTGuM']): 
    legend = 'inferred mean fold-change' 
    if legend not in labels:
        labels.append(legend)
    else:
        legend = '__nolegend__'
    mu = d[d['parameter']=='fc_mu'] 
    sig = d[d['parameter']=='fc_sigma']
    bohr = d[d['parameter']=='empirical_bohr']

    if (mu['median'].values[0] < sig['median'].values[0]) | (1 - mu['median'].values[0] < sig['median'].values[0]):
        _c = colors['purple']
        face = 'w'
        label = 'µ < σ or 1 - μ < σ '
    else:
        _c = colors['orange']
        face = 'w'
        label= 'µ > σ and 1 - µ > σ'
    if label not in labels:
        labels.append(label)
    else:
        label = '__nolegend__'
    if INF == 1: 
        ax[0].plot(mu['pred_bohr'], mu['median'], '.', ms=4, color=_c,
                markerfacecolor=face, markeredgewidth=0.75, label=label)
        ax[0].vlines(mu['pred_bohr'], mu['hpd_min'], mu['hpd_max'], lw=0.75,
        color=_c, label='__nolegend__')
        ax[1].plot(bohr['pred_bohr'], bohr['median'],
                    '.', ms=4, color=_c, markerfacecolor=face,
                    markeredgewidth=0.75, label=label)
        ax[1].vlines(bohr['pred_bohr'], bohr['hpd_min'], bohr['hpd_max'], lw=0.75,
                    color=_c, label = '__nolegend__')

# ##############################################################################
# LEGENDS AND SAVING
# ##############################################################################
ax[0].legend( loc='upper left')
ax[1].legend(loc='upper left')

plt.tight_layout()
if INF == 1:
    plt.savefig('../figs/wt_empirical_data.pdf', bbox_inches='tight')
else:
    plt.savefig('../figs/wt_empirical_inf.pdf', bbox_inches='tight')

