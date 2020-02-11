# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.bayes
colors = phd.viz.phd_style()

# Load the data and sampling statistics
data = pd.read_csv('../data/pathological_F_data.csv')
stats = pd.read_csv('../data/pathological_F_stats.csv')
stats['draw'] = stats.groupby(['true_bohr']).ngroup()

TRUE_DATA = 1
SIM_DATA = 1
INF_DATA = 1
DEL_DATA = 1
DIFF_COLOR = 1

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 1, figsize=(3, 4))
ax[0].set_xlabel('free energy [$k_BT$]')
ax[0].set_ylabel('fold-change')
ax[1].set_xlabel('true free energy [$k_BT$]')
ax[1].set_ylabel('inferred free energy [$k_BT$]')
ax[0].set_ylim([-0.25, 1.45])

# ##############################################################################
# AGREEMENT AND MASTER CURVE
# ##############################################################################
bohr_range = np.linspace(-10, 10, 200)
ax[0].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', 
                        label='__nolegend__')
ax[1].plot(bohr_range, bohr_range, 'k-', label='perfect agreement')

# ##############################################################################
# SIMULATED DATA
# ##############################################################################
for g, d in data.groupby(['draw']):
    if g%10 == 0:
        if TRUE_DATA == 1:
            if g == 0:
                label = 'true fold-change'
            else:
                label = '__nolegend__'
            ax[0].plot(d['bohr'], d['fc_mu'], '.', ms=5, 
                         color=colors['black'], label=label, zorder=100,
                         markerfacecolor=colors['gray'], markeredgewidth=0.75)
        if SIM_DATA == 1:
            if g == 0:
                label = 'simulated data'
            else:
                label = '__nolegend__'
            ax[0].plot(d['bohr'], d['fold_change'], '.', ms=3, 
                        color=colors['purple'], zorder=99,
                        label=label)

# ##############################################################################
# INFERRED DATA
# ##############################################################################
labels = []
for g, d in stats.groupby('draw'):
    _fc = d[d['parameter']=='fc_mu']
    _sig = d[d['parameter']=='fc_sigma']
    _F = d[d['parameter']=='empirical_bohr']
    if g%10 == 0:
        if g == 0:
            legend = 'inferred fold-change'
        else:
            legend = '__nolegend__'
        if INF_DATA == 1:
            ax[0].plot(d['true_bohr'].values[0], _fc['median'], '.', ms=5, 
                        color=colors['orange'], label=legend, zorder=1001,
                        markerfacecolor=colors['light_orange'],
                        markeredgewidth=0.75)
            ax[0].vlines(d['true_bohr'].values[0], _fc['hpd_min'], _fc['hpd_max'], lw=0.75, color=colors['orange'], label='__nolegend__')



        if DEL_DATA==1:
            if DIFF_COLOR == 1:
                if (_fc['median'].values[0] < _sig['median'].values[0]):
                    _color = colors['blue']
                    face = colors['light_blue']
                    label = 'µ < σ'
                elif (1 - _fc['median'].values[0] < _sig['median'].values[0]):
                    _color = colors['green']
                    face = colors['light_green']
                    label = '1 - µ < σ'
                else:
                    _color = colors['orange']
                    face = colors['light_orange']
                    label = 'µ > σ; 1 - μ > σ'
            else:
                face = colors['light_orange']
                _color = colors['orange']
                label = 'inferred free energy'
            if label not in labels:
                labels.append(label)
            else:
                label = '__nolegend__'

            ax[1].plot(d['true_bohr'].values[0], _F['median'], '.', 
                        color=_color, label=label, ms=5, markerfacecolor=face, markeredgewidth=0.75)
            ax[1].vlines(d['true_bohr'].values[0], _F['hpd_min'], 
                        _F['hpd_max'], lw=0.75, color=_color,     label='__nolegend__') 


# ##############################################################################
# LEGEND AND SAVING
# ##############################################################################
ax[0].legend(loc='upper left')
ax[1].legend(loc='upper left')
plt.tight_layout()

if (TRUE_DATA == 1) & (SIM_DATA == 0) &\
    (INF_DATA == 0) & (DEL_DATA == 0) & (DIFF_COLOR==0) :
    plt.savefig('../figs/pathological_F_truth.pdf')
elif (TRUE_DATA == 1) & (SIM_DATA == 1) &\
    (INF_DATA == 0) & (DEL_DATA == 0) & (DIFF_COLOR==0) :
    plt.savefig('../figs/pathological_F_sim.pdf')
elif (TRUE_DATA == 1) & (SIM_DATA == 1) &\
    (INF_DATA == 1) & (DEL_DATA == 0) & (DIFF_COLOR==0) :
    plt.savefig('../figs/pathological_F_inf.pdf')
elif (TRUE_DATA == 1) & (SIM_DATA == 1) &\
    (INF_DATA == 1) & (DEL_DATA == 1) & (DIFF_COLOR==0) :
    plt.savefig('../figs/pathological_F_disagreement.pdf')
elif (TRUE_DATA == 1) & (SIM_DATA == 1) &\
    (INF_DATA == 1) & (DEL_DATA == 1) & (DIFF_COLOR==1) :
    plt.savefig('../figs/pathological_F_reason.pdf')


