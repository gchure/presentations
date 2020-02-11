# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz 
import phd.stats
import phd.thermo
colors = phd.viz.phd_style()

# Load some of the sampling information
epRA_samples = pd.read_csv('../data/DNA_binding_energy_samples.csv')
allo_samples = pd.read_csv('../data/KaKi_epAI_samples.csv')

# Restrict to a single mutant
epRA_samples = epRA_samples[(epRA_samples['mutant']=='Q21A') & 
                            (epRA_samples['repressors']==260)]
allo_samples = allo_samples[(allo_samples['mutant']=='Q294K') &
                            (allo_samples['operator']=='O2')]



# Draw the random samples 
ep_draws = np.random.choice(epRA_samples['ep_RA'].values, size=int(1E4), replace=True)
allo_draws = allo_samples.sample(n=int(1E4), replace=True)
ka_draws = allo_draws['Ka']
ki_draws = allo_draws['Ki']
epAI_draws = allo_draws['ep_AI']

# Select some of the draws
ep_sel = ep_draws[:10]
ka_sel = ka_draws[:10]
ki_sel = ki_draws[:10]
epAI_sel = epAI_draws[:10]

# Define the range of concentrations
c_range = np.logspace(-2, 4, 200)
c, ka, ki ,ep, ep_AI = np.meshgrid(c_range, ka_sel, ki_sel, ep_sel, epAI_sel)
arch = phd.thermo.SimpleRepression(R=260, ep_r=ep, ka=ka, ki=ki, ep_ai=ep_AI,
                                   effector_conc=c).fold_change()
# ##############################################################################
# FIGURE INSTANTIATION -- DRAWS 
# ############################################################################## 
SEL_DRAWS = 0
ALL_DRAWS = 1
fig, ax = plt.subplots(2, 1, figsize=(2, 3))
ax[0].set_xlabel('DNA binding energy [$k_BT$]')
ax[0].set_yticks([])
ax[0].set_ylim([0, 6.5])
ax[0].spines['left'].set_visible(False)
ax[1].set_xlabel('$K_A$ [µM]')
ax[1].set_ylabel('$K_I$ [µM]')

#  Hist the DNA binding samples
bins = 20
ax[0].hist(epRA_samples['ep_RA'], bins=bins, color=colors['light_purple'], 
           edgecolor=colors['dark_purple'], density=True)

# Make a hexbin of Ka and Ki
ax[1].hexbin(allo_samples['Ka'],allo_samples['Ki'],xscale='log', yscale='log', 
             gridsize=50, cmap='magma')

if SEL_DRAWS == 1:
    ep = ep_sel
    ka = ka_sel
    ki = ki_sel
    ep_AI = epAI_sel
elif ALL_DRAWS == 1:
    ep = ep_draws
    ka = allo_draws['Ka']
    ki = allo_draws['Ki']
    ep_AI = allo_draws['ep_AI']

if (SEL_DRAWS == 1) | (ALL_DRAWS == 1):
    # Plot the epRA draws
    ax[0].plot(ep, 6 * np.ones(len(ep)) * (1 + np.random.normal(0, 0.005, len(ep))), '.',
               color='k', ms=4, markerfacecolor='w', markeredgewidth=0.75,
               alpha=0.75, rasterized=True)
    # Plot the KaKi_draws
    ax[1].plot(ka, ki, '.', color='w', ms=5, markeredgecolor='k', 
              markeredgewidth=0.75, alpha=0.75, rasterized=True)
plt.tight_layout()
if SEL_DRAWS == 1:
    plt.savefig('../figs/dbl_sel.pdf', bbox_inches='tight')
elif (SEL_DRAWS ==0) & (ALL_DRAWS==0):
    plt.savefig('../figs/dbl_dists.pdf', bbox_inches='tight')
else:
    plt.savefig('../figs/dbl_samples.pdf', bbox_inches='tight')


# ##############################################################################
# CURVE FIG
# ##############################################################################
SEL_DRAWS = 0
ALL_DRAWS = 1
HPD = 0
fig, ax = plt.subplots(1,1, figsize=(3.5, 3.5))
ax.set_xscale('log')
ax.set_ylim([0.8, 1.1])
ax.set_xlabel('IPTG [µM]')
ax.set_ylabel('fold-change')


if SEL_DRAWS == 1:
    # Select some of the draws
    ep_sel = ep_draws[:10]
    ka_sel = ka_draws[:10]
    ki_sel = ki_draws[:10]
    epAI_sel = epAI_draws[:10]

    # Define the range of concentrations
    c_range = np.logspace(-2, 4, 200)
    c, ka, ki ,ep, ep_AI = np.meshgrid(c_range, ka_sel, ki_sel, ep_sel, epAI_sel)
    arch = phd.thermo.SimpleRepression(R=260, ep_r=ep, ka=ka, ki=ki, ep_ai=ep_AI,
                                   effector_conc=c).fold_change()
    for i in range(len(ep_sel)):
        ax.plot(c_range, arch[i, :, i, i, i], '-', lw=0.1, 
                color=colors['dark_purple'])

if ALL_DRAWS == 1:
    # Define the range of concentrations
    for i in range(int(len(ep_draws)/10)):
        arch = phd.thermo.SimpleRepression(R=260, ep_r=ep_draws[i], 
                ka=ka_draws.values[i], ki=ki_draws.values[i], 
                ep_ai=epAI_draws.values[i],
                effector_conc=c_range).fold_change()
        ax.plot(c_range, arch, '-', lw=0.01, 
                color=colors['dark_purple'], alpha=0.8)


if HPD==1:
    cred_region = np.zeros((2, len(c_range)))
    for i, c in enumerate(c_range):
        arch = phd.thermo.SimpleRepression(R=260, ep_r=ep_draws, 
                ka=ka_draws.values, ki=ki_draws.values, 
                ep_ai=epAI_draws.values,
                effector_conc=c).fold_change()
        cred_region[:, i] = phd.stats.compute_hpd(arch, 0.95)
    ax.fill_between(c_range, cred_region[0, :], cred_region[1, :],  color=colors['purple'],
                    alpha=0.5)
plt.tight_layout()
if (SEL_DRAWS == 1):
    plt.savefig('../figs/DBL_prediction_sel.pdf', bbox_inches='tight')
elif (ALL_DRAWS==1):
    plt.savefig('../figs/DBL_prediction_all.pdf', bbox_inches='tight')
elif HPD==1:
    plt.savefig('../figs/DBL_prediction_HPD.pdf', bbox_inches='tight')


