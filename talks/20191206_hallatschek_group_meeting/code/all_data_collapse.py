# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, color_list = phd.viz.phd_style()
constants = phd.thermo.load_constants()
constants['Oid'] = -17
mut_colors = phd.viz.color_selector('mut')

# Load the data from the mutants work.
data = pd.read_csv('../data/summarized_data.csv')
data = data[data['mutant']!='wt']
epRA_stats = pd.read_csv('../data/DNA_binding_energy_summary.csv')
epRA_stats = epRA_stats[epRA_stats['repressors']==260]
allo_stats = pd.read_csv('../data/KaKi_epAI_summary.csv')
allo_stats = allo_stats[allo_stats['operator']=='O2']

# Load the data from the old gods
old_gods = pd.read_csv('../data/Garcia2011_Brewster2014_data.csv')
new_gods = pd.read_csv('../data/RazoMejia2018_data.csv')
new_gods['repressors'] *= 2
new_gods = new_gods[new_gods['repressors'] > 0]
new_gods.rename(columns={'IPTG_uM':'IPTGuM', 'fold_change_A':'fold_change'},
    inplace=True)
new_gods = new_gods.groupby(['repressors', 'operator', 'IPTGuM']).agg(('mean', 'sem')).reset_index()

# Define plotting constants
OLD_GODS = 0 
NEW_GODS = 0
DNA_MUTS = 0
IND_MUTS = 0
DBL_MUTS = 0
bohr_range = np.linspace(-10, 10, 200)

# Define colors and glyphs
glyphs = {'garcia':'X', 'brewster':'s'}
_color = {'garcia': colors['dark_red'], 'brewster':colors['dark_blue']}
legend = {'garcia': 'Garcia & Phillips 2011', 'brewster':'Brewster et al. 2014'}
# Define the repressor colors. 
edge_reps = {22:colors['dark_purple'],  60:colors['dark_orange'],  
            124: colors['black'], 260: colors['dark_red'], 
            1220: colors['dark_blue'], 1740: colors['dark_green']}
fill_reps = {22:colors['light_purple'],  60:colors['light_orange'],   
            124: colors['light_grey'], 260: colors['light_red'], 
            1220: colors['light_blue'], 1740: colors['light_green']}
#
op_glyphs = {'O1':'s', 'O2':'o', 'O3':'^'}

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(1, 1, figsize=(4, 2.5))
ax.set_ylim([-0.05, 1.1])
ax.set_xlabel('free energy [$k_BT$]')
ax.set_ylabel('fold-change')

# ##############################################################################
# COLLAPSE CURVE
# ##########################3333333#3333########################################
ax.plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', lw=1)

# ##############################################################################
# OLD_GODS
# ##############################################################################
if OLD_GODS == 1:
    for g, d in old_gods.groupby('author'):
        ops = [constants[o] for o in d['operator'].values]
        bohr = phd.thermo.SimpleRepression(R=d['repressor'], ep_r=ops,
                                           ka=constants['Ka'], ki=constants['Ki'],
                     ep_ai=constants['ep_AI'], effector_conc=0).bohr_parameter()
        ax.plot(bohr, d['fold_change'], marker=glyphs[g], color=_color[g], linestyle='none',
        label=legend[g], markeredgewidth=0.75, alpha=0.75,
        ms=3)
# ##############################################################################
# NEW GODS
# ##############################################################################
if NEW_GODS == 1:
    for g, d in new_gods.groupby(['repressors', 'operator']):
        ops = [constants[o] for o in new_gods['operator'].values]
        bohr = phd.thermo.SimpleRepression(R=d['repressors'], ep_r=constants[g[1]],
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
            effector_conc=d['IPTGuM']).bohr_parameter() 
        ax.errorbar(bohr, d['fold_change']['mean'], 
                    d['fold_change']['sem'], fmt=op_glyphs[g[1]],  color=edge_reps[g[0]],
                     markeredgewidth=0.75, markerfacecolor=fill_reps[g[0]], 
                     linestyle='none', lw=0.75, capsize=1, ms=4, alpha=0.75)

# ##############################################################################
# DNA MUTS
# ##############################################################################
op_glyphs = {'O1':'^', 'O2':'v', 'O3':'D'}
if DNA_MUTS == 1:
    for g, d in data[data['class']=='DNA'].groupby(['mutant']):
        if g == 'Q21M':
            label = 'Q18M'
        elif g == 'Q21A':
            label = 'Q18A'
        elif g == 'Y20I':
            label = 'Y17I'
        else:
            label = g
        ep_RA = epRA_stats[(epRA_stats['mutant']==g) &
                (epRA_stats['parameter']=='ep_RA')]['median'].values[0]
        bohr = phd.thermo.SimpleRepression(R=d['repressors'], ep_r=ep_RA, 
                                           ka=constants['Ka'], ki=constants['Ki'],
                                           ep_ai=constants['ep_AI'],
                                           effector_conc=d['IPTGuM']).bohr_parameter()
        ax.errorbar(bohr, d['mean'], d['sem'], fmt='^', color=mut_colors[g],
                label=label, ms=5, markeredgewidth=0.75, markerfacecolor='w')

if IND_MUTS == 1:
    for g, d in data[data['class']=='IND'].groupby(['mutant']):
        if g == 'F164T':
            label = 'F161T'
        elif g == 'Q294K':
            label = 'Q291K'
        elif g == 'Q294V':
            label = 'Q291V'
        elif g == 'Q294R':
            label = 'Q291R'
        else: 
            label = g
        _stats = allo_stats[(allo_stats['mutant']==g)]
        ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
        ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
        ep_AI = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
        ops = [constants[o] for o in d['operator'].values]
        bohr = phd.thermo.SimpleRepression(R=d['repressors'], ep_r=ops, 
                                           ka=ka, ki=ki,
                                           ep_ai=ep_AI,
                                           effector_conc=d['IPTGuM']).bohr_parameter()
        ax.errorbar(bohr, d['mean'], d['sem'], fmt='p', color=mut_colors[g],
                label=label, ms=5, markeredgewidth=1,  markerfacecolor='w')

if DBL_MUTS == 1:
    for g, d in data[data['class']=='DBL'].groupby(['mutant']):
        dna, ind = g.split('-')
        if dna == 'Q21M':
            label = 'Q18M'
        elif dna == 'Q21A':
            label = 'Q18A'
        elif dna == 'Y20I':
            label = 'Y17I'
        else:
            label = g
        if ind == 'F164T':
            label += '-F161T'
        elif ind == 'Q294K':
            label += '-Q291K'
        elif ind == 'Q294V':
            label += '-Q291V'
        elif g == 'Q294R':
            label += '-Q291R'
        else: 
            label = g 
        
        ep_RA = epRA_stats[(epRA_stats['mutant']==g.split('-')[0]) & 
                           (epRA_stats['parameter']=='ep_RA')]['median'].values[0]
        _stats = allo_stats[(allo_stats['mutant']==g.split('-')[1])]
        ka = _stats[_stats['parameter']=='Ka']['median'].values[0]
        ki = _stats[_stats['parameter']=='Ki']['median'].values[0]
        ep_AI = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
        bohr = phd.thermo.SimpleRepression(R=d['repressors'], ep_r=ep_RA, 
                                           ka=ka, ki=ki,
                                           ep_ai=ep_AI,
                                           effector_conc=d['IPTGuM']).bohr_parameter()
        ax.errorbar(bohr, d['mean'], d['sem'], fmt='*', color=mut_colors[g],
                label=label, ms=6, markeredgewidth=0.75,  markerfacecolor='w')

# ax.legend(loc='upper left')

plt.tight_layout()
if (OLD_GODS==0) & (NEW_GODS==0) & (DNA_MUTS==0) & (IND_MUTS==0) & (DBL_MUTS==0):
    plt.savefig('../figs/collapse_only.pdf', bbox_inches='tight', facecolor='white')
if (OLD_GODS==1) & (NEW_GODS==0) & (DNA_MUTS==0) & (IND_MUTS==0) & (DBL_MUTS==0):
    plt.savefig('../figs/old_gods_collapse.pdf', bbox_inches='tight', facecolor='white')
if (OLD_GODS==0) & (NEW_GODS==1) & (DNA_MUTS==0) & (IND_MUTS==0) & (DBL_MUTS==0):
    plt.savefig('../figs/new_gods_collapse.pdf', bbox_inches='tight', facecolor='white')
if (OLD_GODS==1) & (NEW_GODS==1) & (DNA_MUTS==1) & (IND_MUTS==0) & (DBL_MUTS==0):
    plt.savefig('../figs/new_gods_DNA_collapse.pdf', bbox_inches='tight', facecolor='white')
if (OLD_GODS==1) & (NEW_GODS==1) & (DNA_MUTS==1) & (IND_MUTS==1) & (DBL_MUTS==0):
    plt.savefig('../figs/new_gods_DNA_IND_collapse.pdf', bbox_inches='tight', facecolor='white')
if (OLD_GODS==1) & (NEW_GODS==1) & (DNA_MUTS==1) & (IND_MUTS==1) & (DBL_MUTS==1):
    plt.savefig('../figs/new_gods_DNA_IND_DBL_collapse.pdf', bbox_inches='tight', facecolor='white')



       




#%%
